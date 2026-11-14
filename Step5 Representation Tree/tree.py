#!/usr/bin/env python3
"""Cluster ProstT5 embeddings and plot the representation tree and UMAP views."""

import argparse
from pathlib import Path
import pickle
import sys

import matplotlib

matplotlib.use('Agg')  # Save figures without a display or notebook kernel.
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import pdist

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TARGET_COLOR = '#00DFA2'
REFERENCE_COLOR = '#FA6868'
CANDIDATE_COLOR = '#4D96FF'
CLUSTER_COLOR = '#FFCC00'


def load_embeddings(path):
    """Load a trusted Step 4 pickle, keeping its row-to-sequence mapping intact."""
    with path.open('rb') as handle:
        data = pickle.load(handle)
    keys = {'labels', 'sequences', 'representations'}
    if not isinstance(data, dict) or not keys.issubset(data):
        raise ValueError(f'{path}: expected a dictionary containing {sorted(keys)}.')
    matrix = np.asarray(data['representations'])
    if matrix.ndim != 2 or not all(matrix.shape) or matrix.dtype.kind not in 'fiu':
        raise ValueError(f'{path}: representations must be a nonempty numeric 2-D array.')
    if not np.isfinite(matrix).all():
        raise ValueError(f'{path}: representations must contain only finite values.')
    labels, sequences = data['labels'], data['sequences']
    if not isinstance(labels, (list, tuple, np.ndarray)) or not isinstance(sequences, (list, tuple, np.ndarray)):
        raise ValueError(f'{path}: labels and sequences must be lists or arrays.')
    if len(labels) != len(matrix) or len(sequences) != len(matrix):
        raise ValueError(f'{path}: labels, sequences, and representations need matching lengths.')
    if any(not isinstance(label, str) or not label.strip() or
           any(char in label for char in '\t\n\r') for label in labels):
        raise ValueError(f'{path}: labels must be nonempty strings without tabs or newlines.')
    labels = [label.strip() for label in labels]  # Historical IDs have trailing spaces.
    if len(set(labels)) != len(labels):
        raise ValueError(f'{path}: duplicate labels after stripping surrounding whitespace.')
    if any(not isinstance(seq, str) or not seq or any(char.isspace() for char in seq)
           for seq in sequences):
        raise ValueError(f'{path}: sequences must be nonempty strings without whitespace.')
    return labels, list(sequences), matrix


def read_candidate_ids(path):
    """Read one accession per line; empty lines and # comments are ignored."""
    if path is None:
        return set()
    return {line.strip() for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith('#')}


def cluster_embeddings(matrix, threshold, coarse_threshold):
    """Match the notebook: Euclidean distances and Ward linkage, without scaling."""
    pairs = len(matrix) * (len(matrix) - 1) // 2
    print(f'Clustering {matrix.shape[0]:,} proteins ({matrix.shape[1]} dimensions).', flush=True)
    print(f'Pairwise distances alone require approximately {pairs * 8 / 1024**3:.2f} GiB; '
          'linkage needs additional memory.', flush=True)
    distances = pdist(matrix, metric='euclidean')
    tree = linkage(distances, method='ward')
    return (tree, fcluster(tree, threshold, criterion='distance'),
            fcluster(tree, coarse_threshold, criterion='distance'))


def fit_umap(matrix, args):
    # Import only when needed; this is the umap-learn package, not the umap package.
    from umap import UMAP

    neighbors = min(args.n_neighbors, len(matrix) - 1)
    if neighbors != args.n_neighbors:
        print(f'Using n_neighbors={neighbors} for {len(matrix)} input proteins.', flush=True)
    print('Computing UMAP...', flush=True)
    return UMAP(n_neighbors=neighbors, min_dist=args.min_dist, n_components=2,
                metric='euclidean', spread=args.spread, random_state=args.seed,
                n_jobs=1, init='random').fit_transform(matrix)


def cluster_colors(clusters, target_cluster):
    colors = plt.get_cmap('Pastel1')((np.asarray(clusters) - 1) % 7)
    colors[np.asarray(clusters) == target_cluster] = matplotlib.colors.to_rgba(CLUSTER_COLOR)
    return colors


def save_figure(fig, output, name, dpi, pdf=False):
    fig.savefig(output / f'{name}.png', dpi=dpi, bbox_inches='tight')
    if pdf:
        fig.savefig(output / f'{name}.pdf', bbox_inches='tight')
    plt.close(fig)


def plot_tree(tree, table, target_index, output, dpi):
    """Draw the notebook's 300-degree radial tree, cluster ring, and annotations."""
    layout = dendrogram(tree, no_plot=True)
    leaves = np.asarray(layout['leaves'])
    x = np.asarray(layout['icoord'])
    distances = np.asarray(layout['dcoord'])
    span, start = 5 * np.pi / 3, -np.pi / 3
    leaf_width = span / len(leaves)
    # Put leaves at sector centers so the ring keeps its gap even for tiny inputs.
    theta = (x - x.min()) / (x.max() - x.min()) * (span - leaf_width) + start + leaf_width / 2
    max_distance = distances.max()
    radii = 1 - distances / (max_distance if max_distance > 0 else 1)
    # Positive radii avoid backend-dependent clipping of negative polar coordinates.
    segments = []
    for angles, radius in zip(theta, radii):
        smooth_angles = np.r_[angles[0], np.linspace(angles[1], angles[2], 100), angles[3]]
        smooth_radii = np.r_[radius[0], np.linspace(radius[1], radius[2], 100), radius[3]]
        segments.append(np.column_stack((smooth_angles, smooth_radii)))
    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'projection': 'polar'})
    ax.add_collection(LineCollection(segments, colors='tab:grey', linewidths=0.5,
                                     transform=ax.transData))
    leaf_angles = np.linspace(start + leaf_width / 2, start + span - leaf_width / 2, len(leaves))
    clusters = table['hierarchy_cluster'].to_numpy()[leaves]
    target_cluster = table.loc[target_index, 'hierarchy_cluster']
    colors = cluster_colors(clusters, target_cluster)
    # One ring sector per consecutive cluster, rather than one artist per protein.
    boundaries = np.r_[0, np.flatnonzero(np.diff(clusters)) + 1, len(leaves)]
    for left, right in zip(boundaries[:-1], boundaries[1:]):
        ax.bar((leaf_angles[left] + leaf_angles[right - 1]) / 2, 0.05,
               width=leaf_width * (right - left), bottom=1, color=colors[left],
               edgecolor='none', align='center')
    references = table['is_reference'].to_numpy()[leaves]
    candidates = table['is_candidate'].to_numpy()[leaves]
    target = leaves == target_index
    for mask, radius, color in [(references & ~target, 1.13, REFERENCE_COLOR),
                                (candidates, 1.07, CANDIDATE_COLOR),
                                (target, 1.13, TARGET_COLOR)]:
        ax.scatter(leaf_angles[mask], np.full(mask.sum(), radius), s=100, marker='^',
                   color=color, edgecolor='#6b6b6b', linewidth=0.5, zorder=3)
    ax.set_ylim(0, 1.18)
    ax.set_axis_off()
    save_figure(fig, output, 'dendrogram', dpi, pdf=True)
    table.iloc[leaves][['label', 'hierarchy_cluster']].to_csv(output / 'ordered_labels.csv', index=False)


def detail_bounds(points):
    lower, upper = points.min(axis=0), points.max(axis=0)
    width = max(float(np.max(upper - lower)), 1.0) * 1.2
    center = (lower + upper) / 2
    return center - width / 2, center + width / 2


def plot_umap(xy, table, target_index, output, dpi):
    clusters = table['hierarchy_cluster'].to_numpy()
    target_cluster = clusters[target_index]
    in_target = clusters == target_cluster
    reference = table['is_reference'].to_numpy().copy()
    reference[target_index] = False
    candidates = table['is_candidate'].to_numpy()
    colors = cluster_colors(clusters, target_cluster)
    lower, upper = detail_bounds(xy[in_target])
    for detail in (False, True):
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.scatter(xy[:, 0], xy[:, 1], c=colors, s=8 if detail else 3,
                   alpha=0.5, edgecolors='none', rasterized=True)
        ax.scatter(xy[in_target, 0], xy[in_target, 1], color=CLUSTER_COLOR,
                   label='Target cluster', s=20 if detail else 5, edgecolors='none')
        size = 80 if detail else 35
        for mask, color, label in [(reference, REFERENCE_COLOR, 'Other references'),
                                    (candidates, CANDIDATE_COLOR, 'Candidates')]:
            if mask.any():
                ax.scatter(xy[mask, 0], xy[mask, 1], color=color, label=label,
                           marker='^', s=size, edgecolor='#6b6b6b', linewidth=0.5)
        ax.scatter(*xy[target_index], color=TARGET_COLOR, label=table.loc[target_index, 'label'],
                   marker='^', s=size, edgecolor='#6b6b6b', linewidth=0.5, zorder=4)
        if detail:
            ax.set_xlim(lower[0], upper[0])
            ax.set_ylim(lower[1], upper[1])
            ax.legend(loc='best', fontsize=8)
        else:
            ax.add_patch(Rectangle(lower, *(upper - lower), fill=False,
                                   edgecolor='#6b6b6b', linewidth=1))
        ax.set_xlabel('UMAP 1')
        ax.set_ylabel('UMAP 2')
        ax.set_aspect('equal', adjustable='box')
        fig.tight_layout()
        save_figure(fig, output, 'umap_detail' if detail else 'umap', dpi)


def export_tables(table, target_index, output):
    table.drop(columns='sequence').to_csv(output / 'sequence_clusters.tsv', sep='\t', index=False)
    target_cluster = table.loc[target_index, 'hierarchy_cluster']
    selected = table[(table.hierarchy_cluster == target_cluster) & ~table.is_reference]
    with (output / 'target_cluster.fasta').open('w') as handle:
        for row in selected.itertuples(index=False):
            handle.write(f'>{row.label}\n')
            for start in range(0, len(row.sequence), 80):
                handle.write(row.sequence[start:start + 80] + '\n')
    print(f'Target reference {table.loc[target_index, "label"]}: cluster {target_cluster}, '
          f'{len(selected):,} discovery sequences exported.', flush=True)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--references', type=Path, default=ROOT / 'Step4 ProstT5 Representation/cas12_refseq.pkl',
                        help='Trusted reference embedding pickle')
    parser.add_argument('--discovered', type=Path, default=ROOT / 'Step4 ProstT5 Representation/AsCas12f_NR_s50.pkl',
                        help='Trusted discovery embedding pickle')
    parser.add_argument('--output-dir', type=Path, default=HERE / 'outputs', help='Directory for generated results')
    parser.add_argument('--target-reference', default='8J12-Cas12f', help='Reference ID defining the highlighted cluster')
    candidates = parser.add_mutually_exclusive_group()
    candidates.add_argument('--candidate-ids', type=Path, help='Candidate IDs, one per line; default: Step 6 panel if present')
    candidates.add_argument('--no-candidates', action='store_true', help='Disable candidate annotations')
    parser.add_argument('--threshold', type=float, default=10, help='Ward distance cut for plots and target selection')
    parser.add_argument('--coarse-threshold', type=float, default=20, help='Second, coarser Ward distance cut for the table')
    parser.add_argument('--n-neighbors', type=int, default=100, help='UMAP neighborhood size (capped at N - 1)')
    parser.add_argument('--min-dist', type=float, default=0.1, help='UMAP minimum distance')
    parser.add_argument('--spread', type=float, default=3, help='UMAP spread')
    parser.add_argument('--seed', type=int, default=42, help='UMAP random seed')
    parser.add_argument('--dpi', type=int, default=300, help='PNG resolution')
    return parser


def run(args):
    for path in (args.references, args.discovered):
        if not path.is_file():
            raise FileNotFoundError(f'Embedding file not found: {path}. '
                                    'Generate it in Step 4 or pass --references / --discovered.')
    ref_labels, ref_sequences, ref_matrix = load_embeddings(args.references)
    found_labels, found_sequences, found_matrix = load_embeddings(args.discovered)
    if ref_matrix.shape[1] != found_matrix.shape[1]:
        raise ValueError('Reference and discovery embedding dimensions must match.')
    labels = ref_labels + found_labels
    if len(set(labels)) != len(labels):
        raise ValueError('Reference and discovery labels must be unique across both inputs.')
    if len(labels) < 3:
        raise ValueError('At least three proteins in total are required for UMAP.')
    target_label = args.target_reference.strip()
    if target_label not in ref_labels:
        raise ValueError(f'Target reference {target_label!r} is not in the reference pickle.')
    target_index = ref_labels.index(target_label)
    candidate_path = args.candidate_ids
    if candidate_path is None:
        default_candidates = ROOT / 'Step6 Selected Candidates/candidate_ids.txt'
        candidate_path = default_candidates if default_candidates.is_file() else None
    candidate_ids = set() if args.no_candidates else read_candidate_ids(candidate_path)
    missing = candidate_ids - set(found_labels)
    if missing:
        print('Warning: candidate IDs absent from discovery embeddings (not annotated): '
              + ', '.join(sorted(missing)), file=sys.stderr)
    matrix = np.concatenate((ref_matrix, found_matrix))
    table = pd.DataFrame({'label': labels, 'sequence': ref_sequences + found_sequences,
                          'is_reference': [True] * len(ref_labels) + [False] * len(found_labels)})
    table['length_aa'] = table.sequence.str.len()
    table['is_candidate'] = table.label.isin(candidate_ids) & ~table.is_reference
    tree, clusters, coarse_clusters = cluster_embeddings(matrix, args.threshold, args.coarse_threshold)
    table['hierarchy_cluster'] = clusters
    table['hierarchy_cluster_coarse'] = coarse_clusters
    print(f'Clusters at distance {args.threshold:g}: {len(np.unique(clusters))}; '
          f'at {args.coarse_threshold:g}: {len(np.unique(coarse_clusters))}.', flush=True)
    xy = fit_umap(matrix, args)
    table[['UMAP_1', 'UMAP_2']] = xy
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.size': 12, 'font.family': 'DejaVu Sans'})
    print('Saving figures and tables...', flush=True)
    plot_tree(tree, table, target_index, args.output_dir, args.dpi)
    plot_umap(xy, table, target_index, args.output_dir, args.dpi)
    export_tables(table, target_index, args.output_dir)
    print(f'Results saved to {args.output_dir.resolve()}', flush=True)


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if not np.isfinite(args.threshold) or args.threshold <= 0:
        parser.error('--threshold must be finite and positive.')
    if not np.isfinite(args.coarse_threshold) or args.coarse_threshold < args.threshold:
        parser.error('--coarse-threshold must be finite and at least --threshold.')
    if args.n_neighbors < 2:
        parser.error('--n-neighbors must be at least 2.')
    if not np.isfinite(args.spread) or args.spread <= 0 or not np.isfinite(args.min_dist) or not 0 <= args.min_dist <= args.spread:
        parser.error('UMAP requires finite values with 0 <= --min-dist <= --spread and --spread > 0.')
    if not 0 <= args.seed < 2**32:
        parser.error('--seed must be between 0 and 2**32 - 1.')
    if args.dpi < 1:
        parser.error('--dpi must be positive.')
    try:
        run(args)
    except (OSError, ValueError, TypeError, KeyError, EOFError, pickle.UnpicklingError) as error:
        parser.exit(1, f'Error: {error}\n')
    except ImportError as error:
        parser.exit(1, f'Missing Python dependency: {error}\nInstall the project requirements.txt.\n')


if __name__ == '__main__':
    main()
