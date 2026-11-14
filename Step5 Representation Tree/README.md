# Step 5: Representation tree

[`tree.py`](tree.py) is a standalone version of [`tree.ipynb`](tree.ipynb). It combines reference and discovery ProstT5 embeddings, performs hierarchical clustering, and saves a radial dendrogram and overview/detail UMAP plots. The notebook is left unchanged; no configuration file or Jupyter session is required.

## Requirements and inputs

Use Python 3.10+ with the packages in the project-root `requirements.txt`:

```bash
# If these packages are not already installed in your chosen environment:
python -m pip install -r requirements.txt
```

No GPU, model download, or sequence search is needed at this step. Only load **trusted pickle files**: unpickling can execute code.

Both input pickles must use the Step 4 format:

```text
labels           N unique protein ID strings
sequences        N amino-acid sequence strings, in the same order
representations  finite numeric NumPy array with shape (N, D)
```

Use the same embedding model and preprocessing for both inputs. Their feature dimensions must match; IDs must be unique across both files. Surrounding whitespace in labels is stripped, including the trailing spaces used in the notebook. At least three proteins in total are required.

Default inputs, resolved relative to the script rather than the working directory:

- References: `../Step4 ProstT5 Representation/cas12_refseq.pkl`.
- Discovery proteins: `../Step4 ProstT5 Representation/AsCas12f_NR_s50.pkl`, as in the notebook.
- Optional candidate annotations: `../Step6 Selected Candidates/candidate_ids.txt`, if present.

**Both default embedding files are included:** 17,202 discovery proteins and 17 references, for **17,219 proteins in total**. No embedding generation is needed to run this example. This script does not download missing data or substitute a different representative set.

## Run

From the project root, using the supplied embeddings:

```bash
python "Step5 Representation Tree/tree.py"
```

To use discovery embeddings generated in Step 4 instead:

```bash
python "Step5 Representation Tree/tree.py" \
  --discovered "Step4 ProstT5 Representation/outputs/representatives.pkl"
```

For a new dataset, supply your references and the reference ID defining the cluster of interest:

```bash
python "Step5 Representation Tree/tree.py" \
  --references /path/to/references.pkl \
  --discovered /path/to/discovered.pkl \
  --target-reference my_reference_id \
  --no-candidates \
  --output-dir "Step5 Representation Tree/outputs/my_analysis"
```

Explicit relative paths are resolved from your current working directory. Figures are saved without opening a window. Existing result files with the same names in the output directory are overwritten; use a separate directory for each analysis.

Candidate IDs are annotations only, not automatic candidate selection. Use `--candidate-ids /path/to/ids.txt` to supply a different panel (one ID per line; blank lines and `#` comments are ignored), or `--no-candidates` to disable it. IDs absent from the discovery embeddings produce a warning and are skipped. No accession renaming or sequence replacement is performed.

```bash
python "Step5 Representation Tree/tree.py" --help
```

### Analysis settings

| Option | Default | Purpose |
|---|---|---|
| `--threshold` | `10` | Ward distance cut used for colors and target-cluster selection |
| `--coarse-threshold` | `20` | Second distance cut recorded in the membership table |
| `--target-reference` | `8J12-Cas12f` | Reference whose cluster is highlighted and exported |
| `--n-neighbors` | `100` | UMAP neighbors; capped at the number of proteins minus one |
| `--min-dist` | `0.1` | UMAP minimum distance |
| `--spread` | `3` | UMAP spread |
| `--seed` | `42` | Random seed for UMAP |
| `--dpi` | `300` | PNG resolution |
| `--output-dir` | `Step5 Representation Tree/outputs/` | Result directory, script-relative by default |

Clustering uses Euclidean pairwise distances and Ward linkage on the **original embeddings**, with no normalization or dimensionality reduction first. UMAP is a separate, two-dimensional visualization using the notebook's Euclidean metric and neighborhood settings.

Pairwise clustering scales quadratically: for the notebook's 17,219 proteins, the condensed distance array alone occupies about **1.10 GiB**. Allow several GB of available RAM for embeddings, linkage, UMAP, and plotting. The script prints its distance-memory estimate before clustering. Large analyses can take substantial time.

## Outputs

Run the script to generate the following files in `outputs/` by default, or in your chosen `--output-dir`. The notebook retains its original embedded figures, which are separate from these CLI outputs. The inputs and notebook are not modified.

| File | Contents |
|---|---|
| `dendrogram.png`, `dendrogram.pdf` | Radial Ward tree with cluster-color ring and reference/candidate markers |
| `umap.png` | UMAP overview with a rectangle marking the detail view |
| `umap_detail.png` | UMAP zoom fitted to the target cluster |
| `ordered_labels.csv` | Labels and `hierarchy_cluster` IDs in dendrogram leaf order |
| `sequence_clusters.tsv` | All proteins in input order: label, reference/candidate flags, sequence length, fine/coarse cluster IDs, and UMAP coordinates |
| `target_cluster.fasta` | Discovery sequences in the target reference's fine cluster; references are excluded |

The membership table lists references first, followed by discovery proteins. The FASTA preserves discovery-input order and is empty if no discovery proteins share the target cluster. Cluster IDs are dataset-dependent labels, not biological classifications.

Colors follow the notebook: yellow for the target cluster, green for the target reference, coral for other references, and blue for candidates. Other cluster colors cycle through seven pastel colors, so color alone does not uniquely identify a cluster; use the tables for exact membership.

## Differences from the notebook

- The target cluster is looked up using `8J12-Cas12f` (or `--target-reference`), rather than assuming cluster `40` or reference row `14`.
- All references are annotated by their actual input membership; no hard-coded row exclusions or special border for row `13` are used.
- The supplied candidate-ID file replaces repeated hard-coded lists. Its full panel, including `MBP5514557.1` when present in the embeddings, is used consistently in all plots.
- UMAP uses seed `42`, one worker, and random initialization. This supports small datasets and repeatable runs in the same software environment, but does not reproduce the unseeded notebook's coordinates. Results can differ across package versions.
- Overview and detail limits adapt to the current data instead of using fixed notebook coordinates. DejaVu Sans avoids requiring Arial.
- The membership TSV and target-cluster FASTA are additional exports for downstream review.

Embedding proximity is a prioritization aid, not evidence of nuclease activity. Review genomic context and validate candidates experimentally before drawing functional conclusions.
