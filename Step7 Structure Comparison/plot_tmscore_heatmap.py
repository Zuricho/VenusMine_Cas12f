"""Plot precomputed TM-align scores for PDB files in the structures directory."""

import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

# Resolve paths relative to this script, not the current working directory.
BASE_DIR = Path(__file__).resolve().parent
STRUCTURES_DIR = BASE_DIR / 'structures'
INPUT_TSV = BASE_DIR / 'tmalign_pairwise.tsv'
OUTPUT_PNG = BASE_DIR / 'outputs' / 'tmscore_heatmap_clustered.png'

# Plot config
FIGSIZE = (8, 8)
FONT_FAMILY = 'Arial'
FONT_SIZE = 12
V_MIN = 0
V_MAX = 1
COLORBAR_LABEL = 'TM-score'
COLORBAR_WIDTH = '2.5%'
COLORBAR_PAD = 0.1

IMPORTANT = {
    'DAG32506.1',
    'MBQ9950927.1',
    'MBS6459277.1',
    'MCI7585022.1',
    'MBE5919695.1',
}

# Plot style
plt.rcParams['font.family'] = FONT_FAMILY
plt.rcParams['font.size'] = FONT_SIZE

# Load PDB names from structures/, preserving basenames for matching TSV rows.
if not STRUCTURES_DIR.is_dir():
    raise FileNotFoundError(f'Structure directory not found: {STRUCTURES_DIR}')
names = sorted(path.name for path in STRUCTURES_DIR.glob('*.pdb') if path.is_file())
if len(names) < 2:
    raise ValueError(f'At least two PDB files are required in {STRUCTURES_DIR}.')
labels = [Path(name).stem for name in names]
idx = {n: i for i, n in enumerate(names)}

n = len(names)
mat = np.full((n, n), np.nan, dtype=float)

# Diagonal
for i in range(n):
    mat[i, i] = 1.0

# Fill matrix from TSV
with open(INPUT_TSV, newline='') as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        a = row['A']
        b = row['B']
        if a not in idx or b not in idx:
            continue
        try:
            tma = float(row['TMscore_norm_A']) if row['TMscore_norm_A'] else np.nan
            tmb = float(row['TMscore_norm_B']) if row['TMscore_norm_B'] else np.nan
        except ValueError:
            tma = tmb = np.nan
        vals = [v for v in (tma, tmb) if not np.isnan(v)]
        if vals:
            avg = sum(vals) / len(vals)
            i, j = idx[a], idx[b]
            mat[i, j] = avg
            mat[j, i] = avg

# Clustering
# distance = 1 - TMscore (avg). Unknown -> 1.0

dist = 1.0 - mat
np.fill_diagonal(dist, 0.0)
dist = np.where(np.isnan(dist), 1.0, dist)
condensed = squareform(dist, checks=False)
Z = linkage(condensed, method='average')
order = leaves_list(Z)

mat_ord = mat[order][:, order]
labels_ord = [labels[i] for i in order]

# Custom colormap: white until 0.3, gradient to blue at 0.7, blue beyond
cmap = LinearSegmentedColormap.from_list(
    'white_to_blue_thresh',
    [(0.0, '#ffffff'), (0.3, '#ffffff'), (0.7, '#1f4fbf'), (1.0, '#1f4fbf')]
)

fig, ax = plt.subplots(figsize=FIGSIZE)
img = ax.imshow(mat_ord, cmap=cmap, vmin=V_MIN, vmax=V_MAX)

ax.set_xticks(range(n))
ax.set_yticks(range(n))
ax.set_xticklabels(labels_ord, rotation=90, fontsize=FONT_SIZE)
ax.set_yticklabels(labels_ord, fontsize=FONT_SIZE)

# Label colors
for tick in ax.get_xticklabels():
    txt = tick.get_text()
    if txt in IMPORTANT:
        tick.set_color('blue')
    elif txt.endswith('.1'):
        tick.set_color('red')

for tick in ax.get_yticklabels():
    txt = tick.get_text()
    if txt in IMPORTANT:
        tick.set_color('blue')
    elif txt.endswith('.1'):
        tick.set_color('red')

# Colorbar aligned to plot height

divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size=COLORBAR_WIDTH, pad=COLORBAR_PAD)
cb = fig.colorbar(img, cax=cax)
cb.set_label(COLORBAR_LABEL)

plt.tight_layout()
OUTPUT_PNG.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(OUTPUT_PNG, dpi=300)
plt.close(fig)
print(f"Saved {OUTPUT_PNG}")
