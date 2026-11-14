# Step 7: Structure comparison

Visualize pairwise TM-align similarity between candidate and reference structures.

## Files

- `structures/`: 33 PDB files (17 candidates and 16 references).
- `tmalign_pairwise.tsv`: scores for all 528 unique pairs.
- `tmscore_heatmap_clustered.png`: supplied example heatmap.
- `plot_tmscore_heatmap.py`: generates a clustered heatmap from the existing scores.

## Run

Requires Python 3.10+, NumPy, SciPy, and Matplotlib. From the project root:

```bash
# Install dependencies if needed:
python -m pip install numpy scipy matplotlib
python "Step7 Structure Comparison/plot_tmscore_heatmap.py"
```

The script reads `structures/` and `tmalign_pairwise.tsv` relative to its own location. The figure is saved to `Step7 Structure Comparison/outputs/tmscore_heatmap_clustered.png`; rerunning overwrites this output, not the supplied example.

## Notes

- Cells show the mean of the two length-normalized TM-scores. Ordering uses average-linkage clustering on `1 − mean TM-score`.
- The script **does not run TM-align**. For new structures, select the intended protein chains, calculate pairwise scores separately, and update the TSV. Its `A` and `B` columns must match PDB filenames, including `.pdb`.
- Edit the script's path constants, font settings, or `IMPORTANT` IDs to customize the plot. Structural similarity alone does not establish biological activity.
