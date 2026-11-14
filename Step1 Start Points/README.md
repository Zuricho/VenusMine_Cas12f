# Step 1: Start points

## Requirements

- A protein structure in PDB format for the Foldseek query.
- A protein FASTA containing reference sequences for embedding analysis.

## Included examples

| File | Use |
|---|---|
| [`8J12.pdb`](8J12.pdb) | Original complex for inspecting the starting structure without downloading it |
| `AsCas12f.pdb` | Prepared protein-chain query for the AsCas12f example |
| `cas12_refseq.fasta` | 17 reference sequences; `8J12-Cas12f` is the target reference used in Step 5 |

The original complex is supplied from [RCSB PDB: 8J12](https://www.rcsb.org/structure/8J12) alongside the prepared protein query. For the Cas12f example, no structure download or preparation is required: use `AsCas12f.pdb` for a new Foldseek search, or use the bundled hit tables in Step 2 to skip the search entirely.

## Prepare your own starting point

1. Select the protein chain of interest from an experimental or predicted structure.
2. Remove unrelated chains, solvent, nucleic acids, and ligands before searching.
3. Check for missing residues and verify that the structure represents the intended protein/domain.
4. Prepare a reference FASTA with unique IDs in the first whitespace-delimited header field. Use the same IDs in the Step 4 embeddings, and select your target in Step 5 with `tree.py --target-reference YOUR_ID`.

Upload the prepared PDB in [Step 2](../Step2%20FoldSeek%20Search/README.md). Use the reference FASTA in [Step 4](../Step4%20ProstT5%20Representation/README.md).

Foldseek query positions refer to the searched sequence, not necessarily native PDB residue numbering. Choose alignment-filter coordinates for your query rather than copying the AsCas12f example blindly.
