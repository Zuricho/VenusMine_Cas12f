# Step 6: Selected candidates

Use sequence-cluster and embedding results to prioritize proteins for further review. This step is **not an automated activity predictor**.

## Use the supplied candidate example

| File | Use |
|---|---|
| [`selected_candidates.fasta`](selected_candidates.fasta) | 17 candidate protein sequences with accession IDs and descriptions |
| [`candidate_ids.txt`](candidate_ids.txt) | The same 17 IDs, one per line, for sequence extraction or accession lookup |

The panel includes `MBP5514557.1`. **Use the FASTA directly without running any preceding step or retrieving sequences from NCBI.** All 17 IDs and sequences match entries in the supplied Step 4 discovery pickle. Step 5's `tree.py` reads `candidate_ids.txt` by default for annotations only; it does not select candidates or read the candidate FASTA. Selection alone does not imply that every listed protein has experimentally confirmed activity.

To generate Step 5 embedding-cluster memberships and the target-cluster FASTA from the supplied embeddings, run from the project root:

```bash
python "Step5 Representation Tree/tree.py"
```

Inspect `Step5 Representation Tree/outputs/sequence_clusters.tsv` (including its `is_candidate` and `hierarchy_cluster` columns) and `target_cluster.fasta` in the same output directory. If [SeqKit](https://bioinf.shenwei.me/seqkit/) is already installed, extract panel IDs that occur in the generated target cluster:

```bash
mkdir -p "Step6 Selected Candidates/outputs"
seqkit grep -f "Step6 Selected Candidates/candidate_ids.txt" \
  "Step5 Representation Tree/outputs/target_cluster.fasta" \
  > "Step6 Selected Candidates/outputs/example_candidates.fasta"
```

The extracted FASTA follows target-cluster order and has accession-only headers. Only IDs present in that cluster are returned: annotating a candidate does not force it into the target cluster. You do not need to run these commands to use the supplied `selected_candidates.fasta`.

For a new analysis, keep the `outputs/hits_s50_cluster.tsv` produced by MMseqs2 in Step 3 if you want to inspect related sequence-cluster members. The embeddings describe representatives only; retrieving nonrepresentative sequences is a separate task.

## Review your own candidates

1. Run Step 5 with your `--references`, `--discovered`, and `--target-reference` options, plus `--no-candidates` to disable the supplied panel.
2. Inspect `sequence_clusters.tsv`, `ordered_labels.csv`, and `target_cluster.fasta` in your Step 5 output directory.
3. Use your Step 3 `outputs/hits_s50_cluster.tsv` mapping, if available, to inspect other members of promising sequence clusters. Sequence-cluster memberships are distinct from the embedding clusters in Step 5.
4. Check full-length sequences, domain architecture, conserved residues, taxonomic assignments, and genomic context. For Cas systems, assess nearby CRISPR-associated features where genome/contig information is available.
5. Select a nonredundant panel for your application and export those sequences to a new FASTA. IDs must be unique first header tokens.

For a new panel, create your own `my_candidate_ids.txt` in the project root with one chosen ID per line. If [SeqKit](https://bioinf.shenwei.me/seqkit/) is already available:

```bash
mkdir -p "Step6 Selected Candidates/outputs"
seqkit grep -f my_candidate_ids.txt \
  "Step5 Representation Tree/outputs/target_cluster.fasta" \
  > "Step6 Selected Candidates/outputs/my_candidates.fasta"
```

Rerun your Step 5 command with `--candidate-ids my_candidate_ids.txt` instead of `--no-candidates`, keeping the same embedding paths and target reference. Use your actual output directory in the extraction command if it differs from the default above. The script matches IDs only, warns about absent IDs, and does not validate an external candidate FASTA against the embeddings; verify sequence consistency yourself. If selecting a nonrepresentative cluster member, update the representative set and compute its embedding first; changing only a plot label is not sufficient.

[Step 7](../Step7%20Structure%20Comparison/README.md) provides structures for all 17 supplied candidates and a pairwise similarity heatmap. Experimental validation is separate from this computational workflow. Genomic-neighborhood screening, accession retrieval, and laboratory assays are not implemented here.
