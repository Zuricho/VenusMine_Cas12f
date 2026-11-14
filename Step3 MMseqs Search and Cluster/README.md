# Step 3: MMseqs2 search and clustering

Expand the Foldseek hits against a sequence database, then reduce redundancy with sequence clustering.

## Use the supplied representatives

For the supplied Cas12f analysis, skip MMseqs2 and model inference: [Step 4](../Step4%20ProstT5%20Representation/README.md) includes `AsCas12f_NR_s50.pkl` with 17,202 representative IDs, sequences, and embeddings. Use that file directly in Step 5, or export its sequences as shown below. A separate copy of the representative FASTA is unnecessary.

This folder provides `tsv2fasta.py` and instructions for new searches. Historical member-to-representative mappings are intentionally outside the package's scope and cannot be inferred from embeddings alone. For new analyses, generate and retain `hits_s50_cluster.tsv` using the commands below; a new MMseqs2 run may choose different representatives.

Use the supplied pickle's labels and sequences as-is, including `MBP5514557.1`. For a new analysis, verify membership in your own `hits_s50_cluster.tsv` before choosing a different representative; update both the mapping and the representative FASTA, then compute the new protein's embedding. Do not simply rename an old sequence or embedding.

## Requirements for a new search

- [MMseqs2](https://github.com/soedinglab/MMseqs2) available as `mmseqs`.
- A local MMseqs2 protein database, such as NCBI NR. Follow the MMseqs2 documentation to download/build it; databases are not bundled here.
- `Step2 FoldSeek Search/outputs/foldseek_hits.fasta`.
- Python 3.10+ for `tsv2fasta.py` (standard library only).

Large searches may require a compute server, substantial RAM, and temporary disk space. Adjust thread counts for your machine.

## 1. Search

From the project root, replace the database path below:

```bash
DB="/path/to/mmseqs/protein_database"
OUT="Step3 MMseqs Search and Cluster/outputs"
mkdir -p "$OUT"

mmseqs easy-search \
  "Step2 FoldSeek Search/outputs/foldseek_hits.fasta" "$DB" \
  "$OUT/search.tsv" "$OUT/tmp_search" \
  --format-output "query,target,tseq" \
  --threads 8 --num-iterations 2 --max-seqs 2500 --search-type 1
```

These are starting parameters for a new analysis, not a record of historical jobs. Review MMseqs2 sensitivity, E-value, and coverage settings for your target.

## 2. Export full target sequences

```bash
python "Step3 MMseqs Search and Cluster/tsv2fasta.py" \
  "$OUT/search.tsv" "$OUT/search_hits.fasta"
```

The converter expects exactly three headerless columns: `query`, `target`, `tseq`. It keeps one full sequence per target ID, rejects conflicting sequences for the same ID, and does not silently discard different IDs with identical sequences.

## 3. Cluster at 50% identity

```bash
mmseqs easy-cluster "$OUT/search_hits.fasta" "$OUT/hits_s50" \
  "$OUT/tmp_cluster" --min-seq-id 0.5 --threads 8
```

Important outputs:

- `hits_s50_rep_seq.fasta`: cluster representatives for Step 4.
- `hits_s50_cluster.tsv`: representative/member mapping for tracing candidates back to search hits.
- `hits_s50_all_seqs.fasta`: clustered sequences.

Identity alone does not define full-length similarity; select coverage settings appropriate to your data if needed. Keep the representative/member mapping in your own analysis, especially if you choose a different representative from a cluster.

Proceed to [Step 4](../Step4%20ProstT5%20Representation/README.md) with your new `outputs/hits_s50_rep_seq.fasta`, or with `outputs/representatives.fasta` after the optional export above. Recomputing embeddings for the supplied sequences is not required to run Step 5. Preserve `outputs/hits_s50_cluster.tsv` from any new clustering run so you can trace representatives back to their member sequences.
