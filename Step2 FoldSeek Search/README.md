# Step 2: Foldseek search

Find structurally related proteins and export their full sequences for MMseqs2.

## Requirements

- Python 3.10+; the filtering script uses only the standard library.
- For a new search only: a prepared PDB from Step 1 and access to the [Foldseek web server](https://search.foldseek.com/search).

## 1. Use the supplied Cas12f results

The original AsCas12f hit tables are included in `foldseek_results/`. **You can proceed directly to filtering below without running Foldseek or downloading results yourself.**

```text
Step2 FoldSeek Search/foldseek_results/
├── alis_afdb-proteome.m8
├── alis_afdb-swissprot.m8
├── alis_afdb50.m8
├── alis_cath50.m8
├── alis_gmgcl_id.m8
├── alis_mgnify_esm30.m8
└── alis_pdb100.m8
```

These are example results from the supplied AsCas12f query, not searches against today's database releases. `alis_cath50.m8` is intentionally empty: that search returned no hits. The script accepts empty tables and processes all `alis_*.m8` files in the input directory.

**Input format:** headerless Foldseek web-server tables with 19 or 21 tab-separated fields. The script reads target ID (field 2), query start/end (7/8), query length (14), and full target sequence (19), using 1-based field numbers. A generic 12-column alignment table is not sufficient. Check the export schema if the web server changes; do not substitute aligned fragments for full target sequences.

## 2. Filter and export

From the project root:

```bash
python "Step2 FoldSeek Search/filter_foldseek_hits.py" \
  --input-dir "Step2 FoldSeek Search/foldseek_results" \
  --output "Step2 FoldSeek Search/outputs/foldseek_hits.fasta" \
  --protected-start 217 --protected-end 409
```

This keeps hits with `qstart <= 217` and `qend >= 409`. These are example AsCas12f query positions, **not general-purpose thresholds or a claim about historical search settings**. The endpoint test does not require gap-free coverage within that span, and it does not add an E-value or identity cutoff. Review hit quality and adjust the filter for your application.

Identical target sequences are deduplicated across tables. Output IDs are unique numbered IDs, with the database/accession retained in the header description. This avoids ambiguous FASTA IDs when the same accession appears more than once.

With the bundled results and default filter, the script reads **2,046 hits** and exports **1,673 unique sequences**.

**Output:** `outputs/foldseek_hits.fasta`, used as the query set in [Step 3](../Step3%20MMseqs%20Search%20and%20Cluster/README.md). This uniquely named, deduplicated FASTA is generated locally.

A historical export, [`FoldSeek_AsCas12f_foldseek.fasta`](FoldSeek_AsCas12f_foldseek.fasta), is also bundled. It contains **1,682 records but only 1,673 unique sequences**, matching the sequence set produced by the current default filter. Its headers reuse database names as their first tokens, so they are not unique protein IDs. Use the generated `outputs/foldseek_hits.fasta` for downstream commands rather than passing the historical FASTA unchanged.

## Optional: search for your own protein

Upload your prepared PDB to the [Foldseek web server](https://search.foldseek.com/search), select appropriate structure databases, and download the hit tables. Put the extracted `alis_*.m8` files in a separate folder such as `Step2 FoldSeek Search/hits/`, then change `--input-dir` to that folder and choose filter coordinates appropriate to your query. Do not mix new results with the bundled Cas12f example.

Database names and availability can change. All seven example databases are not required for a new search; check that downloaded tables follow the input format above.
