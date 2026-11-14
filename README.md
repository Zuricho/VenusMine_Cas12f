# Cas12f Discovery

A structure-search and protein-language-model workflow for exploring Cas12f candidates, organized in the step-by-step style of [VenusMine](https://github.com/ai4protein/VenusMine).

**Structure → Foldseek hits → MMseqs2 search and clustering → ProstT5 embeddings → representation tree → candidate review → structure comparison**

Please refer to the README in each step:

| Step | Purpose |
|---|---|
| [1. Start points](Step1%20Start%20Points/README.md) | Prepare a query structure and reference sequences |
| [2. Foldseek search](Step2%20FoldSeek%20Search/README.md) | Find structural homologs and export filtered sequences |
| [3. MMseqs2 search and clustering](Step3%20MMseqs%20Search%20and%20Cluster/README.md) | Expand the sequence set and select representatives |
| [4. ProstT5 representations](Step4%20ProstT5%20Representation/README.md) | Generate protein embeddings |
| [5. Representation tree](Step5%20Representation%20Tree/README.md) | Cluster embeddings, plot UMAP, and export sequences |
| [6. Selected candidates](Step6%20Selected%20Candidates/README.md) | Review candidates and use the supplied sequence panel |
| [7. Structure comparison](Step7%20Structure%20Comparison/README.md) | Compare candidate and reference structures using TM-align scores |
