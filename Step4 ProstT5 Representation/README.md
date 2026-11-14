# Step 4: ProstT5 representations

Convert protein sequences into fixed-length vectors using [Rostlab/ProstT5](https://huggingface.co/Rostlab/ProstT5).

## Use the supplied embeddings

The two current embedding files are included as ready-to-use Cas12f examples. **You can skip Steps 1–3, model installation, model downloads, and embedding computation**, and go directly to Step 5 with:

| File | Contents |
|---|---|
| [`AsCas12f_NR_s50.pkl`](AsCas12f_NR_s50.pkl) | 17,202 search proteins: IDs, sequences, and 1,024-dimensional embeddings |
| `cas12_refseq.pkl` | 17 reference sequences and their embeddings |

The supplied discovery pickle is ready to use as-is, including `MBP5514557.1`. Its stored labels and sequences define the authoritative row mapping; no additional representative mapping or override is needed for Step 5. Do not assume a new MMseqs2 run selects identical representatives.

Both supplied embedding arrays use **float16**. The discovery array has shape **(17,202, 1,024)**; sequence lengths range from **26 to 1,000 amino acids**, with median **202**. All discovery labels have trailing whitespace, which `tree.py` strips when matching IDs. The pickle contains no taxonomy, cluster-member mapping, functional validation, or structure coordinates.

## Example sequences for embedding generation

The supplied [`../Step1 Start Points/cas12_refseq.fasta`](../Step1%20Start%20Points/cas12_refseq.fasta) contains the same 17 reference sequences as the reference pickle. After installing the model dependencies below, you can try embedding generation on this small input, from the project root:

```bash
python "Step4 ProstT5 Representation/prostt5_extract.py" \
  "Step1 Start Points/cas12_refseq.fasta" \
  "Step4 ProstT5 Representation/outputs/example_references.pkl"
```

The discovery sequences are already stored in the pickle. To use them as FASTA, follow the optional export in [Step 3](../Step3%20MMseqs%20Search%20and%20Cluster/README.md) to create `Step3 MMseqs Search and Cluster/outputs/representatives.fasta`.

Recomputing embeddings is optional and still requires downloading/loading ProstT5. Use the existing pickles if you only want to continue the Cas12f analysis. New outputs go under `outputs/` so the supplied inputs remain untouched.

## Generate embeddings for new sequences

Use Python 3.10+. In an environment intended for embedding generation, install the additional requirements:

```bash
python -m pip install -r "Step4 ProstT5 Representation/requirements.txt"
```

For CUDA, install a PyTorch build appropriate for your system using the [PyTorch instructions](https://pytorch.org/get-started/locally/). A GPU is recommended; CPU inference uses full precision and can be slow and memory-intensive. The first run downloads the multi-GB model to the Hugging Face cache. Set `HF_HOME` to a directory with sufficient space if necessary.

From the project root:

```bash
python "Step4 ProstT5 Representation/prostt5_extract.py" \
  "Step3 MMseqs Search and Cluster/outputs/hits_s50_rep_seq.fasta" \
  "Step4 ProstT5 Representation/outputs/representatives.pkl"

python "Step4 ProstT5 Representation/prostt5_extract.py" \
  "Step1 Start Points/cas12_refseq.fasta" \
  "Step4 ProstT5 Representation/outputs/references.pkl"
```

Options: `--device cpu`, `--device cuda:0`, `--batch-size 1`, or `--model /path/to/local/ProstT5`. Smaller batches reduce memory use. No sequences are silently truncated or excluded by length; exceptionally long proteins may exhaust memory.

Input FASTA IDs must be unique first header tokens. Gaps, stop symbols, and empty sequences are rejected. Rare residues `U/Z/O/B` are mapped to `X` for model input only; the original uppercase sequences are retained in the output.

The encoder receives the `<AA2fold>` prefix and spaced amino acids. Residue embeddings are mean-pooled, excluding the prefix, end token, and padding. New outputs use float32 vectors; CUDA model inference uses float16. New computations are not expected to match the supplied embeddings bit-for-bit.

## Output format

Each pickle contains a dictionary:

```text
labels           list of N unique protein IDs
sequences        list of N amino-acid sequences in the same order
representations  NumPy array of shape (N, 1024)
```

For each index `i`, `labels[i]`, `sequences[i]`, and `representations[i]` describe the same protein. Use the same model and preprocessing for references and search representatives. Pass new paths to [Step 5](../Step5%20Representation%20Tree/README.md) using `tree.py --references PATH --discovered PATH`; no JSON configuration is needed. **Only load trusted pickle files**, because unpickling can execute code.
