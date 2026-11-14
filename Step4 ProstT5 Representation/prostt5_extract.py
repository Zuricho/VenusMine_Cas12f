#!/usr/bin/env python3
"""Compute mean-pooled amino-acid ProstT5 embeddings from a protein FASTA."""
import argparse
from pathlib import Path
import pickle
import re


def read_fasta(path):
    records, label, parts = {}, None, []

    def store():
        if label is None:
            return
        sequence = ''.join(parts).upper()
        if label in records:
            raise ValueError(f'Duplicate FASTA ID: {label}')
        if not sequence or any(aa not in 'ACDEFGHIKLMNPQRSTVWYXBZOU' for aa in sequence):
            raise ValueError(f'{label}: expected a nonempty amino-acid sequence (no gaps or stop symbols).')
        records[label] = sequence

    with Path(path).open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                store()
                if not line[1:].strip():
                    raise ValueError('Empty FASTA header.')
                label, parts = line[1:].split()[0], []
            else:
                if label is None:
                    raise ValueError('Sequence before first FASTA header.')
                parts.append(line)
    store()
    if not records:
        raise ValueError('Empty FASTA.')
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--device', default='auto', help='auto, cpu, or a PyTorch device such as cuda:0')
    parser.add_argument('--batch-size', type=int, default=1)
    parser.add_argument('--model', default='Rostlab/ProstT5', help='Hugging Face model ID or local model directory')
    args = parser.parse_args()
    if args.batch_size < 1:
        parser.error('--batch-size must be positive.')
    if args.input.resolve() == args.output.resolve():
        parser.error('Input and output must be different files.')
    records = read_fasta(args.input)

    # Keep --help and FASTA validation usable without downloading model packages.
    import numpy as np
    import torch
    from transformers import T5EncoderModel, T5Tokenizer

    device = torch.device(('cuda:0' if torch.cuda.is_available() else 'cpu') if args.device == 'auto' else args.device)
    tokenizer = T5Tokenizer.from_pretrained(args.model, do_lower_case=False)
    model = T5EncoderModel.from_pretrained(args.model).to(device)
    model = model.half() if device.type == 'cuda' else model.float()
    model.eval()
    labels, sequences = list(records), list(records.values())
    vectors = []
    for start in range(0, len(sequences), args.batch_size):
        batch = sequences[start:start+args.batch_size]
        # Keep original sequences in the output; normalize rare residues only for the model.
        texts = ['<AA2fold> ' + ' '.join(re.sub('[UZOB]', 'X', sequence)) for sequence in batch]
        tokens = tokenizer(texts, add_special_tokens=True, padding=True, return_tensors='pt').to(device)
        expected_lengths = [len(sequence)+2 for sequence in batch]  # prefix + residues + EOS
        if tokens.attention_mask.sum(dim=1).tolist() != expected_lengths:
            raise ValueError('Unexpected tokenization; check the model/tokenizer and protein alphabet.')
        with torch.inference_mode():
            hidden = model(input_ids=tokens.input_ids, attention_mask=tokens.attention_mask).last_hidden_state
            for i, sequence in enumerate(batch):
                vectors.append(hidden[i, 1:len(sequence)+1].float().mean(dim=0).cpu().numpy())
        print(f'Embedded {min(start+args.batch_size, len(sequences))}/{len(sequences)}', flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('wb') as handle:
        pickle.dump({'labels': labels, 'sequences': sequences, 'representations': np.stack(vectors)}, handle)
    print(f'Saved {len(labels)} embeddings to {args.output}.')


if __name__ == '__main__':
    main()
