#!/usr/bin/env python3
"""Convert MMseqs2 query,target,tseq output to FASTA, keeping one row per target."""
import argparse
from pathlib import Path


def convert(source, destination):
    sequences = {}
    with Path(source).open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) != 3:
                raise ValueError(f'Line {line_number}: expected query,target,tseq (three tab-separated fields).')
            _, accession, sequence = fields
            if not accession or any(c.isspace() for c in accession):
                raise ValueError(f'Line {line_number}: target must be a nonempty ID without spaces.')
            sequence = sequence.strip().upper()
            if not sequence or not sequence.isascii() or not sequence.isalpha():
                raise ValueError(f'Line {line_number}: invalid full target sequence.')
            if accession in sequences and sequences[accession] != sequence:
                raise ValueError(f'Conflicting sequences for {accession}.')
            sequences[accession] = sequence
    if not sequences:
        raise ValueError('No sequences found.')
    destination = Path(destination)
    if Path(source).resolve() == destination.resolve():
        raise ValueError('Input and output must be different files.')
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open('w') as handle:
        for accession, sequence in sequences.items():
            handle.write(f'>{accession}\n{sequence}\n')
    return len(sequences)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    print(f'Wrote {convert(args.input, args.output)} targets to {args.output}.')


if __name__ == '__main__':
    main()
