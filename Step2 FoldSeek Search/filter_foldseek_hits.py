#!/usr/bin/env python3
"""Filter Foldseek web-server hit tables and write a FASTA for MMseqs2."""
import argparse
from pathlib import Path


def filter_hits(input_dir, protected_start=217, protected_end=409):
    if not 1 <= protected_start <= protected_end:
        raise ValueError('Require 1 <= protected_start <= protected_end.')
    paths = sorted(Path(input_dir).glob('alis_*.m8'))
    if not paths:
        raise ValueError(f'No alis_*.m8 files in {input_dir}.')
    records, seen, total = [], set(), 0
    for path in paths:
        database = path.stem.removeprefix('alis_')
        with path.open() as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.strip():
                    continue
                fields = line.rstrip('\r\n').split('\t')
                if len(fields) not in (19, 21):
                    raise ValueError(f'{path}:{line_number}: expected 19 or 21 fields, including full target sequence.')
                qstart, qend, qlen = map(int, (fields[6], fields[7], fields[13]))
                if not 1 <= qstart <= qend <= qlen:
                    raise ValueError(f'{path}:{line_number}: invalid query coordinates.')
                total += 1
                if qstart > protected_start or qend < protected_end:
                    continue
                sequence = fields[18].strip().upper()
                if not sequence or not sequence.isascii() or not sequence.isalpha():
                    raise ValueError(f'{path}:{line_number}: invalid target sequence in field 19.')
                if sequence in seen:
                    continue
                seen.add(sequence)
                # A numbered first token stays unique even when accessions repeat.
                accession = fields[1].split()[0]
                records.append((f'hit_{len(records)+1:06d} {database}:{accession}', sequence))
    if not records:
        raise ValueError('No hits passed the endpoint filter; check coordinates and input tables.')
    return records, total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--protected-start', type=int, default=217)
    parser.add_argument('--protected-end', type=int, default=409)
    args = parser.parse_args()
    records, total = filter_hits(args.input_dir, args.protected_start, args.protected_end)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('w') as handle:
        for label, sequence in records:
            handle.write(f'>{label}\n{sequence}\n')
    print(f'Read {total} hits; wrote {len(records)} unique sequences to {args.output}.')


if __name__ == '__main__':
    main()
