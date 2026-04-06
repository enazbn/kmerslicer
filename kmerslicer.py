"""
Protein-focused k-mer slicer with streaming, filters, and multiple outputs.

Usage:
  python kmerslicer3.py input.fasta output.txt -k 9 --index median --format among

Formats:
  - among: kmer;accession; positions (positions are zero-padded, ranges compressed)
  - csv:   kmer,accession,positions (positions are space-separated or ranges)
  - json:  one JSON object per line (NDJSON) with {kmer, accession, positions}
"""
from __future__ import annotations
import argparse, sys, gzip, io, json
from typing import Dict, List, Iterable, Optional
from collections import defaultdict
from Bio import SeqIO

ALLOWED_AA = set("ACDEFGHIKLMNPQRSTVWY") # 20 standard amino acids
def _open_maybe_gz(path: str, mode:str):
    if path == "-":
        return sys.stdin if "r" in mode else sys.stdout
    if path.endswith(".gz"):
        return gzip.open(path, mode+"t")
    return open(path, mode)

def compute_index(start_zero_based: int, kmer_size: int, mode: str) -> int:
    #compute the index of the k-mer
    if mode == "left":
        return start_zero_based + 1
    if mode == "right":
        return start_zero_based + kmer_size
    # median (default)
    return start_zero_based + (kmer_size // 2) + 1


def generate_kmers_with_positions(sequence: str, kmer_size: int, index_mode: str, skip_ambiguous: bool) -> Dict[str, List[int]]:
    #generate the k-mers with their positions
    kmers_to_positions: Dict[str, List[int]] = defaultdict(list)
    seq_len = len(sequence)
    if kmer_size > seq_len:
        return {}
    for i in range(seq_len - kmer_size + 1):
        kmer = sequence[i:i + kmer_size]
        if skip_ambiguous and any((c not in ALLOWED_AA) for c in kmer):
            continue
        kmers_to_positions[kmer].append(compute_index(i, kmer_size, index_mode))
    return kmers_to_positions


def compress_positions(positions: List[int], pad_width: int = 7) -> List[str]:
    #compress the positions
    if not positions:
        return []
    positions_sorted = sorted(positions)
    ranges, start, prev = [], positions_sorted[0], positions_sorted[0]
    for p in positions_sorted[1:]:
        if p == prev + 1:
            prev = p
            continue
        
        ranges.append(f"{start:0{pad_width}d}" if start == prev else f"{start:0{pad_width}d}-{prev:0{pad_width}d}")
        start, prev = p, p
    ranges.append(f"{start:0{pad_width}d}" if start == prev else f"{start:0{pad_width}d}-{prev:0{pad_width}d}")
    return ranges



def write_csv_header(out_handle):
    #write the header of the csv file
    out_handle.write("kmer,accession,positions\n")


def write_csv(out_handle, accession: str, kmers_dict: Dict[str, List[int]], 
min_freq: int, max_freq: Optional[int], pad_width: int, sort_kmers:bool): 
    import csv
    writer = csv.writer(out_handle)
    it=sorted(kmers_dict.keys()) if sort_kmers else kmers_dict.keys()
    for kmer in it:
        positions = kmers_dict[kmer]
        if len(positions) < min_freq or (max_freq is not None and len(positions) > max_freq):
            continue
        writer.writerow([kmer, accession, " ".join(compress_positions(positions, pad_width))])
    #write the k-mers with their positions in csv format


def write_jsonl(out_handle, accession: str, kmers_dict: Dict[str, List[int]], 
min_freq: int, max_freq: Optional[int], pad_width: int, sort_kmers:bool): 
    #write the k-mers with their positions in jsonl format
    it = sorted(kmers_dict.keys()) if sort_kmers else kmers_dict.keys()
    for kmer in it:
        positions = kmers_dict[kmer]
        if len(positions) < min_freq or (max_freq is not None and len(positions) > max_freq):
            continue
        obj = {
            "kmer": kmer,
            "accession": accession,
            "positions": " ".join(compress_positions(positions, pad_width))
        }
        out_handle.write(json.dumps(obj, separators=(",", ":")) + "\n")


def process_stream(input_path: str, output_path: str, kmer_size: int, 
index_mode: str, min_freq: int, max_freq: Optional[int], out_format: str, 
skip_ambiguous: bool, show_progress: bool, pad_width: int, sort_kmers:bool):
    #streams the FASTA file and writes the k-mers with their positions in the specified format
    if kmer_size <= 0:
        raise ValueError("k-mer size must be positive")
    with _open_maybe_gz(output_path, "w") as out_f:
        if out_format == "csv":
            write_csv_header(out_f)
        with _open_maybe_gz(input_path, "r") as handle:
            iterator = SeqIO.parse(handle, "fasta")
            if show_progress:
                try:
                    from tqdm import tqdm  # type: ignore
                    iterator = tqdm(iterator, unit="seq", desc="Processing")
                except Exception:
                    # tqdm not available; continue without progress bar
                    pass
            for record in iterator:
                seq = str(record.seq)
                if len(seq) < kmer_size:
                    continue
                kmers = generate_kmers_with_positions(seq, kmer_size, index_mode, skip_ambiguous)
                if out_format == "csv":
                    write_csv(out_f, record.id, kmers, min_freq, max_freq, pad_width=pad_width, sort_kmers=sort_kmers)
                elif out_format == "json":
                    write_jsonl(out_f, record.id, kmers, min_freq, max_freq, pad_width=pad_width, sort_kmers=sort_kmers)
                else:
                    raise ValueError("unsupported format: " + out_format)


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    #parse the arguments
    p = argparse.ArgumentParser(description="Protein k-mer slicer with streaming and multiple outputs")
    p.add_argument("input", help="Input FASTA file")
    p.add_argument("output", help="Output file path")
    p.add_argument("-k", "--kmer-size", type=int, default=9, help="k-mer size (default: 9)")
    p.add_argument("--index", choices=["left", "median", "right"], default="median", help="index position mode (default: median)")
    p.add_argument("--min-frequency", type=int, default=1, help="minimum frequency to include (default: 1)")
    p.add_argument("--max-frequency", type=int, default=0, help="maximum frequency to include (0=unlimited)")
    p.add_argument("--format", choices=["csv", "json"], default="csv", help="output format (default: csv)")
    p.add_argument("--skip-ambiguous", action="store_true", help="Skip k-mers containing non-standard amino acids (filters B, J, O, U, X, Z, etc.)")
    p.add_argument("--progress", action="store_true", help="Show a progress bar (requires tqdm; otherwise ignored)")
    p.add_argument("--pad-width", type=int, default=7, help="padding width for positions (default: 7)")
    p.add_argument("--sort", action="store_true", help="sort k-mers alphabetically (default: False)")
    return p.parse_args(list(argv))


def main(argv: Iterable[str]) -> int:
    #main function
    args = parse_args(argv)
    max_freq = None if args.max_frequency == 0 else args.max_frequency
    try:
        process_stream(
            input_path=args.input,
            output_path=args.output,
            kmer_size=args.kmer_size,
            index_mode=args.index,
            min_freq=args.min_frequency,
            max_freq=max_freq,
            out_format=args.format,
            skip_ambiguous=args.skip_ambiguous,
            show_progress=args.progress,
            pad_width=args.pad_width,
            sort_kmers=args.sort,
        )
    except Exception as e:
        print(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))


