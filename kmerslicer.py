from __future__ import annotations
import argparse, sys, gzip, json, csv
from typing import Iterable, Optional
from Bio import SeqIO

ALLOWED_AA = set("ACDEFGHIKLMNPQRSTVWY")  # 20 standard amino acids


def _open_maybe_gz(path: str, mode: str):
    if path == "-":
        return sys.stdin if "r" in mode else sys.stdout
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode)


def compute_index(start_zero_based: int, kmer_size: int, mode: str) -> int:
    if mode == "left":
        return start_zero_based + 1
    if mode == "right":
        return start_zero_based + kmer_size
    return start_zero_based + (kmer_size // 2) + 1  # median


def write_csv_header(out_handle):
    writer = csv.writer(out_handle)
    writer.writerow(["kmer", "accession", "position"])


def process_stream(
    input_path: str,
    output_path: str,
    kmer_size: int,
    index_mode: str,
    out_format: str,
    skip_ambiguous: bool,
    show_progress: bool,
    sort_kmers: bool,
):
    if kmer_size <= 0:
        raise ValueError("k-mer size must be positive")

    with _open_maybe_gz(output_path, "w") as out_f:
        csv_writer = None
        if out_format == "csv":
            csv_writer = csv.writer(out_f)
            csv_writer.writerow(["kmer", "accession", "position"])

        with _open_maybe_gz(input_path, "r") as handle:
            iterator = SeqIO.parse(handle, "fasta")

            if show_progress:
                try:
                    from tqdm import tqdm  # type: ignore
                    iterator = tqdm(iterator, unit="seq", desc="Processing")
                except Exception:
                    pass

            for record in iterator:
                accession = record.id
                seq = str(record.seq)
                seq_len = len(seq)

                if seq_len < kmer_size:
                    continue

                rows = []
                for i in range(seq_len - kmer_size + 1):
                    kmer = seq[i:i + kmer_size]

                    if skip_ambiguous and any(c not in ALLOWED_AA for c in kmer):
                        continue

                    pos = compute_index(i, kmer_size, index_mode)
                    rows.append((kmer, accession, pos))

                if sort_kmers:
                    rows.sort(key=lambda x: (x[0], x[2]))

                if out_format == "csv":
                    for kmer, accession, pos in rows:
                        csv_writer.writerow([kmer, accession, pos])

                elif out_format == "json":
                    for kmer, accession, pos in rows:
                        obj = {
                            "kmer": kmer,
                            "accession": accession,
                            "position": pos,
                        }
                        out_f.write(json.dumps(obj, separators=(",", ":")) + "\n")

                else:
                    raise ValueError("unsupported format: " + out_format)


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Protein k-mer slicer with expanded one-row-per-position output"
    )
    p.add_argument("input", help="Input FASTA file")
    p.add_argument("output", help="Output file path")
    p.add_argument("-k", "--kmer-size", type=int, default=9, help="k-mer size (default: 9)")
    p.add_argument(
        "--index",
        choices=["left", "median", "right"],
        default="median",
        help="index position mode (default: median)",
    )
    p.add_argument(
        "--format",
        choices=["csv", "json"],
        default="csv",
        help="output format (default: csv)",
    )
    p.add_argument(
        "--skip-ambiguous",
        action="store_true",
        help="Skip k-mers containing non-standard amino acids",
    )
    p.add_argument(
        "--progress",
        action="store_true",
        help="Show a progress bar (requires tqdm; otherwise ignored)",
    )
    p.add_argument(
        "--sort",
        action="store_true",
        help="Sort rows by k-mer then position within each accession",
    )
    return p.parse_args(list(argv))


def main(argv: Iterable[str]) -> int:
    args = parse_args(argv)
    try:
        process_stream(
            input_path=args.input,
            output_path=args.output,
            kmer_size=args.kmer_size,
            index_mode=args.index,
            out_format=args.format,
            skip_ambiguous=args.skip_ambiguous,
            show_progress=args.progress,
            sort_kmers=args.sort,
        )
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
