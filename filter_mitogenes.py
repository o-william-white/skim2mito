#!/usr/bin/env python3
"""
filter_mitogenes.py

Filters sequences from mitochondrial gene FASTA files based on four criteria
(applied in order; a sequence must pass all four to be retained):

  1. Minimum length: sequence must be >= half the average insect CDS/gene length.
  2. Non-fragmented: sequences with a trailing _N suffix on the gene name in the
     header (e.g. rrnS_2, rrnS_1) are considered fragmented and excluded.
  3. Multi-contig duplicate: if the same gene appears on more than one contig
     for the same sample after the length filter, all copies are excluded.
     Running this after the length filter means a short copy on one contig is
     already gone before duplicates are assessed, correctly resolving the common
     case where one copy is too short to be useful.
  4. Taxonomy: the contig from which the sequence was annotated must have
     Class == 'Insecta' in the summary TSV. Sequences from contigs with any
     other Class value (e.g. 'no-hit', 'Malacostraca') are excluded.
     Only the Class column is checked; higher-rank no-hits are ignored.

Expected FASTA header format:
    >{sample}_{contig};{start}-{end};{strand};{gene}[_{fragment_index}]
  e.g.
    >UK001-G01_contig0;545-592;+;rrnS_2   <- fragmented, excluded
    >UK001-G01_contig0;714-741;+;rrnS     <- kept if passes other filters

Summary TSV format (tab-separated, produced by assess_mitogenomes.py):
    ID  Index  Contig  GC  Length  Coverage  ...  Class  ...  Genes list
  'Genes list' holds comma-separated gene names for each contig.
  Rows where Index is 'NA' (no assembly) are skipped.

Usage:
    python filter_mitogenes.py \\
        --input   /path/to/fasta/dir \\
        --genes   cox1,cox2,cob,rrnL,rrnS \\
        --summary /path/to/summary.tsv \\
        [--output /path/to/output/dir]
"""

import argparse
import os
import re
import sys
from collections import defaultdict


# ---------------------------------------------------------------------------
# Reference data
# ---------------------------------------------------------------------------

INSECT_GENE_LENGTHS = {
    "atp6":  681,
    "atp8":  156,
    "cob":   1141,
    "cox1":  1539,
    "cox2":  690,
    "cox3":  789,
    "nad1":  939,
    "nad2":  1041,
    "nad3":  357,
    "nad4":  1341,
    "nad4l": 297,
    "nad5":  1740,
    "nad6":  522,
    "rrnl":  1300,
    "rrns":  790,
}

CANONICAL_NAMES = {
    "atp6":  "atp6",
    "atp8":  "atp8",
    "cob":   "cob",
    "cox1":  "cox1",
    "cox2":  "cox2",
    "cox3":  "cox3",
    "nad1":  "nad1",
    "nad2":  "nad2",
    "nad3":  "nad3",
    "nad4":  "nad4",
    "nad4l": "nad4l",
    "nad5":  "nad5",
    "nad6":  "nad6",
    "rrnl":  "rrnL",
    "rrns":  "rrnS",
}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Filter mitochondrial gene FASTAs by minimum length, fragmentation, "
            "multi-contig duplicate detection, and taxonomy (Class == Insecta)."
        )
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Directory containing per-gene multi-FASTA files (e.g. cox1.fasta, cob.fasta)."
    )
    parser.add_argument(
        "--genes", "-g",
        required=True,
        help="Comma-separated list of gene names to process (e.g. cox1,cox2,cob,rrnL,rrnS)."
    )
    parser.add_argument(
        "--summary", "-s",
        required=True,
        help=(
            "Path to the assembly summary TSV file (produced by assess_mitogenomes.py). "
            "Used for taxonomy filtering (Class == Insecta)."
        )
    )
    parser.add_argument(
        "--output", "-o",
        default=None,
        help="Output directory for filtered FASTAs. Defaults to <input_dir>/filtered."
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# FASTA I/O
# ---------------------------------------------------------------------------

def read_fasta(filepath):
    """Parse a FASTA file; yield (header, sequence) tuples."""
    header    = None
    seq_parts = []
    with open(filepath, "r") as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header    = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        yield header, "".join(seq_parts)


def reformat_header(header):
    """Replace the underscore between sample ID and contig ID with a pipe.

    e.g. UK001-G01_contig0;545-592;+;rrnS  ->  UK001-G01|contig0;545-592;+;rrnS
    """
    parsed = parse_header(header)
    if parsed is None:
        return header
    old_prefix = f"{parsed['sample']}_{parsed['contig']}"
    new_prefix  = f"{parsed['sample']}|{parsed['contig']}"
    return header.replace(old_prefix, new_prefix, 1)


def write_fasta(sequences, filepath):
    """Write (header, sequence) tuples to a FASTA file with reformatted headers."""
    with open(filepath, "w") as fh:
        for header, seq in sequences:
            fh.write(f">{reformat_header(header)}\n{seq}\n")


# ---------------------------------------------------------------------------
# Header parsing
# ---------------------------------------------------------------------------

def parse_header(header):
    """Parse a FASTA header of the form:
        {sample}_{contig};{start}-{end};{strand};{gene_field}

    Returns a dict with keys:
        sample      - e.g. "UK001-G01"
        contig      - e.g. "contig0"
        gene_field  - raw gene field, e.g. "rrnS_2" or "rrnS"
        gene_base   - gene name with trailing _<int> stripped, lowercased
        is_fragment - True if gene_field ends with _<integer>
        sample_id   - full sample+contig string, e.g. "UK001-G01_contig0"

    Returns None if the header does not match the expected format.
    """
    parts = header.split(";")
    if len(parts) < 4:
        return None

    sample_contig = parts[0]
    gene_field    = parts[3]

    match = re.match(r'^(.+)_(contig\d+|circular)$', sample_contig, re.IGNORECASE)
    if not match:
        return None

    sample = match.group(1)
    contig = match.group(2)

    is_fragment = bool(re.search(r'_\d+$', gene_field))
    gene_base   = re.sub(r'_\d+$', '', gene_field).lower()

    return {
        "sample":      sample,
        "contig":      contig,
        "sample_id":   sample_contig,
        "gene_field":  gene_field,
        "gene_base":   gene_base,
        "is_fragment": is_fragment,
    }


# ---------------------------------------------------------------------------
# Summary TSV parsing
# ---------------------------------------------------------------------------

def parse_summary(summary_path):
    """Parse the assembly summary TSV and return a per-contig taxonomy lookup.

    Rows where Index == 'NA' (no assembly) are skipped.

    Returns
    -------
    contig_class : dict  contig_id -> class_value
        Maps the full contig ID (e.g. "UK001-A02_contig5") to its Class
        field value (e.g. "Insecta", "no-hit").
    """
    contig_class = {}

    with open(summary_path, "r") as fh:
        col_names = [c.strip() for c in fh.readline().split("\t")]

        try:
            idx_index  = col_names.index("Index")
            idx_contig = col_names.index("Contig")
            idx_class  = col_names.index("Class")
        except ValueError as e:
            print(
                f"ERROR: Summary file is missing expected column: {e}\n"
                f"Expected columns: Index, Contig, Class",
                file=sys.stderr,
            )
            sys.exit(1)

        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")

            while len(fields) <= max(idx_index, idx_contig, idx_class):
                fields.append("")

            index_val = fields[idx_index].strip()
            contig    = fields[idx_contig].strip()
            cls       = fields[idx_class].strip()

            if index_val.upper() == "NA" or not contig or contig.upper() == "NA":
                continue

            contig_class[contig] = cls

    return contig_class


# ---------------------------------------------------------------------------
# Filtering logic
# ---------------------------------------------------------------------------

def apply_annotation_filters(sequences, contig_class):
    """Apply fragmentation, multi-contig duplicate, and taxonomy filters.

    Filter order:

      1 — Fragmented annotations: drop headers where gene_field ends with _N.

      2 — Multi-contig duplicates: runs on the length-filtered pool, so a
          short copy on one contig is already gone. If a sample still has the
          same gene on >1 contig, all copies are genuinely ambiguous and dropped.

      3 — Taxonomy: the source contig must have Class == 'Insecta' in the
          summary TSV. Only the Class column is checked; higher-rank no-hits
          do not cause exclusion. Sequences from contigs not present in the
          summary are passed through with a stderr warning.

    Parameters
    ----------
    sequences    : list of (header, seq)  — already length-filtered
    contig_class : dict  contig_id -> class_value

    Returns
    -------
    passing  : list of (header, seq)
    counts   : dict  fragment / multi_contig / non_insecta
    drop_log : list of dicts  header / sample / contig / gene_field /
               seq_len / reason / class_value
               reason is one of: 'fragment' | 'multi_contig' | 'non_insecta'
               class_value is populated for 'non_insecta' entries only
    """
    dropped_fragment    = 0
    dropped_multicontig = 0
    dropped_non_insecta = 0
    drop_log            = []

    # ---- Filter 1: fragmented annotations --------------------------------
    non_fragment = []
    for header, seq in sequences:
        parsed = parse_header(header)
        if parsed is None:
            non_fragment.append((header, seq, parsed))
            continue
        if parsed["is_fragment"]:
            dropped_fragment += 1
            drop_log.append({
                "header":      header,
                "sample":      parsed["sample"],
                "contig":      parsed["contig"],
                "gene_field":  parsed["gene_field"],
                "seq_len":     len(seq),
                "reason":      "fragment",
                "class_value": "",
            })
        else:
            non_fragment.append((header, seq, parsed))

    # ---- Filter 2: FASTA-level multi-contig per sample -------------------
    sample_contigs = defaultdict(set)
    for header, seq, parsed in non_fragment:
        if parsed is not None:
            sample_contigs[parsed["sample"]].add(parsed["contig"])

    multi_contig_samples = {
        s for s, ctgs in sample_contigs.items() if len(ctgs) > 1
    }

    after_dup_check = []
    for header, seq, parsed in non_fragment:
        if parsed is not None and parsed["sample"] in multi_contig_samples:
            dropped_multicontig += 1
            drop_log.append({
                "header":      header,
                "sample":      parsed["sample"],
                "contig":      parsed["contig"],
                "gene_field":  parsed["gene_field"],
                "seq_len":     len(seq),
                "reason":      "multi_contig",
                "class_value": "",
            })
        else:
            after_dup_check.append((header, seq, parsed))

    # ---- Filter 3: taxonomy — Class must be 'Insecta' --------------------
    passing = []
    for header, seq, parsed in after_dup_check:
        if parsed is None:
            passing.append((header, seq))
            continue

        contig_id = parsed["sample_id"]
        cls = contig_class.get(contig_id)

        if cls is None:
            print(
                f"  WARNING: contig '{contig_id}' not found in summary; "
                f"taxonomy filter skipped for this sequence.",
                file=sys.stderr,
            )
            passing.append((header, seq))
            continue

        if cls == "Insecta":
            passing.append((header, seq))
        else:
            dropped_non_insecta += 1
            drop_log.append({
                "header":      header,
                "sample":      parsed["sample"],
                "contig":      parsed["contig"],
                "gene_field":  parsed["gene_field"],
                "seq_len":     len(seq),
                "reason":      "non_insecta",
                "class_value": cls,
            })

    counts = {
        "fragment":     dropped_fragment,
        "multi_contig": dropped_multicontig,
        "non_insecta":  dropped_non_insecta,
    }
    return passing, counts, drop_log


def find_fasta_file(input_dir, gene_lower):
    """Locate a FASTA file for the given gene name (case-insensitive)."""
    for fname in os.listdir(input_dir):
        name, ext = os.path.splitext(fname)
        if ext.lower() in (".fasta", ".fa", ".fna") and name.lower() == gene_lower:
            return os.path.join(input_dir, fname)
    return None


def filter_fasta(filepath, min_length, contig_class):
    """Apply all filters to a single gene FASTA file.

    Order of operations:
      1. Minimum length  (>= min_length)
      2. Non-fragmented  (no trailing _N on gene field)
      3. Multi-contig duplicate  (FASTA-level, post-length-filter)
      4. Taxonomy        (Class == 'Insecta' in the summary TSV)

    Returns
    -------
    passing  : list of (header, seq)
    counts   : dict  total / too_short / fragment / multi_contig /
                     non_insecta / passing
    drop_log : list of drop-detail dicts for all filters
    """
    all_seqs = list(read_fasta(filepath))
    total    = len(all_seqs)

    # Filter 1: minimum length
    length_passing = []
    too_short_log  = []
    for h, s in all_seqs:
        if len(s) >= min_length:
            length_passing.append((h, s))
        else:
            parsed = parse_header(h)
            too_short_log.append({
                "header":      h,
                "sample":      parsed["sample"]     if parsed else "?",
                "contig":      parsed["contig"]     if parsed else "?",
                "gene_field":  parsed["gene_field"] if parsed else "?",
                "seq_len":     len(s),
                "reason":      "too_short",
                "class_value": "",
            })

    # Filters 2–4
    passing, annotation_counts, annotation_drop_log = apply_annotation_filters(
        length_passing, contig_class
    )

    counts = {
        "total":        total,
        "too_short":    len(too_short_log),
        "fragment":     annotation_counts["fragment"],
        "multi_contig": annotation_counts["multi_contig"],
        "non_insecta":  annotation_counts["non_insecta"],
        "passing":      len(passing),
    }
    return passing, counts, too_short_log + annotation_drop_log


# ---------------------------------------------------------------------------
# Log formatting helpers
# ---------------------------------------------------------------------------

def format_drop_section(drop_log, reason, heading):
    """Return formatted lines for one drop-reason section of the log.

    For 'too_short' entries the observed sequence length is appended.
    For 'non_insecta' entries the actual Class value is appended.
    """
    entries = [d for d in drop_log if d["reason"] == reason]
    if not entries:
        return [f"  {heading}: none"]

    lines = [f"  {heading} ({len(entries)} sequences):"]
    by_sample = defaultdict(list)
    for d in entries:
        by_sample[d["sample"]].append(d)

    for sample in sorted(by_sample):
        for d in sorted(by_sample[sample], key=lambda x: (x["contig"], x["gene_field"])):
            if reason == "too_short":
                note = f"  [{d['seq_len']} bp]"
            elif reason == "non_insecta":
                note = f"  [Class: {d['class_value']}]"
            else:
                note = ""
            lines.append(
                f"    {sample:<20}  {d['contig']:<12}  {d['gene_field']:<12}{note}"
            )
    return lines


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    input_dir = os.path.abspath(args.input)
    if not os.path.isdir(input_dir):
        print(f"ERROR: Input directory not found: {input_dir}", file=sys.stderr)
        sys.exit(1)

    output_dir = (
        os.path.abspath(args.output) if args.output
        else os.path.join(input_dir, "filtered")
    )
    os.makedirs(output_dir, exist_ok=True)

    requested_genes = [g.strip() for g in args.genes.split(",") if g.strip()]
    if not requested_genes:
        print("ERROR: No genes specified via --genes.", file=sys.stderr)
        sys.exit(1)

    unknown_genes = [g for g in requested_genes if g.lower() not in INSECT_GENE_LENGTHS]
    if unknown_genes:
        print(
            f"ERROR: Unknown gene(s): {', '.join(unknown_genes)}\n"
            f"Known genes: {', '.join(CANONICAL_NAMES.values())}",
            file=sys.stderr,
        )
        sys.exit(1)

    # Parse summary TSV
    summary_path = os.path.abspath(args.summary)
    if not os.path.isfile(summary_path):
        print(f"ERROR: Summary file not found: {summary_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Parsing summary : {summary_path}")
    contig_class = parse_summary(summary_path)
    n_non_insecta_contigs = sum(
        1 for cls in contig_class.values() if cls != "Insecta"
    )
    print(f"  Total contigs in summary      : {len(contig_class)}")
    print(f"  Contigs with Class != Insecta : {n_non_insecta_contigs}")
    print()

    header_line = (
        f"{'Gene':<10} {'Ref len':>8} {'Min len':>8} {'Total':>8} "
        f"{'Too short':>10} {'Fragment':>10} {'Multi-ctg':>10} "
        f"{'Non-ins':>8} {'Passing':>8}  Output file"
    )
    divider = "-" * 110

    summary_rows  = []
    all_drop_logs = {}

    for gene in requested_genes:
        gene_lower = gene.lower()
        canonical  = CANONICAL_NAMES[gene_lower]
        ref_length = INSECT_GENE_LENGTHS[gene_lower]
        min_length = ref_length // 2

        fasta_path = find_fasta_file(input_dir, gene_lower)
        if fasta_path is None:
            na = "N/A"
            summary_rows.append(
                f"{canonical:<10} {na:>8} {na:>8} {na:>8} "
                f"{na:>10} {na:>10} {na:>10} {na:>8} {na:>8}  "
                f"[WARNING: no FASTA file found for '{gene}']"
            )
            all_drop_logs[canonical] = []
            continue

        passing_seqs, counts, drop_log = filter_fasta(
            fasta_path, min_length, contig_class
        )
        all_drop_logs[canonical] = drop_log

        out_filename = f"{canonical}_filtered.fasta"
        out_path     = os.path.join(output_dir, out_filename)

        if passing_seqs:
            write_fasta(passing_seqs, out_path)
        else:
            open(out_path, "w").close()

        summary_rows.append(
            f"{canonical:<10} {ref_length:>8} {min_length:>8} {counts['total']:>8} "
            f"{counts['too_short']:>10} {counts['fragment']:>10} "
            f"{counts['multi_contig']:>10} {counts['non_insecta']:>8} "
            f"{counts['passing']:>8}  {out_filename}"
        )

    summary_lines = [header_line, divider] + summary_rows

    print(f"Input directory : {input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Genes requested : {', '.join(requested_genes)}")
    print()
    for line in summary_lines:
        print(line)
    print()
    print("Done.")

    # Write detailed log
    log_path = os.path.join(output_dir, "filter_mitogenes.log")
    with open(log_path, "w") as log:

        log.write("=" * 70 + "\n")
        log.write("FILTER MITO GENES - RUN LOG\n")
        log.write("=" * 70 + "\n\n")
        log.write(f"Input directory : {input_dir}\n")
        log.write(f"Output directory: {output_dir}\n")
        log.write(f"Genes requested : {', '.join(requested_genes)}\n")
        log.write(f"Summary file    : {summary_path}\n")
        log.write(f"Non-Insecta ctg : {n_non_insecta_contigs}\n\n")

        log.write("SUMMARY TABLE\n")
        log.write("-" * 70 + "\n")
        for line in summary_lines:
            log.write(line + "\n")
        log.write("\n")

        log.write("=" * 70 + "\n")
        log.write("DROPPED SEQUENCES - DETAILED BREAKDOWN\n")
        log.write("=" * 70 + "\n\n")
        log.write(
            "Columns: sample  contig  gene_field"
            "  [seq_len bp for too_short | Class: value for non_insecta]\n\n"
        )

        for gene in requested_genes:
            gene_lower = gene.lower()
            canonical  = CANONICAL_NAMES[gene_lower]
            drop_log   = all_drop_logs.get(canonical, [])

            log.write(f"{'─' * 50}\n")
            log.write(f"Gene: {canonical}\n")
            log.write(f"{'─' * 50}\n")

            if not drop_log:
                log.write("  No sequences dropped.\n\n")
                continue

            for line in format_drop_section(drop_log, "too_short",    "Too short"):
                log.write(line + "\n")
            for line in format_drop_section(drop_log, "fragment",     "Fragmented"):
                log.write(line + "\n")
            for line in format_drop_section(drop_log, "multi_contig", "Multi-contig"):
                log.write(line + "\n")
            for line in format_drop_section(drop_log, "non_insecta",  "Non-Insecta"):
                log.write(line + "\n")
            log.write("\n")

    print(f"Log written to  : {log_path}")


if __name__ == "__main__":
    main()
