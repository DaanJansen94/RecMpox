#!/usr/bin/env python3
"""
RecMpox: Recombination flagging of mpox sequences (Ia vs Ib).

Identifies likely recombination breakpoints in mpox consensus genomes by
comparing them to reference lineages (e.g. Ia vs Ib) at diagnostic SNP sites.
This tool works best for flagging recombination events between mpox cases from
sustained human outbreaks (e.g. cocirculation of Ia and Ib). It was initially
designed for investigating recombination between sh2023a and sh2024 during
cocirculation of clades Ia and Ib in Kinshasa.

Pipeline: (1) Align input with Squirrel (fixed coordinates). (2) Diagnostic
SNPs = positions where ref Ia and ref Ib differ (SNPs only, no indels). (3) At
each diagnostic SNP, classify consensus as Ia, Ib, or other. (4) Output: list of
diagnostic SNPs + per-genome counts and percentages (pct Ia, pct Ib, pct
other), highlighting regions where allegiance switches and thus flagging likely
recombination breakpoints.
"""

import argparse
import base64
import json
import logging
import os
import re
import shutil
import ssl
import subprocess
import sys
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

try:
    # Optional: used only for tree-based ancestor distances when --phylogeny is enabled.
    from Bio import Phylo  # type: ignore[import-not-found]

    _HAS_BIOPHYLO = True
except ImportError:
    _HAS_BIOPHYLO = False

from ._version import __version__

logger = logging.getLogger(__name__)

NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# Bundled references for -phylogeny (lives in recmpox/references/ inside the package)
_REFERENCES_DIR = Path(__file__).resolve().parent / "references"
PHYLOGENY_REFS_FASTA = _REFERENCES_DIR / "mpox_references.fasta"
ROOT_TREE_FIGURE_R = _REFERENCES_DIR / "root_tree_figure.R"


LAPIS_MPOX_DETAILS = "https://lapis.pathoplexus.org/mpox/sample/details"
PATHOPLEXUS_FASTA = "https://pathoplexus.org/seq"
MIN_LENGTH_BP = 190_000
PER_GROUP = 5
IA_SH2024_MIN_DATE = "2024-08-19"


def _safe_fasta_id(raw_id: str) -> str:
    """
    Make a FASTA ID safe for external tools (notably Squirrel), which rejects
    some special characters (e.g. ':'). Keep only [A-Za-z0-9_.-], replace the
    rest with '_' and collapse repeats.
    """
    s = (raw_id or "").strip()
    # Prefer the first token (drop long NCBI descriptions)
    s = s.split()[0] if s else "seq"
    s = re.sub(r"[^A-Za-z0-9_.-]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s or "seq"


# === Diagnostic SNP and allegiance logic (inlined from diagnostic_snp.py) ===

PILEUP_SKIP = set(".,$^+-*#")


def _decode_pileup_base(c: str, ref_base: str) -> Optional[str]:
    """Decode one pileup character to base (A/C/G/T), '-' for deletion, or None."""
    c = c.upper()
    if c == "." or c == ",":
        return ref_base.upper()
    if c == "*":
        return "-"  # deletion in read relative to ref
    if c in "ACGTN":
        return c if c != "N" else None
    return None


def run_minimap2_bam(ref_fa: Path, query_fa: Path, out_bam: Path, threads: int = 1) -> bool:
    """Align query to ref with minimap2, output sorted BAM. Returns True on success."""
    out_sam = out_bam.with_suffix(".sam")
    cmd_align = [
        "minimap2",
        "-a",
        "-x",
        "asm5",
        "-t",
        str(threads),
        str(ref_fa),
        str(query_fa),
    ]
    try:
        with open(out_sam, "w") as f:
            subprocess.run(cmd_align, check=True, stdout=f, stderr=subprocess.PIPE, text=True, timeout=600)
        subprocess.run(
            ["samtools", "view", "-b", "-o", str(out_bam.with_suffix(".bam.tmp")), str(out_sam)],
            check=True,
            capture_output=True,
            text=True,
            timeout=60,
        )
        subprocess.run(
            ["samtools", "sort", "-o", str(out_bam), str(out_bam.with_suffix(".bam.tmp"))],
            check=True,
            capture_output=True,
            text=True,
            timeout=120,
        )
        subprocess.run(["samtools", "index", str(out_bam)], check=True, capture_output=True, text=True, timeout=30)
        out_sam.unlink(missing_ok=True)
        Path(out_bam.with_suffix(".bam.tmp")).unlink(missing_ok=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logger.error("minimap2/samtools failed: %s", e)
        return False


def load_ref_sequence(ref_fa: Path) -> Tuple[str, str]:
    """Load first sequence from FASTA as (seq_id, sequence)."""
    name: Optional[str] = None
    seq: List[str] = []
    with open(ref_fa) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    return (name, "".join(seq))
                name = line[1:].split()[0].strip()
                seq = []
            else:
                seq.append(line.strip())
    if name is not None:
        return (name, "".join(seq))
    raise ValueError(f"No sequence in {ref_fa}")


def get_bases_from_mpileup(
    bam: Path,
    ref_fa: Path,
    ref_length: Optional[int] = None,
    include_gaps: bool = True,
) -> Dict[int, str]:
    """
    Run samtools mpileup -f ref bam; return dict ref_pos_1based -> majority base.
    Base is A/C/G/T or '-' (deletion in read relative to ref, from '*' in pileup).
    If ref_length is set, positions with no coverage are filled with '-' (deletion).
    """
    cmd = ["samtools", "mpileup", "-f", str(ref_fa), str(bam), "-Q", "0"]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if result.returncode != 0:
        logger.warning("samtools mpileup failed: %s", result.stderr)
        return {}
    pos_to_bases: Dict[int, List[str]] = {}
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        _, pos_s, ref_base, _, bases = parts[0], parts[1], parts[2], parts[3], parts[4]
        pos = int(pos_s)
        decoded: List[str] = []
        i = 0
        while i < len(bases):
            c = bases[i]
            if c in "^":
                i += 2
                continue
            if c in "$":
                i += 1
                continue
            if c in "+-":
                i += 1
                n = ""
                while i < len(bases) and bases[i].isdigit():
                    n += bases[i]
                    i += 1
                num = int(n) if n else 0
                i += num
                continue
            b = _decode_pileup_base(c, ref_base)
            if b is not None:
                decoded.append(b)
            i += 1
        if not decoded:
            continue
        from collections import Counter

        (base, _) = Counter(decoded).most_common(1)[0]
        pos_to_bases[pos] = base
    if ref_length is not None and include_gaps:
        for pos in range(1, ref_length + 1):
            if pos not in pos_to_bases:
                pos_to_bases[pos] = "-"
    return pos_to_bases


def build_diagnostic_snps_from_alignment(seq_ia: str, seq_ib: str) -> List[Tuple[int, str, str]]:
    """
    Diagnostic SNPs from two aligned sequences (e.g. from Squirrel).
    Position is diagnostic only if both have a base (A/C/G/T) and they differ (SNPs only).
    Returns list of (position_1based, ia_allele, ib_allele).
    """
    if len(seq_ia) != len(seq_ib):
        logger.warning("Aligned Ia length %d != Ib length %d; using min", len(seq_ia), len(seq_ib))
    n = min(len(seq_ia), len(seq_ib))
    diagnostic: List[Tuple[int, str, str]] = []
    both_base = 0
    same_base = 0
    for i in range(n):
        ia_b = seq_ia[i].upper()
        ib_b = seq_ib[i].upper()
        if ia_b not in "ACGT" or ib_b not in "ACGT":
            continue
        both_base += 1
        if ia_b == ib_b:
            same_base += 1
            continue
        diagnostic.append((i + 1, ia_b, ib_b))
    logger.info(
        "Alignment length: %d bp; columns with both refs ACGT: %d; same base: %d; diagnostic SNPs: %d",
        n,
        both_base,
        same_base,
        len(diagnostic),
    )
    return diagnostic


def find_large_indels(
    seq_ia: str,
    seq_ib: str,
    min_size: int = 100,
) -> List[Tuple[int, int, str]]:
    """
    Find runs where one ref has bases and the other has gaps (indels), length >= min_size bp.
    Both directions: (A) ref1 has bases, ref2 has gap → "ia" (deletion in ref2/Ib);
    (B) ref1 has gap, ref2 has bases → "ib" (deletion in ref1/Ia).
    Returns list of (start_1based, end_1based, ref_with_bases).
    Consensus with sequence (no deletion) → count ref_with_bases; consensus with gap (has deletion) → count the other ref.
    """
    if len(seq_ia) != len(seq_ib):
        n = min(len(seq_ia), len(seq_ib))
    else:
        n = len(seq_ia)
    indels: List[Tuple[int, int, str]] = []
    i = 0
    while i < n:
        ia_b = seq_ia[i].upper()
        ib_b = seq_ib[i].upper()
        ia_has_base = ia_b in "ACGT"
        ib_has_base = ib_b in "ACGT"
        ia_has_gap = ia_b in "-"
        ib_has_gap = ib_b in "-"
        if ia_has_base and ib_has_gap:
            start = i
            while i < n and seq_ia[i].upper() in "ACGT" and seq_ib[i].upper() in "-":
                i += 1
            run_len = i - start
            if run_len >= min_size:
                indels.append((start + 1, i, "ia"))
            continue
        if ia_has_gap and ib_has_base:
            start = i
            while i < n and seq_ia[i].upper() in "-" and seq_ib[i].upper() in "ACGT":
                i += 1
            run_len = i - start
            if run_len >= min_size:
                indels.append((start + 1, i, "ib"))
            continue
        i += 1
    logger.info("Found %d large indels (>= %d bp) between ref1 and ref2", len(indels), min_size)
    return indels


def classify_query_at_indel(query_seq: str, start_1based: int, end_1based: int, ref_with_bases: str) -> str:
    """
    Classify query at an indel region (whole region, majority rule). Kept for compatibility.
    ref_with_bases is "ia" or "ib". Query mostly base → ref_with_bases; mostly gap → other ref.
    """
    start_idx = start_1based - 1
    end_idx = end_1based
    if start_idx < 0 or end_idx > len(query_seq):
        return "ambiguous"
    segment = query_seq[start_idx:end_idx]
    n_base = sum(1 for c in segment if c.upper() in "ACGT")
    n_gap = sum(1 for c in segment if c.upper() in "-")
    n_other = len(segment) - n_base - n_gap
    if n_other > len(segment) / 2:
        return "ambiguous"
    if n_base > n_gap:
        return ref_with_bases
    if n_gap > n_base:
        return "ib" if ref_with_bases == "ia" else "ia"
    return "ambiguous"


def classify_query_at_indel_column(query_seq: str, pos_1based: int, ref_with_bases: str) -> str:
    """
    Classify query at a single column in an indel: ref_with_bases has the base, the other ref has the gap.
    Query has base (ACGT) → ref_with_bases; query has gap (-) or N, ., * (missing/poor data) → the other ref (has the deletion).
    So N and - at indel columns are classified as the ref that has the deletion; when using -include-indels,
    good coverage of those regions is recommended for reliable classification.
    """
    idx = pos_1based - 1
    if idx < 0 or idx >= len(query_seq):
        return "ambiguous"
    q = query_seq[idx].upper()
    if q in "ACGT":
        return ref_with_bases
    if q in "-N.*":
        return "ib" if ref_with_bases == "ia" else "ia"
    return "ambiguous"


def build_diagnostic_snps(
    ref_ia_fa: Path,
    ref_ib_bam: Path,
    ref_ia_seq: Optional[Tuple[str, str]] = None,
) -> List[Tuple[int, str, str]]:
    """
    Diagnostic sites = positions where ref_ia and ref_ib both have a base (A/C/G/T) and differ (SNPs only).
    Excludes insertions and deletions (only SNPs count).
    Returns list of (position_1based, ia_allele, ib_allele).
    """
    if ref_ia_seq is None:
        ref_ia_seq = load_ref_sequence(ref_ia_fa)
    _, ref_ia_sequence = ref_ia_seq
    ref_len = len(ref_ia_sequence)
    ref_ib_bases = get_bases_from_mpileup(ref_ib_bam, ref_ia_fa, ref_length=ref_len, include_gaps=True)
    diagnostic: List[Tuple[int, str, str]] = []
    for pos_1based in range(1, ref_len + 1):
        idx = pos_1based - 1
        ia_base = ref_ia_sequence[idx].upper()
        if ia_base not in "ACGT":
            continue
        ib_base = ref_ib_bases.get(pos_1based, "-").strip()
        if ib_base not in "ACGT-*":
            ib_base = "-"
        if ib_base == "*":
            ib_base = "-"
        if ib_base not in "ACGT":
            continue
        if ia_base == ib_base:
            continue
        diagnostic.append((pos_1based, ia_base, ib_base))
    logger.info("Found %d diagnostic SNPs (Ia vs Ib, SNPs only; indels excluded)", len(diagnostic))
    return diagnostic


def load_alignment_fasta(aln_fasta: Path) -> Dict[str, str]:
    """
    Load alignment FASTA: dict seq_id -> sequence (no gaps stripped).
    All sequences must have the same length.
    """
    seqs: Dict[str, List[str]] = {}
    current_id: Optional[str] = None
    with open(aln_fasta) as f:
        for line in f:
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(seqs[current_id])
                current_id = line[1:].split()[0].strip().replace("/", "_")
                seqs[current_id] = []
            else:
                if current_id is not None:
                    seqs[current_id].append(line.strip())
        if current_id is not None:
            seqs[current_id] = "".join(seqs[current_id])
    return seqs


def _load_tree_and_distances(
    tree_path: Path,
) -> Optional[Dict[Tuple[str, str], float]]:
    """
    Load Newick tree and return symmetric patristic distances between all tip pairs.

    Distances are keyed by (a, b) with a < b. Requires Bio.Phylo; returns None
    if Biopython is not available.
    """
    if not _HAS_BIOPHYLO:
        return None
    try:
        tree = Phylo.read(str(tree_path), "newick")  # type: ignore[name-defined]
    except Exception as e:  # pragma: no cover - defensive
        logger.warning("--phylogeny: could not read treefile %s for ancestor labels: %s", tree_path, e)
        return None
    # Midpoint-root if supported; if it fails, continue with original tree.
    if hasattr(tree, "root_at_midpoint"):
        try:
            tree.root_at_midpoint()
        except Exception:  # pragma: no cover - best-effort
            pass
    terminals = list(tree.get_terminals())
    name_to_clade: Dict[str, Any] = {}
    for clade in terminals:
        if getattr(clade, "name", None):
            name_to_clade[clade.name.strip()] = clade
    if not name_to_clade:
        logger.warning("--phylogeny: no tip names found in treefile %s", tree_path)
        return None
    dists: Dict[Tuple[str, str], float] = {}
    names = list(name_to_clade.keys())
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = names[i], names[j]
            try:
                d = float(tree.distance(name_to_clade[a], name_to_clade[b]))
            except Exception:  # pragma: no cover - defensive
                continue
            dists[(a, b)] = d
    return dists


def _patristic_dist(a: str, b: str, dists: Dict[Tuple[str, str], float]) -> Optional[float]:
    """Lookup symmetric patristic distance between two tip names."""
    if a == b:
        return 0.0
    key = (a, b) if a < b else (b, a)
    return dists.get(key)


def get_query_allegiance_from_alignment(
    query_seq: str,
    diagnostic_snps: List[Tuple[int, str, str]],
    ref_length: int,
    diagnostic_indels: Optional[List[Tuple[int, int, str]]] = None,
) -> List[Tuple[int, str]]:
    """
    At each diagnostic (pos, ia_allele, ib_allele), get query base from alignment column (1-based).
    N -> 'other_n'; other non-ACGT/gap -> 'ambiguous'. Same allegiance rules as get_query_allegiance.
    If diagnostic_indels is provided, append one allegiance per column in each indel (each position counts as a site).
    """
    if len(query_seq) < ref_length:
        logger.warning("Alignment sequence length %d < ref length %d", len(query_seq), ref_length)
    result: List[Tuple[int, str]] = []
    for (pos, ia_a, ib_a) in diagnostic_snps:
        idx = pos - 1
        if idx >= len(query_seq):
            result.append((pos, "ambiguous"))
            continue
        q = query_seq[idx].upper()
        if q == "N":
            result.append((pos, "other_n"))
            continue
        if q not in "ACGT-":
            result.append((pos, "ambiguous"))
            continue
        if ia_a in "ACGT" and ib_a in "ACGT":
            if q in "ACGT":
                if q == ia_a:
                    result.append((pos, "ia"))
                elif q == ib_a:
                    result.append((pos, "ib"))
                else:
                    result.append((pos, "ambiguous"))
            else:
                result.append((pos, "ambiguous"))
        elif ia_a in "ACGT" and ib_a == "-":
            if q in "ACGT":
                result.append((pos, "ia"))
            elif q == "-":
                result.append((pos, "ib"))
            else:
                result.append((pos, "ambiguous"))
        elif ia_a == "-" and ib_a in "ACGT":
            if q == "-":
                result.append((pos, "ia"))
            elif q in "ACGT" and q == ib_a:
                result.append((pos, "ib"))
            else:
                result.append((pos, "ambiguous"))
        else:
            result.append((pos, "ambiguous"))
    if diagnostic_indels:
        for (start, end, ref_who) in diagnostic_indels:
            for pos in range(start, end + 1):
                a = classify_query_at_indel_column(query_seq, pos, ref_who)
                result.append((pos, a))
    return result


def get_query_allegiance(
    query_bam: Path,
    ref_ia_fa: Path,
    diagnostic_snps: List[Tuple[int, str, str]],
    ref_length: Optional[int] = None,
) -> List[Tuple[int, str]]:
    """
    At each diagnostic (pos, ia_allele, ib_allele), get query base/gap from pileup.
    Allegiance: 'ia' if query matches Ia, 'ib' if matches Ib, 'ambiguous' else.
    Handles gaps: if Ib has deletion (ib_allele='-'), query base -> ia, query gap -> ib.
    """
    query_bases = get_bases_from_mpileup(query_bam, ref_ia_fa, ref_length=ref_length, include_gaps=True)
    result: List[Tuple[int, str]] = []
    for (pos, ia_a, ib_a) in diagnostic_snps:
        q = query_bases.get(pos, "-").strip()
        if q not in "ACGT-":
            q = "-"
        if q == "*":
            q = "-"
        if ia_a in "ACGT" and ib_a in "ACGT":
            if q in "ACGT":
                if q == ia_a:
                    result.append((pos, "ia"))
                elif q == ib_a:
                    result.append((pos, "ib"))
                else:
                    result.append((pos, "ambiguous"))
            else:
                result.append((pos, "ambiguous"))
        elif ia_a in "ACGT" and ib_a == "-":
            if q in "ACGT":
                result.append((pos, "ia"))
            elif q == "-":
                result.append((pos, "ib"))
            else:
                result.append((pos, "ambiguous"))
        elif ia_a == "-" and ib_a in "ACGT":
            if q == "-":
                result.append((pos, "ia"))
            elif q in "ACGT" and q == ib_a:
                result.append((pos, "ib"))
            else:
                result.append((pos, "ambiguous"))
        else:
            result.append((pos, "ambiguous"))
    return result


def allegiance_summary(positions_allegiances: List[Tuple[int, str]]) -> Tuple[int, int, int, int]:
    """Return (n_ia, n_ib, n_other, n_other_n) over all diagnostic sites. n_other = ambiguous + other_n."""
    n_ia = sum(1 for _, a in positions_allegiances if a == "ia")
    n_ib = sum(1 for _, a in positions_allegiances if a == "ib")
    n_amb = sum(1 for _, a in positions_allegiances if a == "ambiguous")
    n_other_n = sum(1 for _, a in positions_allegiances if a == "other_n")
    n_other = n_amb + n_other_n
    return (n_ia, n_ib, n_other, n_other_n)


def allegiance_summary_snp_only(
    positions_allegiances: List[Tuple[int, str]],
    diagnostic_snp_positions: List[int],
) -> Tuple[int, int, int, int]:
    """
    Return (n_ia, n_ib, n_other, n_other_n) over diagnostic SNP positions only (excludes indel columns).
    Use this to classify consensus as Ia/Ib from SNP percentages so that poor coverage
    in deletion regions (often 'N') does not inflate 'other'.
    """
    snp_positions = set(diagnostic_snp_positions)
    snp_only = [(p, a) for (p, a) in positions_allegiances if p in snp_positions]
    return allegiance_summary(snp_only)


def consensus_from_snp_percentages(
    n_ia_snp: int,
    n_ib_snp: int,
    n_other_snp: int,
    ref1_label: str,
    ref2_label: str,
    pct_threshold: float = 10.0,
) -> str:
    """
    Classify consensus from diagnostic-SNP-only percentages: if >pct_threshold of SNPs
    are ref1 → ref1; if >pct_threshold are ref2 → ref2; else 'other'.
    Tie-break: if both above threshold, assign to the higher percentage.
    Used for Ia vs Ib (Ia = no deletion, Ib = deletion present) and Ib vs IIb.
    """
    total = n_ia_snp + n_ib_snp + n_other_snp
    if total <= 0:
        return "other"
    pct_ia = 100.0 * n_ia_snp / total
    pct_ib = 100.0 * n_ib_snp / total
    if pct_ia > pct_threshold and pct_ia >= pct_ib:
        return ref1_label
    if pct_ib > pct_threshold and pct_ib > pct_ia:
        return ref2_label
    return "other"


def get_runs_and_breakpoints(
    positions_allegiances: List[Tuple[int, str]],
    diagnostic_snp_positions: List[int],
    min_consecutive: int = 1,
    ignore_other: bool = True,
) -> Tuple[List[Tuple[int, int, str, int]], List[Tuple[int, int, str, str]]]:
    """
    Build runs of consecutive ia/ib along diagnostic SNPs.
    When ignore_other=True (default): "other" positions (N, gap, ambiguous base) are
    transparent to the run builder — a tract continues through them and only ends when
    an actual opposing-clade SNP is encountered.  end_pos and n_snps reflect only the
    clade-matching positions; "other" positions inside a tract are excluded from its
    endpoints and SNP count.
    A breakpoint is only called when *both* flanking runs have >= min_consecutive SNPs.
    Returns:
        runs: list of (start_pos, end_pos, clade, n_snps) with clade in ("ia", "ib").
        breakpoints: list of (end_pos_before, start_pos_after, clade_before, clade_after).
    """
    snp_positions = set(diagnostic_snp_positions)
    ordered = [(p, a) for (p, a) in positions_allegiances if p in snp_positions]
    ordered.sort(key=lambda x: x[0])
    runs: List[Tuple[int, int, str, int]] = []
    i = 0
    while i < len(ordered):
        pos, a = ordered[i]
        if a not in ("ia", "ib"):
            i += 1
            continue
        current_clade = a
        start_pos = pos
        last_clade_pos = pos
        n_snps = 1
        i += 1
        while i < len(ordered):
            next_pos, next_a = ordered[i]
            if next_a == current_clade:
                last_clade_pos = next_pos
                n_snps += 1
                i += 1
            elif next_a in ("other", "ambiguous") and ignore_other:
                i += 1
            else:
                break
        runs.append((start_pos, last_clade_pos, current_clade, n_snps))
    breakpoints: List[Tuple[int, int, str, str]] = []
    for j in range(len(runs) - 1):
        _, end_a, clade_a, n_a = runs[j]
        _, start_b, clade_b, n_b = runs[j + 1]
        if clade_a == clade_b:
            continue
        if n_a >= min_consecutive and n_b >= min_consecutive:
            breakpoints.append((end_a, start_b, clade_a, clade_b))
    return (runs, breakpoints)


def rolling_snp_window(
    positions_allegiances: List[Tuple[int, str]],
    window_size: int,
    informative_only: bool = True,
    min_informative: int = 5,
) -> List[Tuple[int, float]]:
    """
    Rolling window of N SNPs. At each step, fraction of SNPs in window that are 'ia'.
    If informative_only: ia_frac = n_ia / (n_ia + n_ib), excluding ambiguous; window
    with fewer than min_informative (ia+ib) SNPs gets previous frac or 0.5 (skip effect).
    Returns list of (position_of_last_snp_in_window, ia_fraction).
    """
    if window_size <= 0 or len(positions_allegiances) < window_size:
        return []
    out: List[Tuple[int, float]] = []
    for i in range(window_size - 1, len(positions_allegiances)):
        window = positions_allegiances[i - window_size + 1 : i + 1]
        n_ia = sum(1 for _, a in window if a == "ia")
        n_ib = sum(1 for _, a in window if a == "ib")
        pos = window[-1][0]
        if informative_only:
            n_inf = n_ia + n_ib
            if n_inf < min_informative:
                ia_frac = 0.5
            else:
                ia_frac = n_ia / n_inf
        else:
            ia_frac = n_ia / len(window)
        out.append((pos, ia_frac))
    return out


def flag_mixed_regions(
    positions_ia_frac: List[Tuple[int, float]],
    min_mixed: float = 0.3,
    max_mixed: float = 0.7,
) -> List[Tuple[int, int]]:
    """Merge consecutive positions where min_mixed <= ia_fraction <= max_mixed into regions (start, end)."""
    if not positions_ia_frac:
        return []
    regions: List[Tuple[int, int]] = []
    in_region = False
    start: Optional[int] = None
    for pos, frac in positions_ia_frac:
        if min_mixed <= frac <= max_mixed:
            if not in_region:
                start = pos
                in_region = True
        else:
            if in_region:
                regions.append((start, pos))
                in_region = False
    if in_region and start is not None:
        regions.append((start, positions_ia_frac[-1][0]))
    return regions


def plot_allegiance(
    positions_ia_frac: List[Tuple[int, float]],
    flagged_regions: List[Tuple[int, int]],
    out_path: Path,
    title: str = "Clade Ia allegiance (rolling SNP window)",
) -> None:
    """Plot position (x) vs ia_fraction (y); shade flagged regions."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not installed; skipping plot")
        return
    if not positions_ia_frac:
        return
    xs = [p[0] for p in positions_ia_frac]
    ys = [p[1] for p in positions_ia_frac]
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(xs, ys, color="black", linewidth=0.8, alpha=0.9)
    for (start, end) in flagged_regions:
        ax.axvspan(start, end, alpha=0.25, color="red")
    ax.axhline(0.5, color="gray", linestyle="--", alpha=0.5)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("Position (reference)")
    ax.set_ylabel("Proportion Ia (rolling window)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    logger.info("Wrote plot %s", out_path)


def _fetch_url_lapis(url: str, params: str = "") -> str:
    full = f"{url}?{params}" if params else url
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    req = urllib.request.Request(full, headers={"User-Agent": "RecMpox/1.0"})
    with urllib.request.urlopen(req, timeout=90, context=ctx) as resp:
        return resp.read().decode("utf-8")


def _lapis_fetch(q: Dict[str, Any], limit: int = 5000) -> List[Dict[str, Any]]:
    """Fetch from LAPIS details; return list of dicts (raw rows)."""
    q = dict(q)
    q["limit"] = limit
    params = urllib.parse.urlencode(q)
    try:
        data = _fetch_url_lapis(LAPIS_MPOX_DETAILS, params)
        obj = json.loads(data)
        data_list = obj.get("data") or []
        return data_list
    except Exception as e:
        logger.warning("Pathoplexus LAPIS failed: %s", e)
        return []


def _lapis_row_to_tuple(r: Dict[str, Any], date_key_priority: Optional[List[str]] = None) -> Optional[Tuple[str, str, int, Optional[str]]]:
    """(accession_version, date_sort_key, length, insdc) or None."""
    acc_ver = (r.get("accessionVersion") or r.get("accession") or "").strip()
    length = r.get("length")
    if not acc_ver or length is None or int(length) < MIN_LENGTH_BP:
        return None
    insdc = (r.get("insdcAccessionFull") or "").strip() or None
    date = None
    for key in date_key_priority or ["sampleCollectionDate", "sampleCollectionDateRangeLower", "sampleCollectionDateRangeUpper"]:
        v = r.get(key)
        if isinstance(v, str) and v.strip():
            date = v.strip()
            break
    if not date:
        date = "9999-99-99"
    return (acc_ver, date, int(length), insdc)


def _fetch_ib_kinshasa(limit: int = 5000) -> List[Tuple[str, str, int, Optional[str]]]:
    rows = _lapis_fetch(
        {
            "geoLocCountry": "Democratic Republic of the Congo",
            "geoLocAdmin1": "Kinshasa",
            "clade": "Ib",
            "lengthFrom": MIN_LENGTH_BP,
        },
        limit=limit,
    )
    out: List[Tuple[str, str, int, Optional[str]]] = []
    for r in rows:
        t = _lapis_row_to_tuple(r, ["sampleCollectionDate"])
        if t and t[1] != "9999-99-99":
            out.append(t)
    out.sort(key=lambda x: (x[1], x[0]))
    return out[:PER_GROUP]


def _fetch_ia_kinshasa_sh2024(limit: int = 5000) -> List[Tuple[str, str, int, Optional[str]]]:
    rows = _lapis_fetch(
        {
            "geoLocCountry": "Democratic Republic of the Congo",
            "geoLocAdmin1": "Kinshasa",
            "clade": "Ia",
            "outbreak": "sh2024",
            "lengthFrom": MIN_LENGTH_BP,
        },
        limit=limit,
    )
    out: List[Tuple[str, str, int, Optional[str]]] = []
    for r in rows:
        t = _lapis_row_to_tuple(r, ["sampleCollectionDate"])
        if t and t[1] != "9999-99-99" and t[1] >= IA_SH2024_MIN_DATE:
            out.append(t)
    out.sort(key=lambda x: (x[1], x[0]))
    return out[:PER_GROUP]


def _fetch_sh2017(limit: int = 5000) -> List[Tuple[str, str, int, Optional[str]]]:
    rows = _lapis_fetch(
        {
            "outbreak": "sh2017",
            "lengthFrom": MIN_LENGTH_BP,
        },
        limit=limit,
    )
    out: List[Tuple[str, str, int, Optional[str]]] = []
    for r in rows:
        t = _lapis_row_to_tuple(r, ["sampleCollectionDate"])
        if t and t[1] != "9999-99-99":
            out.append(t)
    out.sort(key=lambda x: (x[1], x[0]))
    return out[:PER_GROUP]


def _fetch_iia_earliest(limit: int = 5000) -> List[Tuple[str, str, int, Optional[str]]]:
    rows = _lapis_fetch(
        {
            "clade": "IIa",
            "lengthFrom": MIN_LENGTH_BP,
        },
        limit=limit,
    )
    seen_base: set = set()
    out: List[Tuple[str, str, int, Optional[str]]] = []
    for r in rows:
        if r.get("versionStatus") != "LATEST_VERSION":
            continue
        t = _lapis_row_to_tuple(
            r,
            ["sampleCollectionDate", "sampleCollectionDateRangeLower", "sampleCollectionDateRangeUpper"],
        )
        if not t:
            continue
        base = (t[0].split(".")[0], t[1])
        if base in seen_base:
            continue
        seen_base.add(base)
        out.append(t)
    out.sort(key=lambda x: (x[1], x[0]))
    return out[:PER_GROUP]


def _fetch_fasta_pathoplexus(accession_version: str, out_path: Path) -> bool:
    url = f"{PATHOPLEXUS_FASTA}/{accession_version}.fa"
    try:
        data = _fetch_url_lapis(url, "")
        if not data.strip() or "not found" in data.lower()[:200]:
            return False
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(data)
        return True
    except Exception:
        return False


def _align_consensus_group(combined_fa: Path, out_dir: Path, stem: str, use_clade_ii: bool) -> Optional[Path]:
    """Run Squirrel (cladei or cladeii) or mafft; return path to alignment file or None."""
    squirrel_out = out_dir / "squirrel_out" / stem
    squirrel_out.mkdir(parents=True, exist_ok=True)
    expected_aln = squirrel_out / (combined_fa.stem + ".aln.fasta")
    clade = "cladeii" if use_clade_ii else "cladei"
    try:
        subprocess.run(
            ["squirrel", "--clade", clade, str(combined_fa), "-o", str(squirrel_out), "--tempdir", str(squirrel_out / "tmp")],
            check=True,
            capture_output=True,
            text=True,
            timeout=600,
        )
    except FileNotFoundError:
        pass
    except subprocess.CalledProcessError:
        pass
    if expected_aln.exists():
        return expected_aln
    aln_files = list(squirrel_out.glob("*.aln.fasta"))
    if aln_files:
        return aln_files[0]
    try:
        with open(expected_aln, "w") as out:
            subprocess.run(
                ["mafft", "--auto", "--quiet", str(combined_fa)],
                check=True,
                stdout=out,
                text=True,
                timeout=600,
            )
        return expected_aln
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None


def _build_consensus_from_aln(aln_path: Path, consensus_stem: str) -> Optional[Tuple[str, int]]:
    """Return (consensus_content_with_header, len_ungapped) or None."""
    seqs: Dict[str, str] = {}
    with open(aln_path) as f:
        current_id: Optional[str] = None
        current_seq: List[str] = []
        for line in f:
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0].strip().replace("/", "_")
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_id is not None:
            seqs[current_id] = "".join(current_seq)
    if not seqs:
        return None
    aln_len = len(next(iter(seqs.values())))
    consensus: List[str] = []
    for col in range(aln_len):
        counts = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0, "-": 0}
        for seq in seqs.values():
            b = seq[col].upper() if col < len(seq) else "-"
            if b in counts:
                counts[b] += 1
            else:
                counts["N"] += 1
        acgt = {k: counts[k] for k in "ACGT"}
        best = max(acgt.items(), key=lambda x: x[1])
        total_acgt = sum(acgt.values())
        if total_acgt == 0:
            consensus.append("N")
        elif best[1] > total_acgt / 2:
            consensus.append(best[0])
        else:
            consensus.append("N")
    consensus_ungapped = "".join(consensus).replace("-", "")
    body = "\n".join(consensus_ungapped[i : i + 80] for i in range(0, len(consensus_ungapped), 80)) + "\n"
    header = f">{consensus_stem}\n"
    return (header + body, len(consensus_ungapped))


def _build_lapis_consensus_refs(clades: List[str], consensus_output_dir: Path, work_dir: Path) -> Dict[str, Path]:
    """
    Fetch earliest 5 genomes per selected clade from LAPIS/Pathoplexus and build one consensus FASTA per clade.
    Returns mapping from clade label (Ia/Ib/IIa/IIb) to FASTA path inside consensus_output_dir.
    """
    consensus_output_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir = work_dir / "earliest_consensus_tmp"
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir, ignore_errors=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    try:
        fasta_dir = tmp_dir / "fasta"
        fasta_dir.mkdir(parents=True, exist_ok=True)

        all_groups = [
            ("Ib_Kinshasa", _fetch_ib_kinshasa, "ib_kinshasa", False, "Ib", "sh2023Ib"),
            ("Ia_Kinshasa_sh2024", _fetch_ia_kinshasa_sh2024, "ia_kinshasa", False, "Ia", "sh2024Ia"),
            ("sh2017", _fetch_sh2017, "sh2017", True, "IIb", "sh2017IIb"),
            ("IIa", _fetch_iia_earliest, "iia", True, "IIa", "iia"),
        ]
        wanted = set(clades)
        groups = [g for g in all_groups if g[4] in wanted]
        if not groups:
            raise RuntimeError(f"No supported clades selected for consensus build: {clades!r}")

        out_paths: Dict[str, Path] = {}
        for name, fetch_fn, consensus_stem, use_clade_ii, clade_label, header_stem in groups:
            logger.info("Consensus refs: fetching %s ...", name)
            rows = fetch_fn()
            if len(rows) < 2:
                logger.warning("Consensus refs: skip %s (need at least 2 samples, got %d)", name, len(rows))
                continue

            group_fastas: List[Path] = []
            for acc_ver, date, length, insdc in rows:
                safe_id = acc_ver.replace(".", "_").replace("/", "_")
                path = fasta_dir / f"{name}_{safe_id}.fa"
                if not _fetch_fasta_pathoplexus(acc_ver, path) and insdc:
                    fetch_nucleotide_fasta(insdc.split(".")[0], path)
                if path.exists():
                    group_fastas.append(path)
            if len(group_fastas) < 2:
                logger.warning("Consensus refs: skip %s (could not download enough FASTAs)", name)
                continue

            combined_fa = tmp_dir / f"samples_combined_{consensus_stem}.fa"
            with open(combined_fa, "w") as out_f:
                for p in sorted(group_fastas):
                    text = p.read_text()
                    out_f.write(text)
                    if text and not text.endswith("\n"):
                        out_f.write("\n")

            aln_path = _align_consensus_group(combined_fa, tmp_dir, consensus_stem, use_clade_ii)
            if not aln_path or not aln_path.exists():
                logger.warning("Consensus refs: skip %s (alignment failed; install squirrel and/or mafft)", name)
                continue

            result = _build_consensus_from_aln(aln_path, header_stem)
            if not result:
                logger.warning("Consensus refs: skip %s (consensus build failed)", name)
                continue
            content, length_bp = result
            lines = content.split("\n")
            if lines and lines[0].startswith(">"):
                lines[0] = f">{header_stem}"
            content = "\n".join(lines) + ("\n" if content.endswith("\n") else "")
            out_fa = consensus_output_dir / f"{header_stem}.fa"
            out_fa.parent.mkdir(parents=True, exist_ok=True)
            out_fa.write_text(content)
            logger.info("Consensus refs: wrote %s (length %d bp)", out_fa, length_bp)
            out_paths[clade_label] = out_fa

        if not out_paths:
            raise RuntimeError("Consensus refs: no consensus files produced")
        return out_paths
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def _sanitize_fasta_ids(in_path: Path, out_path: Path) -> None:
    """Rewrite FASTA with safe IDs; keep sequences unchanged."""
    seen: Dict[str, int] = {}
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(in_path) as inp, open(out_path, "w") as out:
        for line in inp:
            if line.startswith(">"):
                raw = line[1:].strip()
                base = _safe_fasta_id(raw)
                n = seen.get(base, 0) + 1
                seen[base] = n
                safe = base if n == 1 else f"{base}_{n}"
                out.write(f">{safe}\n")
            else:
                out.write(line)


def fetch_nucleotide_fasta(accession: str, out_path: Path) -> bool:
    """Fetch nucleotide by NCBI accession; save as FASTA. Returns True on success.
    SSL certificate verification is disabled to work in restricted networks (proxy/firewall).
    """
    params = f"db=nucleotide&id={urllib.parse.quote(accession.strip())}&rettype=fasta&retmode=text"
    url = f"{NCBI_EFETCH}?{params}"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "RecMpox/1.0"})
        # Disable SSL verification to work in restricted networks (DRC, proxy, firewall)
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
        with urllib.request.urlopen(req, timeout=60, context=ctx) as resp:
            data = resp.read().decode("utf-8")
        if not data.strip() or "Error" in data[:200]:
            return False
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w") as f:
            f.write(data)
        logger.info("Downloaded %s -> %s", accession, out_path)
        return True
    except Exception as e:
        logger.error("Fetch failed for %s: %s", accession, e)
        return False


def resolve_ref(spec: str, work_dir: Path, label: str) -> Optional[Path]:
    """Resolve ref: existing FASTA path or NCBI accession (download). Returns path or None."""
    spec = spec.strip()
    path = Path(spec)
    if path.is_file():
        logger.info("Using reference %s: %s", label, path)
        return path
    safe = re.sub(r"[^\w.-]", "_", spec)
    out_path = work_dir / f"ref_{label}_{safe}.fa"
    if out_path.is_file():
        logger.info("Using cached reference %s: %s", label, out_path)
        return out_path
    if not fetch_nucleotide_fasta(spec, out_path):
        return None
    return out_path


# Minor ref % = smaller of pct_ref1, pct_ref2 over ALL diagnostic sites. When >= threshold, flag as potential recombinant.
# Default recombinant threshold: 10% for all (intra- and inter-clade)
MINOR_REF_PCT_THRESHOLD = 10.0

# Default NCBI accessions per clade/subclade (used when -ref X,Y is given)
REF_DEFAULTS = {
    "Ia": "OZ254474.1",   # sh2024a, near-complete Clade I
    "Ib": "PP601219.1",   # sh2023b, near-complete Clade I
    "IIa": "OZ287284.1",  # Clade IIa
    "IIb": "NC_063383.1", # Clade IIb RefSeq (sh2017)
}
# Human-readable names for default refs (shown in HTML summary)
REF_DEFAULT_NAMES = {
    "OZ254474.1": "sh2024a",
    "PP601219.1": "sh2023b",
    "OZ287284.1": "Clade IIa",
    "NC_063383.1": "sh2017",
}

HTML_CHUNK_SIZE = 100   # when more than this many genomes, split HTML into one file per chunk (overzichtelijk)


def _recombinant_call_minor_pct(n_ref1: int, n_ref2: int, n_total: int, minor_pct_threshold: float) -> str:
    """Return 'potential recombinant' when minor ref % >= threshold. Minor ref % = smaller of pct_ref1, pct_ref2 over ALL diagnostic sites (Ia+Ib+other)."""
    if n_total <= 0:
        return "no recombinant"
    minor_pct = 100.0 * min(n_ref1, n_ref2) / n_total
    if minor_pct >= minor_pct_threshold:
        return "potential recombinant"
    return "no recombinant"


def _infer_squirrel_clade(ref1_g: Optional[str], ref2_g: Optional[str]) -> Optional[str]:
    """
    Infer Squirrel --clade from ref genotype labels (-ref1_g, -ref2_g).
    Both Clade I (Ia, Ib) -> cladei; both Clade II (IIa, IIb) -> cladeii; mix or unclear -> None (run Squirrel without --clade).
    """
    def looks_clade_i(lbl: str) -> bool:
        l = lbl.lower().strip()
        return l in ("ia", "ib") or (l.startswith("i") and not l.startswith("ii"))
    def looks_clade_ii(lbl: str) -> bool:
        l = lbl.lower().strip()
        return l in ("iia", "iib") or l.startswith("ii")
    l1 = (ref1_g or "").strip()
    l2 = (ref2_g or "").strip()
    if not l1 or not l2:
        return None
    if looks_clade_i(l1) and looks_clade_i(l2):
        return "cladei"
    if looks_clade_ii(l1) and looks_clade_ii(l2):
        return "cladeii"
    return None


def _short_ref_label(spec: str, max_len: int = 24) -> str:
    """Short label for ref (TSV/HTML header): use accession or path stem; alphanumeric + underscore only."""
    s = str(spec).strip()
    p = Path(s)
    if p.exists() and (p.is_file() or p.is_dir()):
        label = p.stem if p.is_file() else p.name
    else:
        label = s
    label = re.sub(r"[^\w.]", "_", label)
    label = label.strip("._") or "ref"
    return label[:max_len]


def _looks_like_accession(s: str) -> bool:
    """True if string looks like an NCBI accession (no path, alphanumeric + dots/underscores)."""
    s = s.strip()
    if not s or "/" in s or "\\" in s:
        return False
    return bool(re.match(r"^[A-Za-z0-9_.]+$", s))


def _normalize_fasta_headers_to_accession(fasta_path: Path) -> None:
    """Rewrite FASTA so each sequence header is only the accession (first word after '>')."""
    with open(fasta_path) as f:
        content = f.read()
    records: List[Tuple[str, str]] = []
    current_header: Optional[str] = None
    current_lines: List[str] = []
    for line in content.splitlines():
        if line.startswith(">"):
            if current_header is not None:
                records.append((current_header[1:].split()[0].strip(), "\n".join(current_lines)))
            current_header = line
            current_lines = []
        else:
            current_lines.append(line)
    if current_header is not None:
        records.append((current_header[1:].split()[0].strip(), "\n".join(current_lines)))
    with open(fasta_path, "w") as f:
        for acc, seq in records:
            f.write(f">{acc}\n{seq}\n")


def _download_accessions_to_fasta(accessions: List[str], work_dir: Path) -> Path:
    """Download one or more NCBI accessions and return path to a single FASTA (concatenated if multiple)."""
    if len(accessions) == 1:
        acc = accessions[0]
        safe = re.sub(r"[^\w.-]", "_", acc)
        out_path = work_dir / f"query_{safe}.fa"
        if out_path.is_file():
            logger.info("Using cached query from accession %s: %s", acc, out_path)
            return out_path
        logger.info("Downloading query from NCBI accession: %s", acc)
        if not fetch_nucleotide_fasta(acc, out_path):
            raise FileNotFoundError(f"Failed to download accession {acc} from NCBI")
        _normalize_fasta_headers_to_accession(out_path)
        return out_path
    safe_parts = [re.sub(r"[^\w.-]", "_", a) for a in accessions]
    combined_name = "query_" + "_".join(safe_parts[:3]) + ("_etc" if len(safe_parts) > 3 else "") + ".fa"
    out_path = work_dir / combined_name
    if out_path.is_file():
        logger.info("Using cached queries from accessions %s: %s", ",".join(accessions), out_path)
        return out_path
    logger.info("Downloading %d query accessions: %s", len(accessions), ",".join(accessions))
    with open(out_path, "w") as combined:
        for acc in accessions:
            safe = re.sub(r"[^\w.-]", "_", acc)
            single_path = work_dir / f"query_single_{safe}.fa"
            if not single_path.is_file():
                if not fetch_nucleotide_fasta(acc, single_path):
                    raise FileNotFoundError(f"Failed to download accession {acc} from NCBI")
                _normalize_fasta_headers_to_accession(single_path)
            with open(single_path) as f:
                first = f.readline()
                acc_from_file = first[1:].split()[0].strip() if first.startswith(">") else acc
                combined.write(f">{acc_from_file}\n{f.read()}")
    return out_path


def resolve_query_input(input_spec: Path, work_dir: Path) -> Path:
    """
    Resolve query input: existing FASTA file, .txt file of accessions (one per line or comma-separated),
    existing directory, or NCBI accession(s) (download).
    Returns path to a FASTA file (single file or downloaded) or to a directory (caller concatenates).
    """
    if input_spec.is_file():
        if input_spec.suffix.lower() == ".txt":
            with open(input_spec) as f:
                accessions = [a.strip() for line in f for a in line.split(",") if a.strip()]
            if not accessions:
                raise FileNotFoundError(f"No accessions found in {input_spec}")
            for acc in accessions:
                if not _looks_like_accession(acc):
                    raise FileNotFoundError(f"Line in {input_spec} does not look like accession(s): {acc!r}")
            logger.info("Using %d accessions from file: %s", len(accessions), input_spec)
            return _download_accessions_to_fasta(accessions, work_dir)
        logger.info("Using query input file: %s", input_spec)
        return input_spec
    if input_spec.is_dir():
        logger.info("Using query input directory: %s", input_spec)
        return input_spec
    raw = str(input_spec).strip()
    accessions = [a.strip() for a in raw.split(",") if a.strip()]
    if not accessions:
        raise FileNotFoundError(f"Input is not an existing file or directory and has no accession(s): {input_spec}")
    for acc in accessions:
        if not _looks_like_accession(acc):
            raise FileNotFoundError(f"Input is not an existing file or directory and does not look like accession(s): {input_spec}")
    return _download_accessions_to_fasta(accessions, work_dir)


def setup_logging(output_dir: Path, verbose: bool = True) -> None:
    """Configure logging: always to hidden file .recmpox.log in output_dir; optionally to console (when verbose)."""
    log_file = output_dir / ".recmpox.log"
    handlers: List[logging.Handler] = [logging.FileHandler(log_file, mode="w")]
    if verbose:
        handlers.append(logging.StreamHandler())
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=handlers,
        force=True,
    )
    logger.info("Log file: %s", log_file)


def get_first_query_id_and_length(fasta: Path) -> Optional[Tuple[str, int]]:
    """(id, length) of first sequence in FASTA."""
    with open(fasta) as f:
        seq_id, length = None, 0
        for line in f:
            if line.startswith(">"):
                if seq_id is not None:
                    return (seq_id, length)
                seq_id = line[1:].split()[0].strip()
                length = 0
            else:
                length += len(line.strip())
        if seq_id is not None:
            return (seq_id, length)
    return None


def concatenate_fasta_dir(input_dir: Path, out_path: Path) -> Path:
    """Concatenate all .fa/.fasta/.fna in input_dir into one multi-FASTA."""
    exts = ("*.fa", "*.fasta", "*.fna")
    files = sorted(set(f for ext in exts for f in input_dir.glob(ext)))
    if not files:
        raise FileNotFoundError(f"No .fa/.fasta/.fna in {input_dir}")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as out:
        for f in files:
            with open(f) as inp:
                out.write(inp.read())
    logger.info("Concatenated %d files -> %s", len(files), out_path)
    return out_path


def _run_squirrel(clade: Optional[str], squirrel_in: Path, squirrel_out: Path, expected_aln: Path) -> None:
    """Run Squirrel on squirrel_in, output to squirrel_out; expect expected_aln to exist afterward. If clade is None, run Squirrel without --clade (mixed I/II)."""
    if not expected_aln.exists():
        # Put Squirrel temp workdir on the same filesystem as outputs to avoid
        # Snakemake mtime/clock-skew issues on some shared filesystems.
        tempdir = squirrel_out / "tmp"
        tempdir.mkdir(parents=True, exist_ok=True)
        if clade == "cladei":
            logger.info("Running Squirrel (--clade cladei) to align...")
            cmd = ["squirrel", "--clade", "cladei", str(squirrel_in), "-o", str(squirrel_out), "--tempdir", str(tempdir)]
        elif clade == "cladeii":
            logger.info("Running Squirrel (Clade II default) to align...")
            cmd = ["squirrel", str(squirrel_in), "-o", str(squirrel_out), "--tempdir", str(tempdir)]
        else:
            logger.info("Running Squirrel (no --clade; mixed or default reference) to align...")
            cmd = ["squirrel", str(squirrel_in), "-o", str(squirrel_out), "--tempdir", str(tempdir)]
        env = os.environ.copy()
        env["PYTHONNOUSERSITE"] = "1"
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600, env=env)
            if result.returncode != 0:
                # Squirrel uses Snakemake; on some shared filesystems you can get a non-zero exit due to
                # "clock skew" / coarse timestamp resolution even though the expected alignment file was
                # successfully written. If the expected output exists and stderr indicates this case,
                # treat it as a warning and continue.
                stderr_txt = (result.stderr or "").strip()
                stdout_txt = (result.stdout or "").strip()
                is_clock_skew = ("clock skew" in stderr_txt.lower()) or ("older modification time" in stderr_txt.lower())
                if expected_aln.exists() and expected_aln.stat().st_size > 0 and is_clock_skew:
                    logger.warning("Squirrel reported a filesystem mtime/clock-skew error (exit %s) but produced %s; continuing.", result.returncode, expected_aln)
                    for line in stderr_txt.splitlines():
                        logger.warning("Squirrel stderr: %s", line)
                    for line in stdout_txt.splitlines():
                        logger.info("Squirrel stdout: %s", line)
                else:
                    logger.error("Squirrel failed (exit %s).", result.returncode)
                    if stderr_txt:
                        for line in stderr_txt.splitlines():
                            logger.error("Squirrel stderr: %s", line)
                    if stdout_txt:
                        for line in stdout_txt.splitlines():
                            logger.info("Squirrel stdout: %s", line)
                    sys.exit(1)
        except FileNotFoundError as e:
            logger.error("Squirrel not found (install with: conda install -c bioconda squirrel): %s", e)
            sys.exit(1)
    if not expected_aln.exists():
        logger.error("Squirrel did not produce %s", expected_aln)
        sys.exit(1)


def _extract_tracts_as_n_full_length(
    aligned_seq: str,
    merged_tracts: List[Tuple[int, int, str, int]],
    keep_clade: str,
) -> str:
    """Return full-length (degapped) sequence with only *keep_clade* tract
    positions keeping their base; all other positions (other clade + other/ambiguous) → N.
    So only Ia or only Ib bases are retained; everything not classified as that clade becomes N.
    Positions in merged_tracts are 1-based.
    """
    seq = list(aligned_seq)
    in_keep = [False] * len(seq)
    for start_pos, end_pos, clade, _ in merged_tracts:
        if clade != keep_clade:
            continue
        for i in range(start_pos - 1, min(end_pos, len(seq))):
            in_keep[i] = True
    for i in range(len(seq)):
        if not in_keep[i] and seq[i] != "-":
            seq[i] = "N"
    return "".join(c for c in seq if c != "-")


def _extract_full_length_non_opposite(
    aligned_seq: str,
    allegiances: List[Tuple[int, str]],
    keep_clade: str,
) -> str:
    """
    Build a full-length (degapped) sequence where positions are kept for a clade
    whenever they are NOT confidently assigned to the opposite clade.

    keep_clade == "ia": keep any position that is not classified as "ib"
    keep_clade == "ib": keep any position that is not classified as "ia"
    Everything else at non-kept positions becomes N. Gaps are removed.
    """
    pos_to_all = {p: a for (p, a) in allegiances}
    seq_out: List[str] = []
    for idx, b in enumerate(aligned_seq):
        if b == "-":
            continue
        pos = idx + 1
        a = pos_to_all.get(pos, "other")
        if keep_clade == "ia":
            keep = (a != "ib")
        elif keep_clade == "ib":
            keep = (a != "ia")
        else:
            keep = True
        seq_out.append(b if keep else "N")
    return "".join(seq_out)

def _extract_tract_sequences(
    out_dir: Path,
    results: List[dict],
    alignments_queries: Dict[str, str],
    ref1_label: str,
    ref2_label: str,
    is_intra_clade: bool,
    min_consecutive: int = 1,
    line_len: int = 80,
    include_indels: bool = False,
) -> None:
    """Write two FASTA files: only Ia or only Ib tract bases kept; everything else → N.

    Extract only positions classified as that clade; positions classified as the
    other clade or as other/ambiguous are converted to N (ambiguous bases).

    When include_indels is False: only samples with recombinant_call == "potential recombinant".
    When include_indels is True: extract for all samples that have two distinct tracts.

    * ``<out_dir>/<ref1_label>_recombinant_ancestral_tract.fa`` — only ref1 (Ia) tract positions keep base; Ib + other → N.
    * ``<out_dir>/<ref2_label>_recombinant_ancestral_tract.fa`` — only ref2 (Ib) tract positions keep base; Ia + other → N.

    Tract boundaries use all positions in each sample's allegiances (including
    indel columns when -include-indels was used).  Clade labels: ``"ia"`` (ref1), ``"ib"`` (ref2).
    When is_intra_clade is True (e.g. Ia vs Ib or IIa vs IIb), we keep "all
    non-opposite" bases for each clade (Ia tract = not Ib; Ib tract = not Ia),
    so ambiguous/other positions remain in both masked sequences.

    When is_intra_clade is False (inter-clade, e.g. I vs II), we only keep
    positions confidently assigned to that clade (strict tracts) using
    _extract_tracts_as_n_full_length, so masked sequences contain only Ia or
    only Ib tract positions and everything else becomes N.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    ref1_out = out_dir / f"{ref1_label}_recombinant_ancestral_tract.fa"
    ref2_out = out_dir / f"{ref2_label}_recombinant_ancestral_tract.fa"
    n_written = 0
    n_skip_not_recombinant = 0
    n_skip_one_tract = 0
    n_skip_no_seq = 0
    with open(ref1_out, "w") as fh1, open(ref2_out, "w") as fh2:
        for r in results:
            if not include_indels and r.get("recombinant_call") != "potential recombinant":
                n_skip_not_recombinant += 1
                continue
            sample_id = r["id"]
            allegiances = r.get("allegiances", [])
            if not allegiances:
                continue
            # Derive positions from allegiances directly — includes both SNP and
            # indel column positions when -include-indels was used.
            all_positions = sorted(set(p for (p, _) in allegiances))
            merged_tracts = _compute_merged_tracts(allegiances, all_positions, min_consecutive)
            clades_present = {c for _, _, c, _ in merged_tracts}
            if len(clades_present) < 2:
                n_skip_one_tract += 1
                continue
            aligned_seq = alignments_queries.get(sample_id)
            if not aligned_seq:
                n_skip_no_seq += 1
                continue

            safe_id = _safe_fasta_id(sample_id)
            # ref1 (Ia) file
            if is_intra_clade:
                # Intra-clade: keep clade Ia bases and all positions not confidently Ib.
                seq1 = _extract_full_length_non_opposite(aligned_seq, allegiances, keep_clade="ia")
            else:
                # Inter-clade: strict Ia tracts only (Ia bases; everything else → N).
                seq1 = _extract_tracts_as_n_full_length(aligned_seq, merged_tracts, keep_clade="ia")
            len1 = len(seq1)
            non_n1 = sum(1 for b in seq1.upper() if b in "ACGT")
            cov1 = (100.0 * non_n1 / len1) if len1 else 0.0
            fh1.write(f">{safe_id}_{ref1_label}_tract_HC_{cov1:.2f}%\n")
            for i in range(0, len1, line_len):
                fh1.write(seq1[i : i + line_len] + "\n")

            # ref2 (Ib/IIb) file
            if is_intra_clade:
                # Intra-clade: keep clade Ib bases and all positions not confidently Ia.
                seq2 = _extract_full_length_non_opposite(aligned_seq, allegiances, keep_clade="ib")
            else:
                # Inter-clade: strict Ib/IIb tracts only.
                seq2 = _extract_tracts_as_n_full_length(aligned_seq, merged_tracts, keep_clade="ib")
            len2 = len(seq2)
            non_n2 = sum(1 for b in seq2.upper() if b in "ACGT")
            cov2 = (100.0 * non_n2 / len2) if len2 else 0.0
            fh2.write(f">{safe_id}_{ref2_label}_tract_HC_{cov2:.2f}%\n")
            for i in range(0, len2, line_len):
                fh2.write(seq2[i : i + line_len] + "\n")

            n_written += 1

    if n_written == 0:
        ref1_out.unlink(missing_ok=True)
        ref2_out.unlink(missing_ok=True)
        if out_dir.exists():
            try:
                out_dir.rmdir()
            except OSError:
                pass
        logger.warning(
            "--extract-tracts: no sample had two distinct tracts (skipped: %d not recombinant, %d single-tract, %d no alignment). No masked FASTA written.",
            n_skip_not_recombinant, n_skip_one_tract, n_skip_no_seq,
        )
    else:
        logger.info("--extract-tracts: wrote %d sample(s) → %s and %s", n_written, ref1_out, ref2_out)


def _load_fasta_dict(path: Path) -> Dict[str, str]:
    """Load FASTA as dict: first token of header -> full sequence (no gap stripping)."""
    return load_alignment_fasta(path)


def _run_phylogeny_pipeline(
    args: Any,
    work_dir: Path,
    ref1_label: str,
    ref2_label: str,
    squirrel_clade: Optional[str],
) -> None:
    """
    Run phylogeny: merge references + Ia/Ib partition FASTAs, align with Squirrel,
    run IQ-TREE (GTR, -bb 10000), midpoint-root, then write tree + PDF into output/phylogeny/.
    All IQ-TREE outputs and the midpoint-rooted treefile + PDF live in output/phylogeny/.
    """
    tracts_dir = args.output / "tracts"
    ref1_fa = tracts_dir / f"{ref1_label}_recombinant_ancestral_tract.fa"
    ref2_fa = tracts_dir / f"{ref2_label}_recombinant_ancestral_tract.fa"
    if not ref1_fa.exists() or not ref2_fa.exists():
        logger.warning("--phylogeny: extracted tract FASTAs not found (%s, %s); skipping phylogeny.", ref1_fa, ref2_fa)
        return

    # Use bundled reference set (no clade-based restriction)
    refs_path = PHYLOGENY_REFS_FASTA
    if not refs_path.exists():
        logger.error("--phylogeny: bundled references not found at %s", refs_path)
        sys.exit(1)

    refs_seqs = _load_fasta_dict(refs_path)
    if not refs_seqs:
        logger.error("--phylogeny: no sequences in %s", refs_path)
        sys.exit(1)

    part1_seqs = _load_fasta_dict(ref1_fa)
    part2_seqs = _load_fasta_dict(ref2_fa)
    if not part1_seqs or not part2_seqs:
        logger.warning("--phylogeny: one or both partition FASTAs empty; skipping phylogeny.")
        return

    # Phylogeny: one alignment + one tree in output/phylogeny/; all intermediates in work_dir (removed later)
    phylogeny_dir = args.output / "phylogeny"
    phylogeny_dir.mkdir(parents=True, exist_ok=True)
    combined_fa = work_dir / "phylogeny_combined.fa"
    line_len = 80
    seen_ref: Dict[str, int] = {}
    with open(combined_fa, "w") as out:
        for rid, seq in refs_seqs.items():
            safe = _safe_fasta_id(rid)
            n = seen_ref.get(safe, 0) + 1
            seen_ref[safe] = n
            header = safe if n == 1 else f"{safe}_{n}"
            out.write(f">{header}\n")
            for i in range(0, len(seq), line_len):
                out.write(seq[i : i + line_len] + "\n")
        for sid, seq in part1_seqs.items():
            # Tract FASTA headers include coverage suffix like
            # "<safe_sample_id>_<ref1_label>_tract_HC_<cov>%". Recover the base
            # sample ID so that phylogeny tract IDs can be mapped back to the
            # per-genome records by stripping the coverage part.
            base_id = sid
            cov_suffix = f"_{ref1_label}_tract_HC_"
            if cov_suffix in base_id:
                base_id = base_id.split(cov_suffix, 1)[0]
            header = _safe_fasta_id(base_id) + f"_{ref1_label}_tracts"
            out.write(f">{header}\n")
            for i in range(0, len(seq), line_len):
                out.write(seq[i : i + line_len] + "\n")
        for sid, seq in part2_seqs.items():
            base_id = sid
            cov_suffix = f"_{ref2_label}_tract_HC_"
            if cov_suffix in base_id:
                base_id = base_id.split(cov_suffix, 1)[0]
            header = _safe_fasta_id(base_id) + f"_{ref2_label}_tracts"
            out.write(f">{header}\n")
            for i in range(0, len(seq), line_len):
                out.write(seq[i : i + line_len] + "\n")

    logger.info("Wrote combined FASTA for phylogeny (refs + %s + %s partitions)", ref1_label, ref2_label)

    # Squirrel alignment in work_dir
    squirrel_out_phy = work_dir / "phylogeny_squirrel"
    squirrel_out_phy.mkdir(parents=True, exist_ok=True)
    expected_aln = squirrel_out_phy / (combined_fa.stem + ".aln.fasta")
    _run_squirrel(squirrel_clade, combined_fa, squirrel_out_phy, expected_aln)
    if not expected_aln.exists():
        logger.error("--phylogeny: Squirrel did not produce %s", expected_aln)
        sys.exit(1)
    # Single alignment in phylogeny/ (logical name)
    aln_in_phylogeny = phylogeny_dir / "phylogeny_alignment.fasta"
    shutil.copy(expected_aln, aln_in_phylogeny)

    # Compute nearest-reference outbreak labels for each tract from the alignment,
    # so we can display ancestors in the HTML when --phylogeny is used. When a
    # treefile is available and Biopython is installed, these labels will later
    # be refined using patristic distances on the tree topology.
    try:
        aln_dict = load_alignment_fasta(aln_in_phylogeny)
    except Exception as e:
        logger.warning("--phylogeny: could not load alignment for ancestor labels: %s", e)
        aln_dict = {}

    if aln_dict:
        def _outbreak_label_from_header(h: str) -> str:
            if "sh2024Ia" in h:
                return "sh2024Ia"
            if "sh2023Ib" in h:
                return "sh2023Ib"
            if "sh2017IIb" in h:
                return "sh2017IIb"
            if "_IIa_" in h or h.startswith("AY603973.1_IIa") or h.startswith("AY741551.1_IIa") or h.startswith("MN346699.1_IIa") or h.startswith("PV982292.1_IIa"):
                return "IIa"
            if "_Ia_" in h:
                return "Ia"
            if "_Ib_" in h:
                return "Ib"
            return "other"

        ref_ids = [hid for hid in aln_dict.keys() if not hid.endswith(f"_{ref1_label}_tracts") and not hid.endswith(f"_{ref2_label}_tracts")]
        tract1_ids = [hid for hid in aln_dict.keys() if hid.endswith(f"_{ref1_label}_tracts")]
        tract2_ids = [hid for hid in aln_dict.keys() if hid.endswith(f"_{ref2_label}_tracts")]

        ref_outbreak = {rid: _outbreak_label_from_header(rid) for rid in ref_ids}

        def _nearest_outbreak(seq_id: str) -> str:
            seq = aln_dict.get(seq_id)
            if not seq:
                return "other"
            best_label = "other"
            best_dist = None
            for rid in ref_ids:
                ref_seq = aln_dict.get(rid)
                if not ref_seq or len(ref_seq) != len(seq):
                    continue
                dist = 0
                n_comp = 0
                for a, b in zip(seq, ref_seq):
                    if a == "-" or b == "-":
                        continue
                    if a.upper() == "N" or b.upper() == "N":
                        continue
                    n_comp += 1
                    if a.upper() != b.upper():
                        dist += 1
                if n_comp == 0:
                    continue
                if best_dist is None or dist < best_dist:
                    best_dist = dist
                    best_label = ref_outbreak.get(rid, "other")
            return best_label

        if not hasattr(args, "_phylogeny_ancestors"):
            args._phylogeny_ancestors = {}

        def _sample_key_from_tract(hid: str, label: str) -> str:
            suffix = f"_{label}_tracts"
            return hid[: -len(suffix)] if hid.endswith(suffix) else hid

        for hid in tract1_ids:
            key = _sample_key_from_tract(hid, ref1_label)
            anc = _nearest_outbreak(hid)
            rec = args._phylogeny_ancestors.get(key, {})
            # Keep both alignment-based and (optionally) tree-based labels; ref1
            # / ref2 will be overwritten by tree-based labels later when a tree
            # is available so HTML shows the topology-aware nearest ancestor.
            rec["ref1_aln"] = anc
            # Default display = alignment; may be overridden by tree-based label.
            rec.setdefault("ref1", anc)
            args._phylogeny_ancestors[key] = rec

        for hid in tract2_ids:
            key = _sample_key_from_tract(hid, ref2_label)
            anc = _nearest_outbreak(hid)
            rec = args._phylogeny_ancestors.get(key, {})
            rec["ref2_aln"] = anc
            rec.setdefault("ref2", anc)
            args._phylogeny_ancestors[key] = rec

    # IQ-TREE in work_dir (so only one treefile ends up in phylogeny/)
    iqtree_dir = work_dir / "phylogeny_iqtree"
    iqtree_dir.mkdir(parents=True, exist_ok=True)
    iqtree_prefix = iqtree_dir / "alignment"
    try:
        # -pre: output file prefix (alignment.treefile, alignment.iqtree, alignment.log, etc.)
        # -czb: collapse zero-length branches into polytomies
        cmd_iqtree = [
            "iqtree",
            "-s", str(aln_in_phylogeny),
            "-bb", "10000",
            "-pre", str(iqtree_prefix),
            "-czb",
        ]
        # Use all available CPUs by default; if user requests >1 threads, honor it.
        if getattr(args, "threads", None) is not None and int(args.threads) > 1:
            cmd_iqtree.extend(["-nt", str(args.threads)])
        else:
            cmd_iqtree.extend(["-nt", "AUTO"])
        logger.info("Running IQ-TREE: %s", " ".join(cmd_iqtree))
        result = subprocess.run(cmd_iqtree, timeout=7200)
        if result.returncode != 0:
            logger.error("IQ-TREE failed (exit %s).", result.returncode)
            sys.exit(1)
    except FileNotFoundError as e:
        logger.error("iqtree not found (install e.g. conda install -c bioconda iqtree2): %s", e)
        sys.exit(1)

    treefile = iqtree_prefix.with_suffix(".treefile")
    if not treefile.exists():
        treefile = Path(str(iqtree_prefix) + ".treefile")
    if not treefile.exists():
        logger.error("--phylogeny: IQ-TREE did not produce treefile at %s", treefile)
        sys.exit(1)

    # Midpoint root and PDF via R only (root_tree_figure.R); IQ-TREE outputs alignment.treefile
    out_tree_path = phylogeny_dir / "phylogeny_tree.treefile"
    pdf_path = phylogeny_dir / "phylogeny_tree.pdf"
    svg_path = phylogeny_dir / "phylogeny_tree.svg"
    r_figure_pdf = phylogeny_dir / "tree_figure.pdf"
    r_figure_svg = phylogeny_dir / "tree_figure.svg"
    r_rooted = phylogeny_dir / "rooted.tree"
    pdf_written = False
    midpoint_rooted = False

    # R expects input tree in phylogeny dir; use temp name then keep only final outputs
    tree_in_phylogeny = phylogeny_dir / "alignment.treefile"
    if not treefile.exists():
        logger.error("--phylogeny: IQ-TREE treefile not found at %s", treefile)
        sys.exit(1)
    shutil.copy(treefile, tree_in_phylogeny)

    if ROOT_TREE_FIGURE_R.exists():
        try:
            logger.info("Running R script for midpoint rooting and tree figure (root_tree_figure.R)")
            r_result = subprocess.run(
                ["Rscript", str(ROOT_TREE_FIGURE_R), str(tree_in_phylogeny)],
                cwd=str(phylogeny_dir),
                timeout=300,
                capture_output=True,
                text=True,
            )
            if r_result.returncode == 0 and r_rooted.exists():
                shutil.copy(r_rooted, out_tree_path)
                midpoint_rooted = True
                if r_figure_pdf.exists():
                    shutil.copy(r_figure_pdf, pdf_path)
                    pdf_written = True
                if r_figure_svg.exists():
                    shutil.copy(r_figure_svg, svg_path)
                logger.info("Midpoint-rooted tree and PDF from R: %s", out_tree_path)
            elif r_result.returncode != 0 and r_result.stderr:
                logger.warning("R tree figure script failed: %s", r_result.stderr[:200])
                shutil.copy(treefile, out_tree_path)
                logger.info("Wrote tree (unrooted): %s", out_tree_path)
        except FileNotFoundError:
            logger.warning("Rscript not found; install R and run root_tree_figure.R for midpoint root and PDF")
            shutil.copy(treefile, out_tree_path)
            logger.info("Wrote tree (unrooted): %s", out_tree_path)
        except subprocess.TimeoutExpired:
            logger.warning("R tree figure script timed out.")
            shutil.copy(treefile, out_tree_path)
            logger.info("Wrote tree (unrooted): %s", out_tree_path)
        except Exception as r_e:
            logger.warning("R tree figure failed: %s", r_e)
            shutil.copy(treefile, out_tree_path)
            logger.info("Wrote tree (unrooted): %s", out_tree_path)
    else:
        logger.warning("root_tree_figure.R not found at %s; writing unrooted tree", ROOT_TREE_FIGURE_R)
        shutil.copy(treefile, out_tree_path)
        logger.info("Wrote tree (unrooted): %s", out_tree_path)

    # Keep only: phylogeny_alignment.fasta, phylogeny_tree.treefile, phylogeny_tree.pdf, phylogeny_tree.svg
    for f in (tree_in_phylogeny, r_rooted, r_figure_pdf, r_figure_svg):
        if f.exists():
            try:
                f.unlink()
            except OSError:
                pass

    # Optional: refine nearest-ancestor labels using patristic distances on the tree
    if aln_dict and _HAS_BIOPHYLO and out_tree_path.exists():
        dists = _load_tree_and_distances(out_tree_path)
        if dists:
            logger.info("--phylogeny: refining ancestors using tree topology (patristic distances)")

            def _nearest_outbreak_tree(seq_id: str) -> Optional[str]:
                best_label: Optional[str] = None
                best_dist: Optional[float] = None
                for rid in ref_ids:
                    d = _patristic_dist(seq_id, rid, dists)
                    if d is None:
                        continue
                    if best_dist is None or d < best_dist:
                        best_dist = d
                        best_label = ref_outbreak.get(rid, "other")
                return best_label

            def _sample_key_from_tract(hid: str, label: str) -> str:
                suffix = f"_{label}_tracts"
                return hid[: -len(suffix)] if hid.endswith(suffix) else hid

            for hid in tract1_ids:
                key = _sample_key_from_tract(hid, ref1_label)
                anc_tree = _nearest_outbreak_tree(hid)
                if not anc_tree:
                    continue
                rec = args._phylogeny_ancestors.get(key, {})
                # Preserve alignment-based label, but show tree-based label in HTML.
                rec["ref1_tree"] = anc_tree
                rec["ref1"] = anc_tree
                args._phylogeny_ancestors[key] = rec

            for hid in tract2_ids:
                key = _sample_key_from_tract(hid, ref2_label)
                anc_tree = _nearest_outbreak_tree(hid)
                if not anc_tree:
                    continue
                rec = args._phylogeny_ancestors.get(key, {})
                rec["ref2_tree"] = anc_tree
                rec["ref2"] = anc_tree
                args._phylogeny_ancestors[key] = rec

    root_label = "midpoint-rooted tree" if midpoint_rooted else "tree (unrooted)"
    if pdf_written:
        print(f"  Phylogeny: {root_label} {out_tree_path}; PDF {pdf_path}")
    else:
        print(f"  Phylogeny: {root_label} {out_tree_path}; PDF not created (would be {pdf_path}; open .treefile in FigTree to export PDF)")


def _write_all_sequences_fasta(
    out_path: Path,
    ref_ia_key: str,
    ref_ib_key: str,
    ref_ia_seq: str,
    ref_ib_seq: str,
    alignments_queries: Dict[str, str],
    ref_keys_to_skip: set,
    line_len: int = 80,
) -> None:
    """Write one FASTA with ref1 + ref2 + all query sequences (aligned)."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        for sid, seq in [(ref_ia_key, ref_ia_seq), (ref_ib_key, ref_ib_seq)]:
            f.write(f">{sid}\n")
            for i in range(0, len(seq), line_len):
                f.write(seq[i : i + line_len] + "\n")
        for qid, qseq in alignments_queries.items():
            if qid in ref_keys_to_skip:
                continue
            f.write(f">{qid}\n")
            for i in range(0, len(qseq), line_len):
                f.write(qseq[i : i + line_len] + "\n")


def _write_ref_alignment(
    out_path: Path,
    seq_ia: str,
    seq_ib: str,
    id_ia: str,
    id_ib: str,
    line_len: int = 80,
) -> None:
    """Write ref1 (Ia) and ref2 (Ib) aligned sequences to FASTA for visual verification."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        f.write(f">ref1_Ia {id_ia}\n")
        for i in range(0, len(seq_ia), line_len):
            f.write(seq_ia[i : i + line_len] + "\n")
        f.write(f">ref2_Ib {id_ib}\n")
        for i in range(0, len(seq_ib), line_len):
            f.write(seq_ib[i : i + line_len] + "\n")


def _write_indel_regions_side_by_side(
    out_path: Path,
    seq_ia: str,
    seq_ib: str,
    indels: List[Tuple[int, int, str]],
    line_len: int = 80,
) -> None:
    """Write each diagnostic indel region with ref1 and ref2 side-by-side for visual verification."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        f.write("# Diagnostic indel regions: ref1 (Ia) vs ref2 (Ib) in alignment coordinates.\n")
        f.write("# ref_with_bases = ref that has sequence (ACGT); deletion_in = ref that has gaps (-).\n\n")
        for start, end, ref_who in indels:
            start0 = start - 1
            end0 = end
            seg_ia = seq_ia[start0:end0]
            seg_ib = seq_ib[start0:end0]
            deletion_in = "ib" if ref_who == "ia" else "ia"
            f.write(f"## Region {start}-{end} (ref_with_bases={ref_who}, deletion_in={deletion_in})\n")
            for i in range(0, len(seg_ia), line_len):
                s_ia = seg_ia[i : i + line_len]
                s_ib = seg_ib[i : i + line_len]
                pos = start + i
                f.write(f"  ref1_Ia {pos:>6}: {s_ia}\n")
                f.write(f"  ref2_Ib {pos:>6}: {s_ib}\n")
            f.write("\n")


def _snp_positions_svg(positions: List[int], genome_length: int, width_units: int = 1000, height: int = 44) -> str:
    """Build an SVG showing diagnostic SNP positions along the genome; each tick has a title for hover.
    kbp scale marks and SNP ticks are drawn below the baseline, labels at the bottom.
    """
    if genome_length <= 0 or not positions:
        return ""
    y_line = 6           # horizontal baseline near the top
    y_kbp_tick = 16      # bottom of kbp tick (10 px below baseline)
    y_kbp_label = 30     # kbp label baseline (below the tick)
    y_snp_bottom = 19    # bottom of SNP tick (13 px below baseline)

    # Same adaptive step logic as _genome_ruler_html
    step_bp = 10000
    for s in [1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000]:
        if genome_length / s <= 20:
            step_bp = s
            break

    hit_width = max(4, width_units // 80)
    parts = [
        f'<svg class="snp-positions-svg" viewBox="0 0 {width_units} {height}" preserveAspectRatio="xMidYMid meet" style="max-width:100%; height:auto;">',
        f'<line x1="0" y1="{y_line}" x2="{width_units}" y2="{y_line}" stroke="#333" stroke-width="1.5"/>',
    ]

    # kbp scale marks (downward from baseline)
    kbp_pos = 0
    while kbp_pos <= genome_length:
        # Skip regular ticks too close to the endpoint (would produce duplicate label)
        if kbp_pos > 0 and (genome_length - kbp_pos) < step_bp * 0.6:
            kbp_pos += step_bp
            continue
        x = (kbp_pos / genome_length) * width_units
        label = "0" if kbp_pos == 0 else f"{kbp_pos // 1000}k"
        parts.append(
            f'<line x1="{x:.1f}" y1="{y_line}" x2="{x:.1f}" y2="{y_kbp_tick}"'
            f' stroke="#adb5bd" stroke-width="0.8" pointer-events="none"/>'
        )
        anchor = "start" if kbp_pos == 0 else "middle"
        parts.append(
            f'<text x="{x:.1f}" y="{y_kbp_label}" font-size="9" fill="#888"'
            f' text-anchor="{anchor}">{label}</text>'
        )
        kbp_pos += step_bp
    # Final tick at exact genome end if not already aligned
    if genome_length % step_bp != 0:
        x_end = width_units
        parts.append(
            f'<line x1="{x_end}" y1="{y_line}" x2="{x_end}" y2="{y_kbp_tick}"'
            f' stroke="#adb5bd" stroke-width="0.8" pointer-events="none"/>'
        )
        parts.append(
            f'<text x="{x_end}" y="{y_kbp_label}" font-size="9" fill="#888"'
            f' text-anchor="end">{genome_length:,}</text>'
        )

    # SNP ticks (downward from baseline)
    for pos in positions:
        x = (pos / genome_length) * width_units
        x = max(0, min(width_units, x))
        rx = max(0, x - hit_width / 2)
        rw = min(hit_width, width_units - rx)
        parts.append(
            f'<g><title>Position: {pos} bp</title>'
            f'<rect x="{rx:.2f}" y="0" width="{rw:.2f}" height="{height}" fill="transparent" class="snp-tick-hit"/>'
            f'<line x1="{x:.2f}" y1="{y_line}" x2="{x:.2f}" y2="{y_snp_bottom}" stroke="#667eea" stroke-width="1" pointer-events="none"/>'
            f'</g>'
        )

    parts.append("</svg>")
    return "\n".join(parts)


def _snp_positions_histogram_bins(
    positions: List[int], genome_length: int, num_bins: int = 60
) -> Tuple[List[str], List[int]]:
    """Bin diagnostic SNP positions along the genome for a bar chart. Returns (labels, counts)."""
    if genome_length <= 0 or not positions or num_bins < 1:
        return [], []
    bin_width = genome_length / num_bins
    counts = [0] * num_bins
    for pos in positions:
        idx = min(int((pos - 1) / bin_width), num_bins - 1) if pos >= 1 else 0
        counts[idx] += 1
    labels = []
    for i in range(num_bins):
        start_bp = int(i * bin_width)
        end_bp = int((i + 1) * bin_width)
        if end_bp >= genome_length:
            end_bp = genome_length
        start_k = start_bp / 1000
        end_k = end_bp / 1000
        labels.append(f"{start_k:.0f}k–{end_k:.0f}k")
    return labels, counts


def _compute_merged_tracts(
    allegiances: List[Tuple[int, str]],
    diagnostic_snp_positions: List[int],
    min_consecutive: int = 1,
) -> List[Tuple[int, int, str, int]]:
    """Return merged, sustained recombination tracts from per-sample allegiances.

    Returns a list of (start_pos, end_pos, clade, n_snps) in alignment-column coordinates
    (1-based). Only tracts with >= min_consecutive SNPs are kept.

    IMPORTANT: Tracts are built from *consecutive* diagnostic positions with the same
    clade call (ia/ib) in the full series of diagnostic sites. Any intervening
    'other'/ambiguous site breaks the tract – we no longer connect Ia…Ia or Ib…I
    across an 'other' position.
    """
    if not allegiances or not diagnostic_snp_positions:
        return []

    # Map position -> allegiance for quick lookup
    pos_to_all = {p: a for (p, a) in allegiances}

    def _norm(a: str) -> str:
        if a == "ia":
            return "ia"
        if a == "ib":
            return "ib"
        # Treat ambiguous / other_n / anything else as "other"
        return "other"

    tracts: List[Tuple[int, int, str, int]] = []
    cur_clade: Optional[str] = None
    cur_start: Optional[int] = None
    cur_end: Optional[int] = None
    cur_count: int = 0

    for pos in sorted(diagnostic_snp_positions):
        a_raw = pos_to_all.get(pos, "other")
        a = _norm(a_raw)
        if a in ("ia", "ib"):
            if cur_clade == a:
                # Extend current tract
                cur_end = pos
                cur_count += 1
            else:
                # Finish previous tract (if any)
                if cur_clade in ("ia", "ib") and cur_start is not None and cur_end is not None:
                    tracts.append((cur_start, cur_end, cur_clade, cur_count))
                # Start new tract
                cur_clade = a
                cur_start = cur_end = pos
                cur_count = 1
        else:
            # 'other' / ambiguous / N → break any current tract
            if cur_clade in ("ia", "ib") and cur_start is not None and cur_end is not None:
                tracts.append((cur_start, cur_end, cur_clade, cur_count))
            cur_clade = None
            cur_start = None
            cur_end = None
            cur_count = 0

    # Finalize last tract
    if cur_clade in ("ia", "ib") and cur_start is not None and cur_end is not None and cur_count > 0:
        tracts.append((cur_start, cur_end, cur_clade, cur_count))

    # Keep only tracts with enough consecutive SNPs
    return [t for t in tracts if t[3] >= min_consecutive]


def _genome_ruler_html(genome_length: int, min_width: int) -> str:
    """Return an HTML ruler div with kbp tick marks proportional to genome_length."""
    if not genome_length:
        return ""
    # Pick a step size that gives 10–20 ticks
    step_bp = 10000
    for s in [1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000]:
        if genome_length / s <= 20:
            step_bp = s
            break
    ticks = []
    pos = 0
    while pos <= genome_length:
        pct = pos / genome_length * 100
        label = "0" if pos == 0 else f"{pos // 1000}k"
        # Skip regular ticks that would crowd the endpoint label (within 60 % of one step)
        if pos > 0 and (genome_length - pos) < step_bp * 0.6:
            pos += step_bp
            continue
        ticks.append(f'<span class="ruler-tick" style="left:{pct:.2f}%">{label}</span>')
        pos += step_bp
    # Always include a tick at the genome end if not already there, showing exact bp
    if genome_length % step_bp != 0:
        ticks.append(f'<span class="ruler-tick" style="left:100%">{genome_length:,}</span>')
    return f'<div class="strip-ruler" style="min-width:{min_width}px">{"".join(ticks)}</div>'


def _write_results_html(
    out_path: Path,
    results: List[Dict[str, Any]],
    ref1_label: str,
    ref2_label: str,
    recombinant_threshold_note: Optional[str] = None,
    other_explanation: Optional[str] = None,
    is_intra_clade: bool = True,
    minor_threshold: float = 10.0,
    breakpoint_min_consecutive_snps: int = 1,
    part_index: Optional[int] = None,
    total_parts: Optional[int] = None,
    n_diagnostic_snps: Optional[int] = None,
    n_indel_columns: Optional[int] = None,
    ref1_spec: Optional[str] = None,
    ref2_spec: Optional[str] = None,
    diagnostic_snp_positions: Optional[List[int]] = None,
    genome_length: Optional[int] = None,
    phylogeny_pdf_path: Optional[Path] = None,
    phylogeny_pdf_base64: Optional[str] = None,
    phylogeny_svg_content: Optional[str] = None,
) -> None:
    """Write Virasign-style HTML: container, gradient header, sortable table, Chart.js stacked bar, recombinant column."""
    import json
    from datetime import datetime
    out_path.parent.mkdir(parents=True, exist_ok=True)
    gen_time = datetime.utcnow().strftime("%Y-%m-%d")
    n_diagnostic_sites = results[0]["n_diagnostic_snps"] if results else 0
    part_html = ""
    if part_index is not None and total_parts is not None and total_parts > 1:
        part_html = f'<div class="stat-box"><div class="number">Part {part_index} of {total_parts}</div><div class="label">This file</div></div>'
    # Summary breakdown: SNPs, diagnostic indels (columns), total (only when indels included)
    sites_breakdown_html = ""
    if n_diagnostic_snps is not None and n_indel_columns is not None:
        sites_breakdown_html = (
            '<div class="stat-box"><div class="number">{n_snps}</div><div class="label">Diagnostic SNPs</div></div>'
            '<div class="stat-box"><div class="number">{n_indel}</div><div class="label">Diagnostic indels</div></div>'
            '<div class="stat-box"><div class="number">{n_total}</div><div class="label">Diagnostic sites (total)</div></div>'
        ).format(n_snps=n_diagnostic_snps, n_indel=n_indel_columns, n_total=n_diagnostic_sites)
    else:
        sites_breakdown_html = (
            '<div class="stat-box"><div class="number">{n_total}</div><div class="label">Diagnostic sites (SNPs only)</div></div>'
        ).format(n_total=n_diagnostic_sites)
    snps_only_note_html = ""
    # When a phylogeny was run and per-sample ancestor labels are available,
    # results entries may contain an \"ancestors\" field with \"ref1\" and \"ref2\" labels.
    has_ancestors = any("ancestors" in r for r in results)
    cols = [
        ("id", "Sample ID", False),
        ("length", "Length (bp)", True),
        ("n_diagnostic_snps", "Diagnostic sites", True),
        (("n_ia", "pct_ia"), f"{ref1_label} (n | %)", True),
        (("n_ib", "pct_ib"), f"{ref2_label} (n | %)", True),
        (("n_other", "pct_other"), "other (n | %)", True),
    ]
    if has_ancestors:
        cols.append(("ancestors", "ancestors", False))
    else:
        cols.append(("consensus_snp", "consensus (SNP)", False))
    cols.append(
        ("recombinant_call", "recombinant", False),
    )
    html_escape = (lambda s: str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;"))
    refs_summary_html = ""
    if ref1_spec or ref2_spec:
        def _ref_display_name(spec: Optional[str]) -> str:
            """Show filename only for file paths, else spec as-is (e.g. accession)."""
            if not spec:
                return spec or ""
            return Path(spec).name if "/" in str(spec) else spec

        def _ref_box(label: str, spec: Optional[str]) -> str:
            if not spec:
                return f'<div class="stat-box"><div class="number">{html_escape(label)}</div><div class="label">Ref</div></div>'
            display_spec = _ref_display_name(spec)
            spec_esc = html_escape(display_spec)
            if _looks_like_accession(spec):
                url = "https://www.ncbi.nlm.nih.gov/nuccore/" + urllib.parse.quote(spec, safe="")
                content = f'<a class="accession-link" href="{url}" target="_blank" rel="noopener">{spec_esc}</a>'
            else:
                content = spec_esc
            default_name = REF_DEFAULT_NAMES.get(spec, "")
            label_str = html_escape(label) + (f", {html_escape(default_name)}" if default_name else "")
            return f'<div class="stat-box"><div class="number">{content}</div><div class="label">Ref ({label_str})</div></div>'
        refs_summary_html = _ref_box(ref1_label, ref1_spec) + _ref_box(ref2_label, ref2_spec)
    threshold_html = ""
    if recombinant_threshold_note:
        threshold_html = f'<p class="threshold-note">{html_escape(recombinant_threshold_note)}</p>'
    if other_explanation:
        threshold_html += f'<p class="threshold-note">{html_escape(other_explanation)}</p>'
    if n_indel_columns is not None:
        threshold_html += (
            '<p class="threshold-note">When including diagnostic indels, good coverage of those regions is recommended for reliable classification (N and gap at indel columns are assigned to the ref that has the deletion).</p>'
        )

    rows_html = []
    for ri, r in enumerate(results):
        cells = []
        for i, (key, _label, is_num) in enumerate(cols):
            cls = ' class="num"' if is_num else ""
            if isinstance(key, tuple):
                key1, key2 = key
                val1 = r.get(key1, "")
                val2 = r.get(key2, "")
                cells.append(f"<td{cls}>{html_escape(str(val1))} | {html_escape(str(val2))}</td>")
            elif key == "id":
                val_str = str(r.get(key, "")).strip()
                display_str = val_str[:20] + ("\u2026" if len(val_str) > 20 else "")
                title_attr = f' title="{html_escape(val_str)}"' if len(val_str) > 20 else ""
                if val_str and _looks_like_accession(val_str):
                    url = "https://www.ncbi.nlm.nih.gov/nuccore/" + urllib.parse.quote(val_str, safe="")
                    cells.append(f'<td{cls}><a class="accession-link" href="{html_escape(url)}" target="_blank" rel="noopener"{title_attr}>{html_escape(display_str)}</a></td>')
                else:
                    cells.append(f"<td{cls}{title_attr}>{html_escape(display_str)}</td>")
            elif key == "recombinant_call":
                rec = str(r.get(key, ""))
                if rec == "recombinant":
                    badge_cls = "recombinant-badge recombinant"
                elif rec == "potential recombinant":
                    badge_cls = "recombinant-badge potential"
                else:
                    badge_cls = "recombinant-badge no"
                cells.append(f'<td{cls}><span class="{badge_cls}">{html_escape(rec)}</span></td>')
            else:
                if key == "ancestors":
                    val = r.get(key, "")
                    # When phylogeny ancestors are available, we store a small dict
                    # with "ref1" and "ref2" labels (and optionally *_aln / *_tree).
                    if isinstance(val, dict):
                        v1 = val.get("ref1", "")
                        v2 = val.get("ref2", "")
                        cells.append(f"<td{cls}>{html_escape(str(v1))} | {html_escape(str(v2))}</td>")
                    else:
                        cells.append(f"<td{cls}>{html_escape(str(val))}</td>")
                else:
                    cells.append(f"<td{cls}>{html_escape(str(r.get(key, '')))}</td>")
        rec = r.get("recombinant_call", "")
        row_cls = " class=\"recombinant\"" if rec == "potential recombinant" else ""
        rows_html.append(f'<tr data-row="{ri}" data-recombinant="{html_escape(rec)}"{row_cls}>' + "".join(cells) + "</tr>")

    th_cells = []
    for i, (_, label, is_num) in enumerate(cols):
        cls = ' class="num sortable"' if is_num else ' class="sortable"'
        th_cells.append(f"<th{cls} data-col=\"{i}\">{html_escape(label)}</th>")
    thead = "<tr>" + "".join(th_cells) + "</tr>"
    tbody = "\n".join(rows_html)

    filter_cells = []
    for i, (key, _label, is_num) in enumerate(cols):
        if i == 0:
            filter_cells.append('<td></td>')
        elif key == "recombinant_call":
            # Three options only: All, recombinant (or potential recombinant when no phylogeny), no recombinant.
            if has_ancestors:
                rec_options = (
                    '<option value="all">All</option>'
                    '<option value="recombinant">Recombinants only</option>'
                    '<option value="no recombinant">No recombinants only</option>'
                )
            else:
                rec_options = (
                    '<option value="all">All</option>'
                    '<option value="potential recombinant">Potential recombinants only</option>'
                    '<option value="no recombinant">No recombinants only</option>'
                )
            filter_cells.append(f'<td><select id="recFilter" class="rec-filter-select">{rec_options}</select></td>')
        elif not is_num:
            filter_cells.append('<td></td>')
        else:
            filter_cells.append(f'<td class="num"><input type="number" step="any" placeholder="min" data-col="{i}" class="min-thresh"/></td>')
    filter_row = "<tr id=\"filterrow\">" + "".join(filter_cells) + "</tr>"

    # Build list of dicts per row for JS (exclude allegiances to keep JSON small)
    # Flatten tuple keys (merged display columns) and always include chart keys
    _chart_keys = {"id", "pct_ia", "pct_ib", "pct_other", "n_ia", "n_ib", "n_other", "n_other_n"}
    _data_keys = {k for c in cols for k in (c[0] if isinstance(c[0], tuple) else (c[0],))} | _chart_keys
    data_list = [{k: r.get(k, "") for k in _data_keys} for r in results]
    data_json = json.dumps(data_list).replace("</", "<\\/")

    # Diagnostic sites per sample: strip (genome position, color Ia/Ib/other) + table for ALL consensus genomes
    genome_length = results[0]["length"] if results else 0
    # Trim strip display to the last diagnostic SNP position (any allegiance, incl. other)
    if diagnostic_snp_positions and genome_length:
        display_length = max(diagnostic_snp_positions)
    else:
        display_length = genome_length
    _first_alle = next((r.get("allegiances", []) for r in results if r.get("allegiances")), [])
    strip_min_w = max(600, len(_first_alle) * 2)
    rec_sites_html = ""
    for ri, r in enumerate(results):
        sample_id = r.get("id", "")
        allegiances = r.get("allegiances", [])
        rec_call = r.get("recombinant_call", "")
        if not allegiances:
            continue
        sorted_alle = sorted(allegiances, key=lambda x: x[0])
        strip_segments = ""
        for (pos, allegiance) in sorted_alle:
            # Use "other-n" for N's so the strip filter can show only N segments when "Other (N's)" is selected
            if allegiance == "ia":
                cls = "ia"
            elif allegiance == "ib":
                cls = "ib"
            elif allegiance == "other_n":
                cls = "other other-n"
            else:
                cls = "other"
            lbl = ref1_label if allegiance == "ia" else (ref2_label if allegiance == "ib" else "other")
            pct = pos / display_length * 100 if display_length else 0
            strip_segments += f'<span class="strip-segment {cls}" title="{pos} bp – {html_escape(lbl)}" style="left:{pct:.3f}%"></span>'
        section_cls = "rec-sites-section" + (" recombinant" if rec_call == "potential recombinant" else "")
        ruler_html = _genome_ruler_html(display_length, strip_min_w)
        rec_sites_html += (
            f'<div class="{section_cls}" data-row="{ri}" data-recombinant="{html_escape(rec_call)}">'
            f'<div class="rec-sites-row">'
            f'<span class="rec-sites-sample-id" title="{html_escape(sample_id)}">{html_escape(sample_id)}</span>'
            f'<div class="strip-cell"><div class="strip-genome" style="min-width:{strip_min_w}px" role="img" aria-label="Diagnostic sites along genome">{strip_segments}</div>{ruler_html}</div>'
            f'</div>'
            f'<details class="rec-sites-details"><summary>Show diagnostic site table (by tract)</summary>'
            f'<p class="threshold-note">Tracts = consecutive diagnostic sites with same classification. One row per tract.</p>'
            f'<table class="rec-sites-table"><thead><tr><th>Start (bp)</th><th>End (bp)</th><th>Clade</th><th>Sites</th></tr></thead><tbody>'
        )
        # Group consecutive positions with same allegiance into tracts (treat other_n and ambiguous as "other")
        def _eff_allegiance(a: str) -> str:
            return "other" if a in ("ambiguous", "other_n") else a

        i = 0
        while i < len(sorted_alle):
            start_pos, allegiance = sorted_alle[i]
            end_pos = start_pos
            eff = _eff_allegiance(allegiance)
            j = i + 1
            while j < len(sorted_alle) and _eff_allegiance(sorted_alle[j][1]) == eff:
                end_pos = sorted_alle[j][0]
                j += 1
            label = ref1_label if eff == "ia" else (ref2_label if eff == "ib" else "other")
            n_sites = j - i
            rec_sites_html += f'<tr><td class="num">{start_pos}</td><td class="num">{end_pos}</td><td>{html_escape(label)}</td><td class="num">{n_sites}</td></tr>'
            i = j
        rec_sites_html += "</tbody></table></details></div>"
    if rec_sites_html:
        sections_html = rec_sites_html
        # Use reference names (labels or accessions) in strip legend so it's clear which ref is blue/orange
        strip_ref1_name = ref1_label if ref1_label not in ("ref1", "ref2") else (ref1_spec or ref1_label)
        strip_ref2_name = ref2_label if ref2_label not in ("ref1", "ref2") else (ref2_spec or ref2_label)
        rec_sites_html = (
            '<details class="collapsible-section diagnostic-strips-chart" open id="diagnosticStripsSection">'
            f'<summary><h2>Classification of diagnostic sites per sample</h2><button class="pdf-btn" onclick="event.stopPropagation();exportStripSvg(\'diagnosticStripsSection\',\'diagnosticStripsContainer\',\'Classification of diagnostic sites per sample\',\'{html_escape(ref1_label)}\',\'{html_escape(ref2_label)}\')">&#8595; Download SVG</button></summary>'
            '<div class="section-inner chart-section">'
            '<div class="strip-classification-filter-wrap">'
            '<label class="strip-filter-check-label"><input type="checkbox" id="stripClassificationFilterCheck" aria-describedby="stripClassificationOptions"> Filter strips by classification</label>'
            '<div id="stripClassificationOptions" class="strip-classification-options" aria-hidden="true">'
            '<label><input type="checkbox" class="strip-class-opt" data-filter="ia"> {ref1}</label>'
            '<label><input type="checkbox" class="strip-class-opt" data-filter="ib"> {ref2}</label>'
            '<label><input type="checkbox" class="strip-class-opt" data-filter="other_all"> Other (all)</label>'
            '<label><input type="checkbox" class="strip-class-opt" data-filter="other_n"> Other (N\'s)</label>'
            '</div>'
            '</div>'
            '<p class="threshold-note">One strip per consensus: each segment = one diagnostic site in genomic order. <span id="stripFilterCount" aria-live="polite"></span></p>'
            '<p class="threshold-note">{ref1} (blue), {ref2} (orange), other (gray).</p>'
            '<div class="strip-legend"><span class="strip-legend-ia"></span> {ref1} &nbsp; <span class="strip-legend-ib"></span> {ref2} &nbsp; <span class="strip-legend-other"></span> other</div>'
            '<div class="strip-strips-container" id="stripScrollWrapper">'
            '<div id="diagnosticStripsContainer">'
        ).format(ref1=html_escape(strip_ref1_name), ref2=html_escape(strip_ref2_name))
        rec_sites_html = rec_sites_html + sections_html + "</div></div></div></details>"

    # Breakpoints per sample: regions (runs) with breakpoints marked (optional consecutive-SNP filtering)
    breakpoints_section_html = ""
    if diagnostic_snp_positions and results:
        min_consecutive = max(1, int(breakpoint_min_consecutive_snps))
        no_regions_placeholder = '<span class="threshold-note">No regions</span>'
        breakpoints_sections_html = ""
        for ri, r in enumerate(results):
            allegiances = r.get("allegiances", [])
            sample_id = r.get("id", "")
            rec_call = r.get("recombinant_call", "")
            merged_tracts = _compute_merged_tracts(allegiances, diagnostic_snp_positions, min_consecutive)
            bp_strip_min_w = max(600, display_length // 150) if display_length else 600
            strip_segments = ""
            for j, (start_pos, end_pos, clade, n_snps) in enumerate(merged_tracts):
                cls = "ia" if clade == "ia" else "ib"
                lbl = ref1_label if clade == "ia" else ref2_label
                left_pct = start_pos / display_length * 100 if display_length else 0
                width_pct = max(0.3, (end_pos - start_pos + 1) / display_length * 100) if display_length else 2
                strip_segments += (
                    f'<span class="strip-segment region-segment {cls}" title="{start_pos}–{end_pos} {html_escape(lbl)} ({n_snps} SNPs)" style="left:{left_pct:.3f}%; width:{width_pct:.3f}%;"></span>'
                )
            n_tracts = len(merged_tracts)
            n_breakpoints = max(0, n_tracts - 1)
            section_cls = "rec-sites-section" + (" recombinant" if rec_call == "potential recombinant" else "")
            # Single tract = genome entirely one clade → show "nothing present" (no recombination)
            if n_tracts == 0:
                strip_display = no_regions_placeholder
                summary_text = "Show recombination tracts (Number of tracts: 0, breakpoints: 0)"
                details_content = f'<p class="threshold-note">No sustained tracts (≥{min_consecutive} consecutive SNPs) in this genome.</p>'
            elif n_tracts == 1:
                strip_display = '<span class="threshold-note">No recombination (genome entirely one clade)</span>'
                summary_text = "No recombination tracts (genome entirely one clade)"
                details_content = '<p class="threshold-note">No recombination detected; genome is entirely one clade.</p>'
            else:
                bp_ruler_html = _genome_ruler_html(display_length, bp_strip_min_w)
                strip_display = (
                    f'<div class="strip-genome breakpoints-strip" style="min-width:{bp_strip_min_w}px" role="img" aria-label="Predicted regions and breakpoints">{strip_segments}</div>'
                    + bp_ruler_html
                )
                summary_text = f"Show recombination tracts (Number of tracts: {n_tracts}, breakpoints: {n_breakpoints})"
                details_content = (
                    f'<table class="rec-sites-table"><thead><tr><th>Tract #</th><th>Beginning of tract (bp)</th><th>End of tract (bp)</th><th>Clade</th></tr></thead><tbody>'
                    + "".join(
                        f'<tr><td class="num">{i}</td><td class="num">{start_pos}</td><td class="num">{end_pos}</td><td>{html_escape(ref1_label if clade == "ia" else ref2_label)}</td></tr>'
                        for i, (start_pos, end_pos, clade, n_snps) in enumerate(merged_tracts, start=1)
                    )
                    + "</tbody></table>"
                )
            breakpoints_sections_html += (
                f'<div class="{section_cls}" data-row="{ri}" data-recombinant="{html_escape(rec_call)}">'
                f'<div class="rec-sites-row">'
                f'<span class="rec-sites-sample-id" title="{html_escape(sample_id)}">{html_escape(sample_id)}</span>'
                f'<div class="strip-cell">{strip_display}</div>'
                f'</div>'
                f'<details class="rec-sites-details"><summary>{summary_text}</summary>'
                f'{details_content}</details></div>'
            )
        strip_ref1_name = ref1_label if ref1_label not in ("ref1", "ref2") else (ref1_spec or ref1_label)
        strip_ref2_name = ref2_label if ref2_label not in ("ref1", "ref2") else (ref2_spec or ref2_label)
        breakpoints_section_html = (
            '<details class="collapsible-section diagnostic-strips-chart" open id="breakpointsStripsSection">'
            f'<summary><h2>Recombination tracts and predicted breakpoints per sample</h2><button class="pdf-btn" onclick="event.stopPropagation();exportStripSvg(\'breakpointsStripsSection\',\'breakpointsStripsContainer\',\'Recombination tracts and predicted breakpoints per sample\',\'{html_escape(ref1_label)}\',\'{html_escape(ref2_label)}\')">&#8595; Download SVG</button></summary>'
            '<div class="section-inner chart-section">'
            '<p class="threshold-note">Each coloured tract spans from the <strong>first to the last diagnostic SNP</strong> unambiguously derived from that clade. The predicted breakpoint lies somewhere in the <strong>uncoloured gap</strong> between adjacent tracts — its exact position cannot be determined because those intervening regions lack clade-informative diagnostic SNPs. Minimum consecutive diagnostic SNPs per tract: <strong>{min_consecutive}</strong>. <span id="breakpointsFilterCount" aria-live="polite"></span></p>'
            '<p class="threshold-note">{ref1} (blue), {ref2} (orange). Grey gaps = predicted breakpoint region (may be widened by ambiguous bases or poorly sequenced areas).</p>'
            '<div class="strip-legend"><span class="strip-legend-ia"></span> {ref1} &nbsp; <span class="strip-legend-ib"></span> {ref2} &nbsp; <span class="strip-legend-gap"></span> predicted breakpoint region (affected by ambiguous bases / poor coverage)</div>'
            '<div class="strip-strips-container" id="breakpointsStripScrollWrapper">'
            '<div id="breakpointsStripsContainer">'
        ).format(ref1=html_escape(strip_ref1_name), ref2=html_escape(strip_ref2_name), min_consecutive=min_consecutive)
        breakpoints_section_html = breakpoints_section_html + breakpoints_sections_html + "</div></div></div></details>"

    # Figure: diagnostic SNP positions on genome – count, histogram, then ruler (exact positions)
    snp_positions_section_html = ""
    snp_histogram_json = "[]"
    if diagnostic_snp_positions and genome_length:
        n_snps = len(diagnostic_snp_positions)
        hist_labels, hist_counts = _snp_positions_histogram_bins(
            diagnostic_snp_positions, genome_length, num_bins=60
        )
        snp_histogram_json = json.dumps({"labels": hist_labels, "counts": hist_counts})
        ruler_svg = _snp_positions_svg(diagnostic_snp_positions, display_length)
        snp_positions_table_rows = "".join(
            f'<tr><td class="num">{i}</td><td class="num">{pos}</td></tr>'
            for i, pos in enumerate(diagnostic_snp_positions, start=1)
        )
        snp_positions_section_html = (
            '<details class="collapsible-section" open id="snpPositionsSection">'
            '<summary><h2>Diagnostic SNP positions between reference genomes</h2>'
            '<span style="display:flex;gap:6px;margin-left:auto;flex-shrink:0;">'
            '<button class="pdf-btn" style="margin-left:0" onclick="event.stopPropagation();exportChartPng(\'chartSnpPositions\',\'snp_positions_histogram.png\')">&#8595; Histogram PNG</button>'
            '<button class="pdf-btn" style="margin-left:0" onclick="event.stopPropagation();exportSnpRulerSvg()">&#8595; Ruler SVG</button>'
            '</span></summary>'
            '<div class="section-inner chart-section">'
            '<p class="threshold-note"><strong>{n_snps} diagnostic SNPs.</strong> Density of diagnostic SNPs along the reference (alignment coordinates). Use this to interpret where recombination breakpoints may fall. Squirrel always builds alignments relative to reference NC_003310 (Clade I) or NC_063383 (Clade II), not the refs you specified.</p>'
            '<div class="snp-positions-wrapper"><div class="chart-container" style="height:220px;"><canvas id="chartSnpPositions"></canvas></div></div>'
            '<h3 class="snp-ruler-title">Exact positions</h3>'
            '<p class="threshold-note">Each tick marks one diagnostic SNP position along the genome (0 to {display_length} bp).</p>'
            '<div class="snp-positions-wrapper snp-ruler-wrapper">' + ruler_svg + '</div>'
            '<details class="rec-sites-details"><summary>Show diagnostic site table</summary>'
            '<table class="rec-sites-table"><thead><tr><th>#</th><th>Position (bp)</th></tr></thead><tbody>'
            + snp_positions_table_rows +
            '</tbody></table></details>'
            '</div></details>'
        ).format(n_snps=n_snps, genome_length=genome_length, display_length=display_length)

    phylogeny_section_html = ""
    if phylogeny_svg_content or phylogeny_pdf_base64:
        # Prefer inline SVG (no PDF viewer frame); fall back to PDF embed if no SVG
        svg_safe = (phylogeny_svg_content or "").replace("</script>", "<\\/script>")
        pdf_data_uri = ("data:application/pdf;base64," + phylogeny_pdf_base64) if phylogeny_pdf_base64 else ""
        download_btn = (
            f'<a class="pdf-btn" href="{pdf_data_uri}" download="phylogeny_tree.pdf" onclick="event.stopPropagation()">&#8595; Download PDF</a>'
            if pdf_data_uri else ""
        )
        if svg_safe:
            body_content = '<div class="phylogeny-svg-container">' + svg_safe + "</div>"
        elif pdf_data_uri:
            body_content = f'<embed src="{pdf_data_uri}" type="application/pdf" width="100%" height="800" style="max-width:100%; border-radius:8px; border:1px solid #dee2e6;" />'
        else:
            body_content = ""
        phylogeny_section_html = (
            '<details class="collapsible-section" open id="phylogenySection">'
            "<summary><h2>Phylogeny of recombinant ancestors</h2>"
            + download_btn
            + "</summary>"
            '<div class="section-inner chart-section">'
            + body_content
            + "</div></details>"
        )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>RecMpox Results</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/xlsx@0.18.5/dist/xlsx.full.min.js"></script>
<style>
* {{ margin: 0; padding: 0; box-sizing: border-box; }}
body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); padding: 40px 20px; min-height: 100vh; }}
.container {{ max-width: 1400px; margin: 0 auto; background: white; border-radius: 12px; box-shadow: 0 20px 60px rgba(0,0,0,0.3); overflow: hidden; }}
.header {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 30px; text-align: center; }}
.header h1 {{ font-size: 2em; margin-bottom: 8px; }}
.header p {{ opacity: 0.9; font-size: 1rem; }}
.content {{ padding: 30px; }}
.summary {{ background: #f8f9fa; padding: 20px; margin-bottom: 25px; border-radius: 8px; border-left: 4px solid #667eea; }}
.generated-note {{ font-size: 0.85em; color: #6c757d; margin-bottom: 10px; }}
.generated-note code {{ background: #e9ecef; padding: 1px 4px; border-radius: 3px; font-size: 0.95em; }}
.summary-stats {{ display: flex; gap: 20px; flex-wrap: wrap; margin-top: 10px; }}
.threshold-note {{ margin-top: 12px; color: #495057; font-size: 0.95em; }}
.filter-applies-note {{ margin-top: 6px; color: #6c757d; font-size: 0.9em; }}
.stat-box {{ background: white; padding: 15px 25px; border-radius: 6px; box-shadow: 0 2px 4px rgba(0,0,0,0.08); }}
.stat-box .number {{ font-size: 1.8em; font-weight: bold; color: #667eea; }}
.stat-box .label {{ color: #6c757d; font-size: 0.9em; margin-top: 4px; }}
.chart-section {{ margin-bottom: 30px; background: #fff; padding: 20px; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.08); }}
.phylogeny-svg-container {{ margin-top: 10px; overflow-x: auto; max-width: 100%; border-radius: 8px; border: 1px solid #dee2e6; }}
.phylogeny-svg-container svg {{ max-width: 100%; height: auto; display: block; }}
.chart-legend {{ display: flex; align-items: center; gap: 8px; flex-wrap: wrap; margin-top: 8px; margin-bottom: 10px; font-size: 0.9em; color: #495057; }}
.chart-legend-item {{ display: inline-block; width: 14px; height: 14px; border-radius: 2px; vertical-align: middle; }}
.snp-positions-wrapper {{ margin-top: 10px; overflow-x: auto; }}
.snp-positions-svg {{ display: block; min-width: 400px; }}
.snp-positions-svg .snp-tick-hit {{ cursor: pointer; }}
.snp-ruler-title {{ margin-top: 24px; margin-bottom: 6px; color: #333; font-size: 1em; }}
.snp-ruler-wrapper {{ margin-top: 6px; margin-bottom: 8px; }}
.diagnostic-strips-chart {{ max-width: 100%; min-width: 0; overflow: hidden; box-sizing: border-box; }}
.chart-section h2 {{ margin-bottom: 15px; color: #333; font-size: 1.2em; }}
.chart-wrapper {{ overflow-x: auto; overflow-y: hidden; margin-top: 10px; border-radius: 6px; }}
.chart-inner {{ position: relative; min-height: 320px; }}
.chart-container {{ position: relative; width: 100%; height: 320px; }}
.table-section {{ overflow-x: auto; margin-top: 20px; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.08); }}
.table-header {{ background: #f8f9fa; padding: 15px 20px; border-bottom: 2px solid #dee2e6; display: flex; justify-content: space-between; align-items: center; }}
.table-header h2 {{ margin: 0; color: #333; font-size: 1.2em; }}
table {{ width: 100%; border-collapse: collapse; background: white; }}
th {{ background: #f8f9fa; padding: 12px 15px; text-align: left; font-weight: 600; color: #495057; border-bottom: 2px solid #dee2e6; position: sticky; top: 0; z-index: 1; }}
th.sortable {{ cursor: pointer; user-select: none; }}
th.sortable:hover {{ background: #e9ecef; }}
th.sortable::after {{ content: ' ↕'; opacity: 0.5; font-size: 0.8em; }}
th.sortable.asc::after {{ content: ' ↑'; opacity: 1; }}
th.sortable.desc::after {{ content: ' ↓'; opacity: 1; }}
td {{ padding: 12px 15px; border-bottom: 1px solid #e9ecef; }}
td.num {{ text-align: right; }}
tr:hover {{ background: #f8f9fa; }}
tr.hidden {{ display: none; }}
.rec-sites-section.hidden {{ display: none; }}
tr.highlight {{ background: #fff3cd !important; }}
tr.recombinant {{ }}
#filterrow {{ background: #ecf0f1; }}
#filterrow .filter-label {{ font-weight: 600; color: #495057; }}
#filterrow input {{ width: 70px; padding: 6px; border: 1px solid #bdc3c7; border-radius: 4px; }}
#filterrow .rec-filter-select {{ padding: 6px 8px; border: 1px solid #bdc3c7; border-radius: 4px; min-width: 120px; }}
.recombinant-badge {{ font-size: 0.85em; padding: 2px 8px; border-radius: 4px; font-weight: 500; }}
.recombinant-badge.recombinant {{ background: #f8d7da; color: #721c24; }}
.recombinant-badge.potential {{ background: #fff3cd; color: #856404; }}
.recombinant-badge.no {{ background: #d4edda; color: #155724; }}
.accession-link {{ color: #1976d2; text-decoration: none; font-weight: 500; }}
.accession-link:hover {{ text-decoration: underline; }}
#diagnosticStripsContainer {{ width: 100%; max-width: 100%; min-width: 0; overflow: hidden; box-sizing: border-box; }}
.strip-strips-container {{ width: 100%; max-width: 100%; min-width: 0; margin-top: 8px; margin-bottom: 12px; border-radius: 6px; border: 1px solid #dee2e6; padding: 8px; overflow: hidden; box-sizing: border-box; }}
.rec-sites-section {{ margin-bottom: 16px; padding: 10px; background: #fafafa; border-radius: 6px; border-left: 3px solid #dee2e6; width: 100%; max-width: 100%; min-width: 0; overflow: hidden; box-sizing: border-box; }}
.rec-sites-section.recombinant {{ border-left-color: #E89B3C; }}
.rec-sites-row {{ display: flex; align-items: center; gap: 12px; flex-wrap: nowrap; min-width: 0; width: 100%; max-width: 100%; }}
.rec-sites-sample-id {{ font-size: 0.95em; font-weight: 600; color: #333; width: 300px; min-width: 300px; max-width: 300px; overflow: visible; white-space: normal; word-break: break-word; flex-shrink: 0; }}
.strip-cell {{ flex: 1 0 0; min-width: 0; overflow-x: auto; overflow-y: hidden; border-radius: 4px; border: 1px solid #e9ecef; -webkit-overflow-scrolling: touch; }}
.strip-genome {{ position: relative; width: 100%; height: 24px; min-width: 200px; border-radius: 4px; overflow: hidden; background: #e9ecef; }}
.strip-genome.breakpoints-strip {{ background: linear-gradient(135deg, #c4d3e0 0%, #b6c8d7 60%, #bfcdd9 100%); height: 32px; border-radius: 6px; box-shadow: inset 0 2px 8px rgba(74,144,217,0.13), inset 0 -1px 3px rgba(0,0,0,0.07); }}
.strip-genome.breakpoints-strip::after {{ content: ''; position: absolute; inset: 0; background: linear-gradient(to bottom, rgba(255,255,255,0.22) 0%, transparent 50%); pointer-events: none; z-index: 5; border-radius: inherit; }}
.strip-segment {{ position: absolute; top: 0; height: 100%; width: 2px; transition: opacity 0.15s; }}
.strip-segment:hover {{ opacity: 0.75; }}
.strip-segment.ia {{ background: #4A90D9; }}
.strip-segment.ib {{ background: #E89B3C; }}
.strip-segment.other {{ background: #95a5a6; }}
.strip-segment.other-n {{ background: #95a5a6; }}
.strip-segment.segment-hidden-by-filter {{ display: none; }}
.strip-segment.region-segment {{ min-width: 4px; border-radius: 3px; }}
.strip-segment.breakpoint-marker {{ width: 4px; background: #c0392b; transform: translateX(-50%); }}
#breakpointsStripsContainer .rec-sites-section {{ background: linear-gradient(135deg, #fafbfc 0%, #f4f6f9 100%); border-left: 4px solid #dee2e6; box-shadow: 0 1px 4px rgba(0,0,0,0.05); transition: box-shadow 0.15s; }}
#breakpointsStripsContainer .rec-sites-section.recombinant {{ border-left-color: #E89B3C; }}
#breakpointsStripsContainer .rec-sites-section:hover {{ box-shadow: 0 3px 10px rgba(0,0,0,0.10); }}
.strip-legend {{ display: flex; align-items: center; gap: 4px; flex-wrap: wrap; margin-bottom: 12px; font-size: 0.9em; color: #495057; }}
.strip-legend-ia {{ display: inline-block; width: 14px; height: 14px; background: #4A90D9; border-radius: 2px; }}
.strip-legend-ib {{ display: inline-block; width: 14px; height: 14px; background: #E89B3C; border-radius: 2px; }}
.strip-legend-other {{ display: inline-block; width: 14px; height: 14px; background: #95a5a6; border-radius: 2px; }}
.strip-legend-gap {{ display: inline-block; width: 14px; height: 14px; background: linear-gradient(135deg, #c4d3e0, #b6c8d7); border: 1px solid #9fb8cc; border-radius: 2px; }}
.strip-classification-filter-wrap {{ display: flex; flex-wrap: wrap; align-items: center; justify-content: flex-end; gap: 12px 20px; margin-bottom: 12px; padding: 10px 12px; background: #f8f9fa; border-radius: 6px; border: 1px solid #e9ecef; }}
.strip-classification-filter-wrap label {{ font-size: 0.95em; color: #495057; cursor: pointer; white-space: nowrap; }}
.strip-classification-filter-wrap input[type="checkbox"] {{ margin-right: 6px; vertical-align: middle; }}
.strip-classification-options {{ display: none; gap: 12px 20px; flex-wrap: wrap; align-items: center; }}
.strip-classification-options[aria-hidden="false"] {{ display: flex; }}
.strip-filter-check-label {{ font-weight: 600; }}
.strip-ruler {{ position: relative; width: 100%; height: 22px; min-width: 200px; margin-top: 3px; }}
.ruler-tick {{ position: absolute; transform: translateX(-50%); font-size: 0.67em; font-weight: 500; color: #6c757d; white-space: nowrap; line-height: 1; padding-top: 6px; letter-spacing: 0.01em; }}
.ruler-tick::before {{ content: ''; display: block; position: absolute; top: 0; left: 50%; transform: translateX(-50%); width: 1px; height: 5px; background: #adb5bd; }}
.rec-sites-details {{ margin-top: 10px; font-size: 0.9em; }}
.rec-sites-details summary {{ cursor: pointer; color: #667eea; font-weight: 500; }}
.rec-sites-table {{ margin-top: 8px; border-collapse: collapse; font-size: 0.9em; max-height: 200px; overflow: auto; }}
.rec-sites-table th, .rec-sites-table td {{ padding: 6px 10px; border: 1px solid #dee2e6; text-align: left; }}
.rec-sites-table th {{ background: #f8f9fa; }}
.rec-sites-table .num {{ text-align: right; }}
.collapsible-section {{ margin-bottom: 20px; border: 1px solid #dee2e6; border-radius: 8px; overflow: hidden; }}
.collapsible-section summary {{ cursor: pointer; padding: 15px 20px; background: #f8f9fa; font-weight: 600; color: #333; list-style: none; display: flex; align-items: center; }}
.collapsible-section summary::-webkit-details-marker {{ display: none; }}
.collapsible-section summary::before {{ content: '▼'; font-size: 0.7em; margin-right: 10px; transition: transform 0.2s; }}
.collapsible-section:not([open]) summary::before {{ transform: rotate(-90deg); }}
.collapsible-section summary h2 {{ margin: 0; font-size: 1.2em; }}
.collapsible-section .section-inner {{ padding: 20px; }}
.pdf-btn {{ margin-left: auto; flex-shrink: 0; font-size: 0.78em; padding: 5px 13px; border: 1.5px solid #667eea; border-radius: 5px; background: white; color: #667eea; cursor: pointer; font-weight: 600; transition: background 0.15s, color 0.15s; white-space: nowrap; }}
.pdf-btn:hover {{ background: #667eea; color: white; }}
.pdf-btn:disabled {{ opacity: 0.6; cursor: default; }}
</style>
</head>
<body>
<!-- Generated: {gen_time} -->
<div class="container">
<div class="header">
<h1>RecMpox Results</h1>
<p>Tool to flag potential recombinant mpox genomes.</p>
</div>
<div class="content">
<div class="summary">
<h2>Summary</h2>
<p class="generated-note">Report generated: {gen_time}</p>
{snps_only_note_html}
<div class="summary-stats">
<div class="stat-box"><div class="number">{len(results)}</div><div class="label">Genomes (this file)</div></div>
{sites_breakdown_html}
{refs_summary_html}
{part_html}
</div>
{threshold_html}
</div>
<details class="collapsible-section" open>
<summary><h2>Per-genome classification</h2><button class="pdf-btn" onclick="event.stopPropagation();exportTableXlsx()">&#8595; Download XLSX</button></summary>
<div class="section-inner table-section">
<table id="t">
<thead>
{thead}
{filter_row}
</thead>
<tbody>
{tbody}
</tbody>
</table>
</div>
</details>
<details class="collapsible-section" open>
<summary><h2>Diagnostic SNPs per genome</h2><button class="pdf-btn" onclick="event.stopPropagation();exportChartPng(\'chartBar\',\'diagnostic_snps_barplot.png\')">&#8595; Download PNG</button></summary>
<div class="section-inner chart-section">
<p class="threshold-note">Stacked percentage per genome: % {html_escape(ref1_label)} (blue), % {html_escape(ref2_label)} (purple), % other (gray).</p>
<div class="chart-legend stacked-bar-legend">
<span class="chart-legend-item" style="background:#667eea;"></span> % {html_escape(ref1_label)} &nbsp;
<span class="chart-legend-item" style="background:#764ba2;"></span> % {html_escape(ref2_label)} &nbsp;
<span class="chart-legend-item" style="background:#c8c8c8;"></span> % other
</div>
<div class="chart-wrapper">
<div id="chartInner" class="chart-inner">
<div class="chart-container"><canvas id="chartBar"></canvas></div>
</div>
</div>
</div>
</details>
{snp_positions_section_html}
{rec_sites_html}
{breakpoints_section_html}
{phylogeny_section_html}
</div>
</div>
<script>
function exportChartPng(canvasId, filename) {{
  var canvas = document.getElementById(canvasId);
  if (!canvas) {{ alert('Chart not ready yet – please wait a moment and try again.'); return; }}
  var a = document.createElement('a');
  a.href = canvas.toDataURL('image/png');
  a.download = filename;
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
}}

function exportSnpRulerSvg() {{
  var svg = document.querySelector('.snp-positions-svg');
  if (!svg) {{ alert('SVG ruler not found.'); return; }}
  var serializer = new XMLSerializer();
  var svgStr = serializer.serializeToString(svg);
  var blob = new Blob([svgStr], {{ type: 'image/svg+xml;charset=utf-8' }});
  var a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = 'diagnostic_snp_positions_ruler.svg';
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
  URL.revokeObjectURL(a.href);
}}

function exportTableXlsx() {{
  if (typeof XLSX === 'undefined') {{ alert('XLSX library not loaded.'); return; }}
  var headers = [];
  var thead = document.querySelector('#t thead tr:first-child');
  if (!thead) return;
  thead.querySelectorAll('th').forEach(function(th) {{
    headers.push(th.textContent.trim().replace(/[\u2195\u2191\u2193\u25b2\u25bc]/g, '').trim());
  }});
  var tbody = document.querySelector('#t tbody');
  var wsData = [headers];
  Array.from(tbody.querySelectorAll('tr:not(.hidden)')).forEach(function(tr) {{
    var row = [];
    Array.from(tr.cells).forEach(function(td) {{
      var val = td.textContent.trim();
      var num = parseFloat(val);
      row.push(isNaN(num) ? val : num);
    }});
    wsData.push(row);
  }});
  var ws = XLSX.utils.aoa_to_sheet(wsData);
  var wb = XLSX.utils.book_new();
  XLSX.utils.book_append_sheet(wb, ws, 'RecMpox Results');
  XLSX.writeFile(wb, 'recmpox_results.xlsx');
}}

function escSvg(s) {{
  return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;');
}}

function exportStripSvg(sectionId, containerId, figTitle, ref1Label, ref2Label) {{
  var section = document.getElementById(sectionId);
  if (section) section.open = true;
  var container = document.getElementById(containerId);
  if (!container) return;

  var rows = Array.from(container.querySelectorAll('.rec-sites-section:not(.hidden)'));
  if (!rows.length) {{ alert('No visible samples to export.'); return; }}

  // Layout constants (px → SVG user units)
  var ml = 16, mr = 20;
  var idW = 190;       // sample-ID column
  var stripW = 1100;   // genome strip width
  var W = ml + idW + stripW + mr;
  var rowH = 28;
  var rowGap = 5;
  var rulerH = 22;
  var titleH = 22;
  var legendH = 26;
  var topPad = titleH + 8;
  var rowsH = rows.length * (rowH + rowGap) - rowGap;
  var rulerLineY = topPad + rowsH + 6;
  var rulerBaseY = rulerLineY + 18;
  var legY = rulerBaseY + 22;
  var H = legY + 14;

  var colorMap = {{ ia: '#4A90D9', ib: '#E89B3C', other: '#95a5a6' }};
  var isBp = !!container.querySelector('.breakpoints-strip');
  var stripBg = isBp ? '#bfcdd9' : '#e9ecef';

  var s = [];
  s.push('<?xml version="1.0" encoding="UTF-8"?>');
  s.push('<svg xmlns="http://www.w3.org/2000/svg" width="' + W + '" height="' + H + '">');
  s.push('<rect width="100%" height="100%" fill="white"/>');

  // Title
  s.push('<text x="' + ml + '" y="' + (titleH - 4) + '" font-family="Arial,sans-serif" font-size="13" font-weight="bold" fill="#333">' + escSvg(figTitle) + '</text>');

  // Sample rows
  rows.forEach(function(row, ri) {{
    var ry = topPad + ri * (rowH + rowGap);

    // Sample ID (use title attr for full name if truncated)
    var idEl = row.querySelector('.rec-sites-sample-id');
    var sid = idEl ? (idEl.getAttribute('title') || idEl.textContent.trim()) : '';
    if (sid.length > 28) sid = sid.slice(0,27) + '\u2026';
    s.push('<text x="' + (ml+idW-4) + '" y="' + (ry+rowH/2) + '" font-family="Arial,sans-serif" font-size="10" fill="#333" text-anchor="end" dominant-baseline="middle">' + escSvg(sid) + '</text>');

    // Strip background
    s.push('<rect x="' + (ml+idW) + '" y="' + ry + '" width="' + stripW + '" height="' + rowH + '" fill="' + stripBg + '" rx="3"/>');

    // Segments
    var genome = row.querySelector('.strip-genome');
    if (genome) {{
      genome.querySelectorAll('.strip-segment').forEach(function(seg) {{
        var cls = seg.classList;
        var color = cls.contains('ia') ? colorMap.ia : cls.contains('ib') ? colorMap.ib : cls.contains('other') ? colorMap.other : null;
        if (!color) return;
        var lp = parseFloat(seg.style.left);
        if (isNaN(lp)) return;
        var wp = seg.style.width ? parseFloat(seg.style.width) : 0;
        var sx = ml + idW + (lp/100) * stripW;
        var sw = wp ? Math.max(2, (wp/100)*stripW) : 2;
        s.push('<rect x="' + sx.toFixed(1) + '" y="' + ry + '" width="' + sw.toFixed(1) + '" height="' + rowH + '" fill="' + color + '"/>');
      }});
      // Strip border
      s.push('<rect x="' + (ml+idW) + '" y="' + ry + '" width="' + stripW + '" height="' + rowH + '" fill="none" stroke="#dee2e6" stroke-width="0.5" rx="3"/>');
    }}
  }});

  // Ruler – below the sample rows, read ticks from first visible row
  var firstRuler = rows[0].querySelector('.strip-ruler');
  if (firstRuler) {{
    s.push('<line x1="' + (ml+idW) + '" y1="' + rulerLineY + '" x2="' + (ml+idW+stripW) + '" y2="' + rulerLineY + '" stroke="#ccc" stroke-width="0.5"/>');
    firstRuler.querySelectorAll('.ruler-tick').forEach(function(tick) {{
      var lp = parseFloat(tick.style.left);
      if (isNaN(lp)) return;
      var tx = ml + idW + (lp/100) * stripW;
      s.push('<line x1="' + tx + '" y1="' + rulerLineY + '" x2="' + tx + '" y2="' + (rulerLineY+6) + '" stroke="#adb5bd" stroke-width="1"/>');
      s.push('<text x="' + tx + '" y="' + rulerBaseY + '" font-family="Arial,sans-serif" font-size="9" fill="#777" text-anchor="middle">' + escSvg(tick.textContent.trim()) + '</text>');
    }});
  }}

  // Legend
  var legX = ml + idW;
  function legItem(x, color, label, isRect) {{
    if (isRect) s.push('<rect x="'+x+'" y="'+(legY-11)+'" width="13" height="13" fill="'+color+'" stroke="#adb5bd" stroke-width="0.5" rx="2"/>');
    else s.push('<rect x="'+x+'" y="'+(legY-11)+'" width="13" height="13" fill="'+color+'" rx="2"/>');
    s.push('<text x="'+(x+17)+'" y="'+legY+'" font-family="Arial,sans-serif" font-size="10" fill="#495057">'+escSvg(label)+'</text>');
  }}
  legItem(legX, '#4A90D9', ref1Label, false); legX += 14+8+ref1Label.length*6+12;
  legItem(legX, '#E89B3C', ref2Label, false); legX += 14+8+ref2Label.length*6+12;
  if (!isBp) {{
    legItem(legX, '#95a5a6', 'other', false);
  }} else {{
    legItem(legX, '#bfcdd9', 'predicted breakpoint region (affected by ambiguous bases / poor coverage)', true);
  }}

  s.push('</svg>');

  var blob = new Blob([s.join('\\n')], {{ type: 'image/svg+xml;charset=utf-8' }});
  var a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = figTitle.replace(/[^a-z0-9]+/gi,'_').toLowerCase() + '.svg';
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
  URL.revokeObjectURL(a.href);
}}
(function() {{
  function run() {{
    if (typeof Chart === "undefined") {{ setTimeout(run, 30); return; }}
  var data = {data_json};
  var snpHistogramData = {snp_histogram_json};
  var table = document.getElementById("t");
  var thead = table.querySelector("thead tr:first-child");
  var tbody = table.querySelector("tbody");
  var rows = Array.from(tbody.querySelectorAll("tr"));
  var sortDir = 1;
  var lastCol = -1;

  function numVal(cell) {{
    var t = cell && cell.textContent ? cell.textContent.trim() : "";
    var n = parseFloat(t);
    return isNaN(n) ? (t || 0) : n;
  }}

  function applyFilters() {{
    var inputs = table.querySelectorAll(".min-thresh");
    var recFilter = document.getElementById("recFilter");
    var recVal = recFilter ? recFilter.value : "all";
    rows.forEach(function(tr) {{
      var show = true;
      inputs.forEach(function(inp) {{
        var col = parseInt(inp.dataset.col, 10);
        var minV = parseFloat(inp.value);
        if (!isNaN(minV) && inp.value !== "") {{
          var cellVal = numVal(tr.cells[col]);
          if (typeof cellVal !== "number" || cellVal < minV) show = false;
        }}
      }});
      if (recVal !== "all") {{
        var dr = tr.getAttribute("data-recombinant");
        if (dr !== recVal) show = false;
      }}
      tr.classList.toggle("hidden", !show);
    }});
    var visibleRows = new Set();
    rows.forEach(function(tr) {{
      if (!tr.classList.contains("hidden")) {{
        var ri = tr.getAttribute("data-row");
        if (ri !== null) visibleRows.add(ri);
      }}
    }});
    var stripFilterCheck = document.getElementById("stripClassificationFilterCheck");
    var stripOpts = document.querySelectorAll(".strip-class-opt:checked");
    // When "Filter strips by classification" is checked: figure is empty until you pick Ia, Ib, etc.
    // When unchecked: show everything.
    var stripClassificationOn = stripFilterCheck && stripFilterCheck.checked;
    function rowMatchesStripClassification(ri) {{
      if (!stripClassificationOn || !data[ri]) return true;
      var r = data[ri];
      var n_ia = (r.n_ia != null) ? (parseInt(r.n_ia, 10) || 0) : 0;
      var n_ib = (r.n_ib != null) ? (parseInt(r.n_ib, 10) || 0) : 0;
      var n_other = (r.n_other != null) ? (parseInt(r.n_other, 10) || 0) : 0;
      var n_other_n = (r.n_other_n != null) ? (parseInt(r.n_other_n, 10) || 0) : 0;
      var hasOtherAll = n_other > 0;
      var hasOtherN = n_other_n > 0;
      for (var i = 0; i < stripOpts.length; i++) {{
        var f = stripOpts[i].getAttribute("data-filter");
        if (f === "ia" && n_ia > 0) return true;
        if (f === "ib" && n_ib > 0) return true;
        if (f === "other_all" && hasOtherAll) return true;
        if (f === "other_n" && hasOtherN) return true;
      }}
      return false;
    }}
    var stripVisibleRows = new Set();
    visibleRows.forEach(function(ri) {{
      if (rowMatchesStripClassification(ri)) stripVisibleRows.add(ri);
    }});
    // When filter is on but no option selected, show nothing; when options selected, show matching rows
    var rowsForStrips = !stripClassificationOn ? visibleRows : (stripOpts.length > 0 ? stripVisibleRows : new Set());
    var container = document.getElementById("diagnosticStripsContainer");
    if (container) {{
      var sections = container.querySelectorAll(".rec-sites-section");
      var showIa = false, showIb = false, showOtherAll = false, showOtherN = false;
      for (var i = 0; i < stripOpts.length; i++) {{
        var f = stripOpts[i].getAttribute("data-filter");
        if (f === "ia") showIa = true;
        if (f === "ib") showIb = true;
        if (f === "other_all") showOtherAll = true;
        if (f === "other_n") showOtherN = true;
      }}
      sections.forEach(function(sec) {{
        var ri = sec.getAttribute("data-row");
        var sectionVisible = rowsForStrips.has(ri);
        sec.classList.toggle("hidden", !sectionVisible);
        var stripGenome = sec.querySelector(".strip-genome:not(.breakpoints-strip)");
        if (stripGenome) {{
          stripGenome.querySelectorAll(".strip-segment").forEach(function(seg) {{
            if (!stripClassificationOn) {{
              seg.classList.remove("segment-hidden-by-filter");
            }} else if (sectionVisible) {{
              var hide = false;
              if (seg.classList.contains("ia")) hide = !showIa;
              else if (seg.classList.contains("ib")) hide = !showIb;
              else if (seg.classList.contains("other-n")) hide = !showOtherAll && !showOtherN;
              else if (seg.classList.contains("other")) hide = !showOtherAll;
              seg.classList.toggle("segment-hidden-by-filter", hide);
            }}
          }});
        }}
      }});
    }}
    var bpContainer = document.getElementById("breakpointsStripsContainer");
    if (bpContainer) {{
      var bpSections = bpContainer.querySelectorAll(".rec-sites-section");
      bpSections.forEach(function(sec) {{
        var ri = sec.getAttribute("data-row");
        sec.classList.toggle("hidden", !rowsForStrips.has(ri));
      }});
    }}
    var countEl = document.getElementById("stripFilterCount");
    if (countEl) {{
      var n = rowsForStrips.size;
      var total = rows.length;
      if (stripClassificationOn)
        countEl.textContent = "Showing " + n + " of " + total + " samples (filtered by classification).";
      else if (n < total && total > 0)
        countEl.textContent = "Showing " + n + " of " + total + " samples.";
      else
        countEl.textContent = "";
    }}
    var bpCountEl = document.getElementById("breakpointsFilterCount");
    if (bpCountEl) {{
      var n = rowsForStrips.size;
      var total = rows.length;
      if (stripClassificationOn)
        bpCountEl.textContent = "Showing " + n + " of " + total + " samples (filtered by classification).";
      else if (n < total && total > 0)
        bpCountEl.textContent = "Showing " + n + " of " + total + " samples.";
      else
        bpCountEl.textContent = "";
    }}
    if (window.chartBar) updateChart();
  }}

  table.querySelectorAll(".min-thresh").forEach(function(inp) {{
    inp.addEventListener("input", applyFilters);
    inp.addEventListener("change", applyFilters);
  }});
  var recFilterEl = document.getElementById("recFilter");
  if (recFilterEl) {{
    recFilterEl.addEventListener("change", applyFilters);
    recFilterEl.addEventListener("input", applyFilters);
  }}
  var stripFilterCheckEl = document.getElementById("stripClassificationFilterCheck");
  if (stripFilterCheckEl) {{
    stripFilterCheckEl.addEventListener("change", function() {{
      var opts = document.getElementById("stripClassificationOptions");
      if (opts) opts.setAttribute("aria-hidden", this.checked ? "false" : "true");
      applyFilters();
    }});
  }}
  document.querySelectorAll(".strip-class-opt").forEach(function(el) {{
    el.addEventListener("change", function() {{
      // If any option is turned on, auto-open + check the master filter to make it obvious
      var anyChecked = document.querySelector(".strip-class-opt:checked") != null;
      var master = document.getElementById("stripClassificationFilterCheck");
      var opts = document.getElementById("stripClassificationOptions");
      if (master) master.checked = anyChecked;
      if (opts) opts.setAttribute("aria-hidden", anyChecked ? "false" : "true");
      applyFilters();
    }});
  }});

  thead.querySelectorAll("th.sortable").forEach(function(th) {{
    var colIndex = parseInt(th.dataset.col, 10);
    th.addEventListener("click", function() {{
      sortDir = (lastCol === colIndex) ? -sortDir : 1;
      lastCol = colIndex;
      thead.querySelectorAll("th").forEach(function(h) {{ h.classList.remove("asc", "desc"); }});
      th.classList.add(sortDir === 1 ? "asc" : "desc");
      rows.sort(function(a, b) {{
        var av = numVal(a.cells[colIndex]);
        var bv = numVal(b.cells[colIndex]);
        if (av < bv) return -sortDir;
        if (av > bv) return sortDir;
        return 0;
      }});
      rows.forEach(function(r) {{ tbody.appendChild(r); }});
      var stripContainer = document.getElementById("diagnosticStripsContainer");
      if (stripContainer) {{
        var order = rows.map(function(r) {{ return r.getAttribute("data-row"); }});
        var byRow = {{}};
        stripContainer.querySelectorAll(".rec-sites-section").forEach(function(sec) {{
          var ri = sec.getAttribute("data-row");
          if (ri !== null) byRow[ri] = sec;
        }});
        order.forEach(function(ri) {{
          if (byRow[ri]) stripContainer.appendChild(byRow[ri]);
        }});
      }}
      var bpStripContainer = document.getElementById("breakpointsStripsContainer");
      if (bpStripContainer) {{
        var order = rows.map(function(r) {{ return r.getAttribute("data-row"); }});
        var byRow = {{}};
        bpStripContainer.querySelectorAll(".rec-sites-section").forEach(function(sec) {{
          var ri = sec.getAttribute("data-row");
          if (ri !== null) byRow[ri] = sec;
        }});
        order.forEach(function(ri) {{
          if (byRow[ri]) bpStripContainer.appendChild(byRow[ri]);
        }});
      }}
      if (window.chartBar) updateChart();
    }});
  }});

  function truncId(id) {{ return id && id.length > 20 ? id.slice(0, 20) + "\u2026" : (id || ""); }}
  var chartBar = null;
  function updateChart() {{
    var visibleRows = rows.filter(function(r) {{ return !r.classList.contains("hidden"); }});
    var labels = visibleRows.map(function(r) {{
      var ri = parseInt(r.getAttribute("data-row"), 10);
      return truncId((data[ri] && data[ri].id) ? data[ri].id : "");
    }});
    var pct1 = visibleRows.map(function(r) {{
      var ri = parseInt(r.getAttribute("data-row"), 10);
      return (data[ri] != null) ? (parseFloat(data[ri].pct_ia) || 0) : 0;
    }});
    var pct2 = visibleRows.map(function(r) {{
      var ri = parseInt(r.getAttribute("data-row"), 10);
      return (data[ri] != null) ? (parseFloat(data[ri].pct_ib) || 0) : 0;
    }});
    var pctOther = visibleRows.map(function(r) {{
      var ri = parseInt(r.getAttribute("data-row"), 10);
      return (data[ri] != null) ? (parseFloat(data[ri].pct_other) || 0) : 0;
    }});
    var chartInner = document.getElementById("chartInner");
    if (chartInner) chartInner.style.minWidth = Math.max(400, labels.length * 56) + "px";
    if (!chartBar) {{
      var ctx = document.getElementById("chartBar").getContext("2d");
      chartBar = new Chart(ctx, {{
        type: "bar",
        data: {{
          labels: labels,
          datasets: [
            {{ label: "% {html_escape(ref1_label)}", data: pct1, stack: "stack1", backgroundColor: "rgba(102,126,234,0.85)", borderColor: "rgba(102,126,234,1)", borderWidth: 1 }},
            {{ label: "% {html_escape(ref2_label)}", data: pct2, stack: "stack1", backgroundColor: "rgba(118,75,162,0.85)", borderColor: "rgba(118,75,162,1)", borderWidth: 1 }},
            {{ label: "% other", data: pctOther, stack: "stack1", backgroundColor: "rgba(200,200,200,0.85)", borderColor: "rgba(160,160,160,1)", borderWidth: 1 }}
          ]
        }},
        options: {{
          responsive: true,
          maintainAspectRatio: false,
          datasets: {{ bar: {{ barPercentage: 0.95, categoryPercentage: 0.95 }} }},
          onClick: function(ev, els) {{
            if (els.length && labels[els[0].index]) {{
              var id = labels[els[0].index];
              rows.forEach(function(r) {{ r.classList.remove("highlight"); }});
              for (var i = 0; i < rows.length; i++) {{
                if (rows[i].cells[0] && rows[i].cells[0].textContent.trim() === id) {{
                  rows[i].classList.add("highlight");
                  rows[i].scrollIntoView({{ block: "nearest", behavior: "smooth" }});
                  break;
                }}
              }}
            }}
          }},
          scales: {{
            x: {{ title: {{ display: true, text: "Sample ID" }}, ticks: {{ maxRotation: 45, minRotation: 45, autoSkip: false, font: {{ size: 11 }} }} }},
            y: {{ title: {{ display: true, text: "Percentage (%)" }}, min: 0, max: 100, ticks: {{ stepSize: 20 }} }}
          }},
          plugins: {{ legend: {{ display: true, position: "top" }} }}
        }}
      }});
      window.chartBar = chartBar;
    }} else {{
      chartBar.data.labels = labels;
      chartBar.data.datasets[0].data = pct1;
      chartBar.data.datasets[1].data = pct2;
      chartBar.data.datasets[2].data = pctOther;
      chartBar.update();
    }}
  }}

  rows.forEach(function(tr) {{
    tr.addEventListener("click", function() {{
      rows.forEach(function(r) {{ r.classList.remove("highlight"); }});
      tr.classList.add("highlight");
    }});
  }});

  if (data.length) updateChart();
  applyFilters();

  if (snpHistogramData && snpHistogramData.labels && snpHistogramData.labels.length && document.getElementById("chartSnpPositions")) {{
    var ctxSnp = document.getElementById("chartSnpPositions").getContext("2d");
    new Chart(ctxSnp, {{
      type: "bar",
      data: {{
        labels: snpHistogramData.labels,
        datasets: [{{ label: "Diagnostic SNPs", data: snpHistogramData.counts, backgroundColor: "rgba(102,126,234,0.7)", borderColor: "rgba(102,126,234,1)", borderWidth: 1 }}]
      }},
      options: {{
        responsive: true,
        maintainAspectRatio: false,
        indexAxis: "x",
        scales: {{
          x: {{ title: {{ display: true, text: "Position along genome (kb)" }}, ticks: {{ maxRotation: 45, autoSkip: true, maxTicksLimit: 20 }} }},
          y: {{ title: {{ display: true, text: "Number of SNPs" }}, beginAtZero: true, ticks: {{ stepSize: 1 }} }}
        }},
        plugins: {{ legend: {{ display: false }} }}
      }}
    }});
  }}
  }}
  run();
}})();
</script>
</body>
</html>
"""
    with open(out_path, "w") as f:
        f.write(html)


def split_fasta_into_seqs(fasta: Path, out_dir: Path) -> List[Path]:
    """Split multi-FASTA into one file per sequence. Returns list of paths."""
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: List[Path] = []
    current_id, current_lines = None, []
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                if current_id is not None:
                    p = out_dir / f"{current_id}.fa"
                    with open(p, "w") as out:
                        out.writelines(current_lines)
                    paths.append(p)
                current_id = line[1:].split()[0].strip().replace("/", "_")
                current_lines = [line]
            else:
                current_lines.append(line)
        if current_id is not None:
            p = out_dir / f"{current_id}.fa"
            with open(p, "w") as out:
                out.writelines(current_lines)
            paths.append(p)
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="recmpox",
        description=f"RecMpox v{__version__}: Flag potential recombination in mpox consensus genomes using diagnostic sites between two reference lineages.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
        epilog="""
Examples:
  recmpox -i fasta/ -o output -ref Ia,Ib
  recmpox -i OZ375330.1 -o output -ref Ib,IIb  # UK recombinant case example
  recmpox -i accessions.txt -o output -ref Ia,Ib

  # Override references: -ref1/-ref2 with optional -ref1_g/-ref2_g
  recmpox -i fasta/ -o output -ref1 NC_003310.1 -ref2 PP601219.1 -ref1_g Ia -ref2_g Ib
  recmpox -i fasta/ -o output -ref Ia,Ib -ref1 /path/to/custom_ia.fa

""",
    )
    required = parser.add_argument_group("required arguments (must specify when running)")
    optional = parser.add_argument_group("optional arguments")
    parser.add_argument("-h", "-help", "--help", action="help", help="show this help message and exit")
    optional.add_argument("-version", action="version", version=f"RecMpox v{__version__}")
    optional.add_argument(
        "-m",
        "-minor-ref-pct",
        dest="minor_ref_pct",
        type=float,
        default=MINOR_REF_PCT_THRESHOLD,
        metavar="",
        help=f"Minor reference %% threshold for calling 'potential recombinant' (default: {MINOR_REF_PCT_THRESHOLD:g})",
    )
    optional.add_argument(
        "-b",
        "-breakpoint-snp",
        dest="breakpoint_min_snps",
        action="store_const",
        const=2,
        default=1,
        help="Ignore single-SNP runs when inferring breakpoints (with -b: minimum consecutive diagnostic SNPs per tract = 2; default without -b: 1)",
    )
    required.add_argument("-i", "-input", dest="input", type=Path, default=None, metavar="", help="FASTA file, directory of .fa/.fasta/.fna, .txt file of accessions (one per line or comma-separated), NCBI accession, or comma-separated accessions (e.g. -i ACC1,ACC2 or -i accessions.txt)")
    required.add_argument("-ref", dest="ref", type=str, default=None, metavar="", help="Reference pair: two comma-separated clades among Ia, Ib, IIa, IIb (e.g. Ia,Ib or Ib,IIb). Builds one consensus per clade from earliest 5 genomes (Pathoplexus; requires Squirrel/mafft). Either -ref or both -ref1 and -ref2 are required")
    required.add_argument("-ref1", type=str, default=None, metavar="", help="First reference: FASTA path or NCBI accession; overrides ref1 when using -ref. Required if -ref is not used")
    required.add_argument("-ref2", type=str, default=None, metavar="", help="Second reference: FASTA path or NCBI accession; overrides ref2 when using -ref. Required if -ref is not used")
    optional.add_argument("-o", "-output", dest="output_dir", type=str, default="output", metavar="", help="Output directory (default: output); path is relative to cwd; always removed and recreated at start of each run")
    optional.add_argument("-ref1_g", type=str, default=None, metavar="", help="Genotype label for ref1 (TSV/HTML column headers; default from -ref or ref1 accession)")
    optional.add_argument("-ref2_g", type=str, default=None, metavar="", help="Genotype label for ref2 (TSV/HTML column headers; default from -ref or ref2 accession)")
    optional.add_argument("-include-indels", action="store_true", dest="include_indels", help="Include diagnostic indels (large indels) in addition to SNPs; by default only diagnostic SNPs are used")
    optional.add_argument("-min-indel-size", type=int, default=100, dest="min_indel_size", metavar="", help="Minimum indel length (bp) for diagnostic indels when using -include-indels (default: 100)")
    optional.add_argument("-t", "-threads", dest="threads", type=int, default=1, metavar="", help="Specify number of threads to use (n=1 by default)")
    optional.add_argument("-q", "-quiet", action="store_true", dest="quiet", help="Log to file only")
    optional.add_argument(
        "-extract-tracts",
        action="store_true",
        dest="extract_tracts",
        help="Split recombinant genomes by ancestry",
    )
    optional.add_argument(
        "-phylogeny",
        action="store_true",
        dest="phylogeny",
        help="Phylogeny of recombinant ancestors",
    )
    # If called with no arguments, show help (same output as --help) instead of erroring.
    if len(sys.argv) == 1:
        parser.print_help()
        return

    args = parser.parse_args()

    if getattr(args, "phylogeny", False):
        args.extract_tracts = True

    if getattr(args, "minor_ref_pct", None) is None:
        args.minor_ref_pct = MINOR_REF_PCT_THRESHOLD
    if args.minor_ref_pct < 0 or args.minor_ref_pct > 100:
        parser.error("-minor-ref-pct must be between 0 and 100")
    # breakpoint_min_snps is a fixed 1 (default) or 2 (when -b/-breakpoint-snp is used)

    if args.input is None:
        parser.error("-i/-input is required")

    # Resolve ref1/ref2: from -ref (defaults) or from -ref1/-ref2 (required if no -ref)
    # Normalize -ref labels: any case (Ia, ia, IA, IIb, iib, etc.) -> canonical Ia, Ib, IIa, IIb
    def _normalize_ref_label(s: str) -> str:
        t = s.strip().lower()
        if t == "ia": return "Ia"
        if t == "ib": return "Ib"
        if t == "iia": return "IIa"
        if t == "iib": return "IIb"
        return s.strip()

    ref1_spec_resolved: Optional[str] = None
    ref2_spec_resolved: Optional[str] = None
    ref1_g_resolved: Optional[str] = None
    ref2_g_resolved: Optional[str] = None

    if getattr(args, "ref", None):
        parts = [p.strip() for p in args.ref.split(",") if p.strip()]
        if len(parts) != 2:
            parser.error("-ref must be two comma-separated labels, e.g. Ia,Ib or IIa,IIb (got %r)" % getattr(args, "ref"))
        L1, L2 = _normalize_ref_label(parts[0]), _normalize_ref_label(parts[1])
        if L1 not in REF_DEFAULTS or L2 not in REF_DEFAULTS:
            parser.error("-ref labels must be among Ia, Ib, IIa, IIb (got %s, %s)" % (L1, L2))
        ref1_spec_resolved = REF_DEFAULTS[L1]
        ref2_spec_resolved = REF_DEFAULTS[L2]
        ref1_g_resolved = L1
        ref2_g_resolved = L2
    if args.ref1 is not None:
        ref1_spec_resolved = args.ref1
    if args.ref2 is not None:
        ref2_spec_resolved = args.ref2
    if getattr(args, "ref1_g", None) is not None:
        ref1_g_resolved = _normalize_ref_label(args.ref1_g) if args.ref1_g.strip().lower() in ("ia", "ib", "iia", "iib") else args.ref1_g.strip()
    if getattr(args, "ref2_g", None) is not None:
        ref2_g_resolved = _normalize_ref_label(args.ref2_g) if args.ref2_g.strip().lower() in ("ia", "ib", "iia", "iib") else args.ref2_g.strip()

    if ref1_spec_resolved is None or ref2_spec_resolved is None:
        parser.error("Either -ref LABEL1,LABEL2 (e.g. Ia,Ib) or both -ref1 and -ref2 are required")
    args.ref1 = ref1_spec_resolved
    args.ref2 = ref2_spec_resolved
    args.ref1_g = ref1_g_resolved if ref1_g_resolved is not None else _short_ref_label(args.ref1)
    args.ref2_g = ref2_g_resolved if ref2_g_resolved is not None else _short_ref_label(args.ref2)

    ref1_label = (getattr(args, "ref1_g", None) or _short_ref_label(args.ref1)).replace(".", "_")[:24]
    ref2_label = (getattr(args, "ref2_g", None) or _short_ref_label(args.ref2)).replace(".", "_")[:24]
    ref1_label = re.sub(r"[^\w]", "_", ref1_label) or "ref1"
    ref2_label = re.sub(r"[^\w]", "_", ref2_label) or "ref2"

    # Infer Squirrel clade from ref genotype labels: both Ia/Ib -> cladei; both IIa/IIb -> cladeii; mix -> None
    squirrel_clade = _infer_squirrel_clade(getattr(args, "ref1_g", None), getattr(args, "ref2_g", None))
    if squirrel_clade is None:
        squirrel_clade = _infer_squirrel_clade(ref1_label, ref2_label)
    is_intra_clade = squirrel_clade is not None
    minor_threshold = float(getattr(args, "minor_ref_pct", MINOR_REF_PCT_THRESHOLD))
    if squirrel_clade == "cladei":
        logger.info("Inferred Squirrel --clade cladei from ref1_g/ref2_g (Clade I)")
    elif squirrel_clade == "cladeii":
        logger.info("Inferred Squirrel Clade II from ref1_g/ref2_g")
    else:
        logger.info("Mixed or unspecified ref labels; Squirrel will run without -clade")

    args.output = Path(args.output_dir).resolve()
    if args.output.exists():
        logger.info("Removing existing output %s for new run", args.output)
        if args.output.is_dir():
            shutil.rmtree(args.output)
        else:
            args.output.unlink()
    args.output.mkdir(parents=True, exist_ok=True)
    logger.info("Output directory: %s", args.output)
    work_dir = args.output / "work"
    work_dir.mkdir(parents=True, exist_ok=True)
    setup_logging(args.output, verbose=not args.quiet)

    # When -ref L1,L2: build consensus refs from Pathoplexus (earliest 5 per clade → one consensus per clade) and use those FASTAs as ref1/ref2
    if getattr(args, "ref", None):
        L1 = ref1_g_resolved
        L2 = ref2_g_resolved
        consensus_dir = work_dir / "ref_consensus"
        try:
            consensus_map = _build_lapis_consensus_refs([L1, L2], consensus_dir, work_dir)
        except Exception as e:
            logger.error("Consensus refs: failed to build consensus references from Pathoplexus: %s", e)
            sys.exit(1)
        if L1 not in consensus_map or L2 not in consensus_map:
            logger.error("Consensus refs: missing consensus FASTA for %s or %s", L1, L2)
            sys.exit(1)
        args.ref1 = str(consensus_map[L1])
        args.ref2 = str(consensus_map[L2])

    ref_ia_path = resolve_ref(args.ref1, work_dir, "1")
    ref_ib_path = resolve_ref(args.ref2, work_dir, "2")
    if ref_ia_path is None:
        logger.error("Could not resolve ref1: %s", args.ref1)
        sys.exit(1)
    if ref_ib_path is None:
        logger.error("Could not resolve ref2: %s", args.ref2)
        sys.exit(1)

    ref_ia_id, ref_ia_unaligned = load_ref_sequence(ref_ia_path)
    ref_ib_id, ref_ib_unaligned = load_ref_sequence(ref_ib_path)
    ref_ia_key = _safe_fasta_id(ref_ia_id.replace("/", "_"))
    ref_ib_key = _safe_fasta_id(ref_ib_id.replace("/", "_"))

    # --- Step 1: Align ref Ia + ref Ib ONLY to find diagnostic SNPs and indels ---
    squirrel_out_refs = work_dir / "squirrel_out_refs"
    squirrel_out_refs.mkdir(parents=True, exist_ok=True)
    squirrel_in_refs = work_dir / "squirrel_input_refs_only.fa"
    with open(squirrel_in_refs, "w") as out:
        out.write(f">{ref_ia_key}\n")
        for i in range(0, len(ref_ia_unaligned), 80):
            out.write(ref_ia_unaligned[i : i + 80] + "\n")
        out.write(f">{ref_ib_key}\n")
        for i in range(0, len(ref_ib_unaligned), 80):
            out.write(ref_ib_unaligned[i : i + 80] + "\n")
    logger.info("Step 1: Built %s (ref Ia + ref Ib only) to find diagnostic sites", squirrel_in_refs)

    aln_refs_stem = squirrel_in_refs.stem + ".aln.fasta"
    squirrel_aln_refs = squirrel_out_refs / aln_refs_stem
    _run_squirrel(squirrel_clade, squirrel_in_refs, squirrel_out_refs, squirrel_aln_refs)
    logger.info("Using refs alignment: %s", squirrel_aln_refs)

    alignments_refs = load_alignment_fasta(squirrel_aln_refs)
    if not alignments_refs:
        logger.error("No sequences in refs alignment %s", squirrel_aln_refs)
        sys.exit(1)

    def find_ref_key(ref_key: str, keys: List[str]) -> Optional[str]:
        if ref_key in keys:
            return ref_key
        for k in keys:
            if k.startswith(ref_key + "_"):
                return k
        return None

    ref_ia_aln_key = find_ref_key(ref_ia_key, list(alignments_refs.keys()))
    ref_ib_aln_key = find_ref_key(ref_ib_key, list(alignments_refs.keys()))
    if ref_ia_aln_key is None:
        logger.error("Ref Ia (%s) not found in refs alignment; keys: %s", ref_ia_key, list(alignments_refs.keys())[:5])
        sys.exit(1)
    if ref_ib_aln_key is None:
        logger.error("Ref Ib (%s) not found in refs alignment; keys: %s", ref_ib_key, list(alignments_refs.keys())[:5])
        sys.exit(1)
    ref_ia_seq = alignments_refs[ref_ia_aln_key]
    ref_ib_seq = alignments_refs[ref_ib_aln_key]
    ref_len = len(ref_ia_seq)

    diagnostic_snps = build_diagnostic_snps_from_alignment(ref_ia_seq, ref_ib_seq)
    diagnostic_indels: Optional[List[Tuple[int, int, str]]] = None
    n_indel_columns: Optional[int] = None
    if getattr(args, "include_indels", False):
        min_indel = getattr(args, "min_indel_size", 100)
        diagnostic_indels = find_large_indels(ref_ia_seq, ref_ib_seq, min_size=min_indel)
        with open(work_dir / "diagnostic_indels.txt", "w") as f:
            f.write("start\tend\tref_with_bases\tdeletion_in\n")
            for start, end, ref_who in diagnostic_indels:
                deletion_in = "ib" if ref_who == "ia" else "ia"
                f.write(f"{start}\t{end}\t{ref_who}\t{deletion_in}\n")
        logger.info("Wrote work/diagnostic_indels.txt (%d large indels >= %d bp)", len(diagnostic_indels), min_indel)
        _write_indel_regions_side_by_side(work_dir / "ref1_ref2_indel_regions.txt", ref_ia_seq, ref_ib_seq, diagnostic_indels)
        n_indel_columns = sum(end - start + 1 for (start, end, _) in diagnostic_indels)
    else:
        logger.info("Using diagnostic SNPs only (default); use -include-indels to add large indels")
    if not diagnostic_snps and not (diagnostic_indels and len(diagnostic_indels) > 0):
        logger.error("No diagnostic sites found; ref Ia and ref Ib may be identical in alignment")
        sys.exit(1)
    n_diagnostic_sites = len(diagnostic_snps) + (n_indel_columns or 0)
    logger.info("Diagnostic sites: %d SNPs%s; %d total sites", len(diagnostic_snps), (f"; %d indel columns" % n_indel_columns) if n_indel_columns is not None else "", n_diagnostic_sites)

    # --- Step 2: Align ref Ia + queries, then classify each query at diagnostic positions ---
    try:
        query_input = resolve_query_input(args.input, work_dir)
    except FileNotFoundError as e:
        logger.error("%s", e)
        sys.exit(1)
    if query_input.is_dir():
        queries_fa = work_dir / "concatenated_input.fa"
        try:
            concatenate_fasta_dir(query_input, queries_fa)
        except FileNotFoundError as e:
            logger.error("%s", e)
            sys.exit(1)
    else:
        queries_fa = query_input
    squirrel_out_queries = work_dir / "squirrel_out_queries"
    squirrel_out_queries.mkdir(parents=True, exist_ok=True)
    squirrel_in_queries = work_dir / "squirrel_input_ref_and_queries.fa"
    queries_fa_sanitized = work_dir / "queries_sanitized.fa"
    _sanitize_fasta_ids(Path(queries_fa), queries_fa_sanitized)
    with open(squirrel_in_queries, "w") as out:
        # Use aligned refs from step 1 (ref_ia_seq, ref_ib_seq) so step 2 alignment has same coordinates
        # IMPORTANT: Squirrel is strict about FASTA IDs (no special characters like ':').
        # The step-1 alignment keys can include the full NCBI description (spaces->underscores),
        # which may contain ':' and break Squirrel. Use safe, short IDs here.
        out.write(f">{ref_ia_key}\n")
        for i in range(0, len(ref_ia_seq), 80):
            out.write(ref_ia_seq[i : i + 80] + "\n")
        out.write(f">{ref_ib_key}\n")
        for i in range(0, len(ref_ib_seq), 80):
            out.write(ref_ib_seq[i : i + 80] + "\n")
        with open(queries_fa_sanitized) as f:
            out.write(f.read())
    logger.info("Step 2: Built %s (ref Ia + ref Ib + consensus genomes) to classify at diagnostic positions", squirrel_in_queries)

    aln_queries_stem = squirrel_in_queries.stem + ".aln.fasta"
    squirrel_aln_queries = squirrel_out_queries / aln_queries_stem
    _run_squirrel(squirrel_clade, squirrel_in_queries, squirrel_out_queries, squirrel_aln_queries)
    logger.info("Using queries alignment: %s", squirrel_aln_queries)

    alignments_queries = load_alignment_fasta(squirrel_aln_queries)
    if not alignments_queries:
        logger.error("No sequences in queries alignment %s", squirrel_aln_queries)
        sys.exit(1)
    ref_ia_in_queries = find_ref_key(ref_ia_key, list(alignments_queries.keys()))
    ref_ib_in_queries = find_ref_key(ref_ib_key, list(alignments_queries.keys()))
    if ref_ia_in_queries is None:
        logger.error("Ref Ia not found in queries alignment; keys: %s", list(alignments_queries.keys())[:5])
        sys.exit(1)
    if ref_ib_in_queries is None:
        logger.error("Ref Ib not found in queries alignment; keys: %s", list(alignments_queries.keys())[:5])
        sys.exit(1)
    aln_ref_len = len(alignments_queries[ref_ia_in_queries])
    if aln_ref_len != ref_len:
        logger.warning("Queries alignment length %d != refs alignment length %d; using refs length for diagnostic positions", aln_ref_len, ref_len)

    diagnostic_snp_positions = [p for (p, _, _) in diagnostic_snps]
    ref_keys_in_queries = {ref_ia_in_queries, ref_ib_in_queries}
    results = []
    query_items = [(qid, alignments_queries[qid]) for qid in alignments_queries if qid not in ref_keys_in_queries]
    for i, (query_id, query_seq) in enumerate(query_items):
        logger.info("Processing %s (%d/%d)", query_id, i + 1, len(query_items))
        allegiances = get_query_allegiance_from_alignment(
            query_seq, diagnostic_snps, ref_len, diagnostic_indels=diagnostic_indels
        )
        if not allegiances:
            logger.warning("Query %s: no diagnostic calls", query_id)
            continue
        n_ia, n_ib, n_other, n_other_n = allegiance_summary(allegiances)
        total = n_ia + n_ib + n_other
        pct_ia = round(100.0 * n_ia / total, 2) if total else 0
        pct_ib = round(100.0 * n_ib / total, 2) if total else 0
        pct_other = round(100.0 * n_other / total, 2) if total else 0
        # SNP-only summary for consensus and deletion present (SNP-based interpretation)
        n_ia_snp, n_ib_snp, n_other_snp, _ = allegiance_summary_snp_only(allegiances, diagnostic_snp_positions)
        total_snp = n_ia_snp + n_ib_snp + n_other_snp
        consensus_snp = consensus_from_snp_percentages(
            n_ia_snp, n_ib_snp, n_other_snp, ref1_label, ref2_label, pct_threshold=10.0
        )
        deletion_present = (consensus_snp == ref2_label) if consensus_snp != "other" else None
        minor_ref_pct = min(pct_ia, pct_ib)
        recombinant_call = _recombinant_call_minor_pct(n_ia, n_ib, total, minor_threshold)
        rec = {
            "id": query_id,
            "length": len(query_seq),
            "n_diagnostic_snps": n_diagnostic_sites,
            "n_ia": n_ia,
            "n_ib": n_ib,
            "n_other": n_other,
            "n_other_n": n_other_n,
            "pct_ia": pct_ia,
            "pct_ib": pct_ib,
            "pct_other": pct_other,
            "consensus_snp": consensus_snp,
            "deletion_present": deletion_present,
            "minor_ref_pct": minor_ref_pct,
            "recombinant_call": recombinant_call,
            "allegiances": allegiances,
        }
        results.append(rec)

    # TSV header: match main table semantics.
    has_ancestors = any("ancestors" in r for r in results)
    header_parts = [
        "id", "length", "n_sites", f"n_{ref1_label}", f"n_{ref2_label}", "n_other",
        f"pct_{ref1_label}", f"pct_{ref2_label}", "pct_other",
    ]
    if has_ancestors:
        header_parts.append("ancestors")
        header_parts.append("recombinant_call")
    else:
        header_parts.append("consensus_SNP")
        header_parts.append("recombinant_call")
    header = "\t".join(header_parts) + "\n"

    def row(r: Dict[str, Any]) -> str:
        parts = [
            r["id"], str(r["length"]), str(r["n_diagnostic_snps"]), str(r["n_ia"]), str(r["n_ib"]), str(r["n_other"]),
            str(r["pct_ia"]), str(r["pct_ib"]), str(r["pct_other"]),
        ]
        if has_ancestors:
            anc = r.get("ancestors", {})
            if isinstance(anc, dict):
                v1 = anc.get("ref1", "")
                v2 = anc.get("ref2", "")
                parts.append(f"{v1}|{v2}")
            else:
                parts.append(str(anc) if anc is not None else "")
            parts.append(r["recombinant_call"])
        else:
            parts.append(r["consensus_snp"])
            parts.append(r["recombinant_call"])
        return "\t".join(str(p) for p in parts) + "\n"

    out_tsv = args.output / "recmpox_results.tsv"
    with open(out_tsv, "w") as f:
        f.write(header)
        for r in results:
            f.write(row(r))

    # Diagnostic sites per potential recombinant: sample_id, diagnostic_site, clade_classification (Ia/Ib/other)
    out_sites_tsv = args.output / "potential_recombinants_diagnostic_sites.tsv"
    rec_samples = [r for r in results if r.get("recombinant_call") == "potential recombinant"]
    if rec_samples:
        with open(out_sites_tsv, "w") as f:
            f.write("sample_id\tdiagnostic_site\tclade_classification\n")
            for r in rec_samples:
                sample_id = r["id"]
                for (pos, allegiance) in r.get("allegiances", []):
                    if allegiance == "ia":
                        label = ref1_label
                    elif allegiance == "ib":
                        label = ref2_label
                    else:
                        label = "other"
                    f.write(f"{sample_id}\t{pos}\t{label}\n")
        logger.info("Wrote %s (%d potential recombinant samples)", out_sites_tsv, len(rec_samples))

    other_explanation = (
        "%% other = diagnostic sites where the query neither matched %s nor %s (different base or gap/N count as other)."
    ) % (ref1_label, ref2_label)

    # Optional: extract per-clade tract sequences for potential recombinants (Ia-only / Ib-only positions)
    if getattr(args, "extract_tracts", False):
        if not getattr(args, "include_indels", False):
            logger.info(
                "--extract-tracts tip: consider adding -include-indels so that large "
                "clade-specific deletions are included in tract boundaries."
            )
        _extract_tract_sequences(
            out_dir=args.output / "tracts",
            results=results,
            alignments_queries=alignments_queries,
            ref1_label=ref1_label,
            ref2_label=ref2_label,
            is_intra_clade=is_intra_clade,
            min_consecutive=int(getattr(args, "breakpoint_min_snps", 1)),
            include_indels=getattr(args, "include_indels", False),
        )

    # Optional: phylogeny from refs + two partition FASTAs (requires extract-tracts)
    if getattr(args, "phylogeny", False):
        _run_phylogeny_pipeline(args, work_dir, ref1_label, ref2_label, squirrel_clade)
        # Attach phylogeny ancestors (if any) to per-genome records for HTML/TSV
        # and refine recombinant_call based on ancestor mismatch.
        phy_anc = getattr(args, "_phylogeny_ancestors", None)
        if phy_anc:
            for rec in results:
                sid = rec.get("id")
                if sid in phy_anc:
                    rec["ancestors"] = phy_anc[sid]
            # Refine calls using ancestors:
            # - If both ancestors are present and identical → no recombinant
            # - If they are different and the sample passed the % threshold → recombinant
            for rec in results:
                anc = rec.get("ancestors")
                if not isinstance(anc, dict):
                    continue
                a1 = anc.get("ref1")
                a2 = anc.get("ref2")
                if not (a1 and a2):
                    continue
                # Same ancestor on both tracts → treat as non-recombinant regardless of % threshold.
                if a1 == a2:
                    rec["recombinant_call"] = "no recombinant"
                # Different ancestors and previously flagged by % threshold → recombinant.
                elif rec.get("recombinant_call") == "potential recombinant":
                    rec["recombinant_call"] = "recombinant"

    # Summary note: when --phylogeny was used and ancestors are present, explain the two-step logic.
    has_ancestors = any("ancestors" in r for r in results)
    if has_ancestors:
        recombinant_threshold_note = (
            f"Recombinant calls use two steps. (1) A {minor_threshold:g}% minor-ref threshold screens for potential recombinants. "
            f"(2) With --phylogeny, the inferred ancestors of the two tracts are compared: if they are the same outbreak (overlap), the sample is classified as no recombinant; if they differ, it is classified as recombinant."
        )
    else:
        recombinant_threshold_note = (
            f"A {minor_threshold:g}% threshold is used for all recombinant calls: "
            f"when minor ref % ≥ {minor_threshold:g}%, the sample is flagged as potential recombinant."
        )

    # HTML: one file if <= HTML_CHUNK_SIZE genomes, else one file per chunk (written after phylogeny so PDF can be included)
    html_files: List[Path] = []
    n_snps = len(diagnostic_snps)
    phylogeny_pdf = (args.output / "phylogeny" / "phylogeny_tree.pdf") if getattr(args, "phylogeny", False) else None
    phylogeny_svg = (args.output / "phylogeny" / "phylogeny_tree.svg") if getattr(args, "phylogeny", False) else None
    phylogeny_pdf_b64: Optional[str] = None
    phylogeny_svg_content: Optional[str] = None
    if phylogeny_pdf and phylogeny_pdf.exists():
        try:
            phylogeny_pdf_b64 = base64.b64encode(phylogeny_pdf.read_bytes()).decode("ascii")
        except OSError:
            pass
    if phylogeny_svg and phylogeny_svg.exists():
        try:
            phylogeny_svg_content = phylogeny_svg.read_text(encoding="utf-8")
        except OSError:
            pass
    if len(results) <= HTML_CHUNK_SIZE:
        out_html = args.output / "recmpox_results.html"
        _write_results_html(out_html, results, ref1_label, ref2_label, recombinant_threshold_note, other_explanation, is_intra_clade, minor_threshold, breakpoint_min_consecutive_snps=int(getattr(args, "breakpoint_min_snps", 1)), n_diagnostic_snps=n_snps, n_indel_columns=n_indel_columns, ref1_spec=args.ref1, ref2_spec=args.ref2, diagnostic_snp_positions=[p for (p, _, _) in diagnostic_snps], genome_length=ref_len, phylogeny_pdf_path=phylogeny_pdf, phylogeny_pdf_base64=phylogeny_pdf_b64, phylogeny_svg_content=phylogeny_svg_content)
        html_files.append(out_html)
        logger.info("Wrote %s", out_html)
    else:
        chunks = [results[i:i + HTML_CHUNK_SIZE] for i in range(0, len(results), HTML_CHUNK_SIZE)]
        for part, chunk in enumerate(chunks, start=1):
            out_html = args.output / f"recmpox_results_{part}.html"
            _write_results_html(out_html, chunk, ref1_label, ref2_label, recombinant_threshold_note, other_explanation, is_intra_clade, minor_threshold, breakpoint_min_consecutive_snps=int(getattr(args, "breakpoint_min_snps", 1)), part_index=part, total_parts=len(chunks), n_diagnostic_snps=n_snps, n_indel_columns=n_indel_columns, ref1_spec=args.ref1, ref2_spec=args.ref2, diagnostic_snp_positions=[p for (p, _, _) in diagnostic_snps], genome_length=ref_len, phylogeny_pdf_path=phylogeny_pdf, phylogeny_pdf_base64=phylogeny_pdf_b64, phylogeny_svg_content=phylogeny_svg_content)
            html_files.append(out_html)
            logger.info("Wrote %s (%d genomes)", out_html, len(chunk))
        # Index page linking to all parts
        index_path = args.output / "recmpox_results_index.html"
        with open(index_path, "w") as f:
            f.write("<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"UTF-8\"><title>RecMpox Results – Index</title></head><body><h1>RecMpox Results</h1><p>%d genomes in %d parts (max %d per file).</p><ul>\n" % (len(results), len(chunks), HTML_CHUNK_SIZE))
            for part in range(1, len(chunks) + 1):
                f.write('<li><a href="recmpox_results_%d.html#diagnosticStripsSection">Part %d of %d</a> (table + bar chart + diagnostic strips)</li>\n' % (part, part, len(chunks)))
            f.write("</ul></body></html>")
        logger.info("Wrote %s", index_path)
        html_files.insert(0, index_path)

    # Zip HTML report for easier sharing (e.g. email)
    zip_path = args.output / "recmpox_results.zip"
    try:
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
            for p in html_files:
                if p.exists():
                    zf.write(p, arcname=p.name)
        logger.info("Created %s (HTML report archive)", zip_path)
    except OSError as e:
        logger.warning("Could not create HTML zip %s: %s", zip_path, e)

    with open(args.output / "diagnostic_snps.txt", "w") as f:
        f.write("position\tia_allele\tib_allele\n")
        for pos, ia_a, ib_a in diagnostic_snps:
            f.write(f"{pos}\t{ia_a}\t{ib_a}\n")

    # Remove work dir (intermediate files not needed)
    if work_dir.exists():
        shutil.rmtree(work_dir)
        logger.info("Removed intermediate files (work/)")

    if n_indel_columns is not None:
        logger.info(
            "Wrote %s (%d genomes). Diagnostic SNPs: %d; indel columns: %d; total diagnostic sites: %d.",
            out_tsv, len(results), n_snps, n_indel_columns, n_diagnostic_sites,
        )
    else:
        logger.info(
            "Wrote %s (%d genomes). Diagnostic SNPs: %d; total diagnostic sites: %d.",
            out_tsv, len(results), n_snps, n_diagnostic_sites,
        )
    html_str = html_files[0].name if len(html_files) == 1 else "index: " + html_files[0].name + " + " + ", ".join(p.name for p in html_files[1:])
    print(f"Done. {len(results)} genomes. Results: {out_tsv}  HTML: {html_str}")
    if n_indel_columns is not None:
        print(f"  Diagnostic SNPs: {n_snps}; indel columns: {n_indel_columns}; total diagnostic sites: {n_diagnostic_sites}")
    else:
        print(f"  Diagnostic SNPs: {n_snps}; total diagnostic sites: {n_diagnostic_sites}")


if __name__ == "__main__":
    main()
