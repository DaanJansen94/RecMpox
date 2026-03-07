"""
Diagnostic-SNP methodology for detecting recombinant mpox (Clade Ia vs Ib).

1. Build diagnostic SNPs: positions where ref_ia and ref_ib consistently differ.
2. For each query: allegiance (ia/ib/ambiguous) at each diagnostic SNP.
3. Rolling window of N SNPs: proportion Ia; flag when mixed (e.g. 30-70%).
4. Output: report + plot (position vs ia_fraction).
"""

import logging
import re
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# Pileup base codes: . , = A C G T N etc; ^] for start; $ for end; [+-]N[ACGT]+ for indel
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
        "minimap2", "-a", "-x", "asm5", "-t", str(threads),
        str(ref_fa), str(query_fa),
    ]
    try:
        with open(out_sam, "w") as f:
            subprocess.run(cmd_align, check=True, stdout=f, stderr=subprocess.PIPE, text=True, timeout=600)
        subprocess.run(
            ["samtools", "view", "-b", "-o", str(out_bam.with_suffix(".bam.tmp")), str(out_sam)],
            check=True, capture_output=True, text=True, timeout=60,
        )
        subprocess.run(
            ["samtools", "sort", "-o", str(out_bam), str(out_bam.with_suffix(".bam.tmp"))],
            check=True, capture_output=True, text=True, timeout=120,
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
    name, seq = None, []
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
        chrom, pos_s, ref_base, _, bases = parts[0], parts[1], parts[2], parts[3], parts[4]
        pos = int(pos_s)
        decoded = []
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
    # Fill positions with no coverage as gap (deletion in aligned seq relative to ref)
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
    diagnostic = []
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
        n, both_base, same_base, len(diagnostic),
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
    # Gap (-) or N, ., * (missing/poor) → count as the ref that has the deletion
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
    ref_id, ref_ia_sequence = ref_ia_seq
    ref_len = len(ref_ia_sequence)
    ref_ib_bases = get_bases_from_mpileup(ref_ib_bam, ref_ia_fa, ref_length=ref_len, include_gaps=True)
    diagnostic = []
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
        # Only SNPs: both must have a base and differ (exclude Ib deletion/insertion)
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
    current_id = None
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


def get_query_allegiance_from_alignment(
    query_seq: str,
    diagnostic_snps: List[Tuple[int, str, str]],
    ref_length: int,
    diagnostic_indels: Optional[List[Tuple[int, int, str]]] = None,
) -> List[Tuple[int, str]]:
    """
    At each diagnostic (pos, ia_allele, ib_allele), get query base from alignment column (1-based).
    N, gap, or other ambiguous -> 'ambiguous'. Same allegiance rules as get_query_allegiance.
    If diagnostic_indels is provided, append one allegiance per column in each indel (each position counts as a site).
    """
    if len(query_seq) < ref_length:
        logger.warning("Alignment sequence length %d < ref length %d", len(query_seq), ref_length)
    result = []
    for (pos, ia_a, ib_a) in diagnostic_snps:
        idx = pos - 1
        if idx >= len(query_seq):
            result.append((pos, "ambiguous"))
            continue
        q = query_seq[idx].upper()
        if q not in "ACGT-":
            result.append((pos, "ambiguous"))
            continue
        if q == "N":
            result.append((pos, "ambiguous"))
            continue
        # Same allegiance logic as get_query_allegiance
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
    result = []
    for (pos, ia_a, ib_a) in diagnostic_snps:
        q = query_bases.get(pos, "-").strip()
        if q not in "ACGT-":
            q = "-"
        if q == "*":
            q = "-"
        # Both alleles are base (SNP)
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
        # Ia has base, Ib has gap (deletion in Ib)
        elif ia_a in "ACGT" and ib_a == "-":
            if q in "ACGT":
                result.append((pos, "ia"))
            elif q == "-":
                result.append((pos, "ib"))
            else:
                result.append((pos, "ambiguous"))
        # Ia has gap, Ib has base (insertion in Ib)
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


def allegiance_summary(positions_allegiances: List[Tuple[int, str]]) -> Tuple[int, int, int]:
    """Return (n_ia, n_ib, n_ambiguous) over all diagnostic SNPs."""
    n_ia = sum(1 for _, a in positions_allegiances if a == "ia")
    n_ib = sum(1 for _, a in positions_allegiances if a == "ib")
    n_amb = sum(1 for _, a in positions_allegiances if a == "ambiguous")
    return (n_ia, n_ib, n_amb)


def allegiance_summary_snp_only(
    positions_allegiances: List[Tuple[int, str]],
    diagnostic_snp_positions: List[int],
) -> Tuple[int, int, int]:
    """
    Return (n_ia, n_ib, n_ambiguous) over diagnostic SNP positions only (excludes indel columns).
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
    endpoints and SNP count.  This makes tracts robust to patchy coverage: a stretch of
    N's inside an otherwise consistent Ia region will not split the tract.

    When ignore_other=False: any "other" position ends the current run (original
    behaviour — more conservative, breaks tracts at every ambiguous site).

    A breakpoint is only called when *both* flanking runs have >= min_consecutive SNPs.

    Returns:
        runs: list of (start_pos, end_pos, clade, n_snps) with clade in ("ia", "ib").
        breakpoints: list of (end_pos_before, start_pos_after, clade_before, clade_after).
    """
    snp_positions = set(diagnostic_snp_positions)
    # Keep only SNP positions, sorted by position
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
        last_clade_pos = pos  # last position that actually matched current_clade
        n_snps = 1
        i += 1
        while i < len(ordered):
            next_pos, next_a = ordered[i]
            if next_a == current_clade:
                last_clade_pos = next_pos
                n_snps += 1
                i += 1
            elif next_a in ("other", "ambiguous") and ignore_other:
                # Transparent: skip without updating last_clade_pos or n_snps
                i += 1
            else:
                # Opposing clade (or "other" when ignore_other=False) → end this run
                break
        runs.append((start_pos, last_clade_pos, current_clade, n_snps))

    breakpoints: List[Tuple[int, int, str, str]] = []
    for j in range(len(runs) - 1):
        start_a, end_a, clade_a, n_a = runs[j]
        start_b, end_b, clade_b, n_b = runs[j + 1]
        if clade_a == clade_b:
            continue
        # Breakpoint only if *both* runs have >= min_consecutive SNPs (ignore single-SNP runs)
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
    out = []
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
    regions = []
    in_region = False
    start = None
    for pos, frac in positions_ia_frac:
        if min_mixed <= frac <= max_mixed:
            if not in_region:
                start = pos
                in_region = True
        else:
            if in_region:
                regions.append((start, pos))
                in_region = False
    if in_region:
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
