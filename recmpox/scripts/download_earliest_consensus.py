#!/usr/bin/env python3
"""
Fetch earliest 5 genomes per selected clade; build one majority-rule consensus
FASTA per clade. Pathoplexus LAPIS (open + restricted), length >= 190 kb.
Requires Squirrel and/or mafft.

Usage:
  python download_earliest_consensus.py [--clades Ia,Ib,IIa,IIb] [--out-dir DIR]
  --clades: which clades (default: all four). Ia, Ib, IIa, IIb (comma or space separated)
  --out-dir: if set, write consensus FASTA files only here (e.g. for RecMpox); otherwise script directory
"""

import argparse
import json
import shutil
import ssl
import subprocess
import sys
import tempfile
import urllib.parse
import urllib.request
from pathlib import Path

LAPIS_MPOX_DETAILS = "https://lapis.pathoplexus.org/mpox/sample/details"
PATHOPLEXUS_FASTA = "https://pathoplexus.org/seq"
NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
MIN_LENGTH_BP = 190_000
PER_GROUP = 5
IA_SH2024_MIN_DATE = "2024-08-19"


def _fetch_url(url: str, params: str = "") -> str:
    full = f"{url}?{params}" if params else url
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    req = urllib.request.Request(full, headers={"User-Agent": "RecMpox/1.0"})
    with urllib.request.urlopen(req, timeout=90, context=ctx) as resp:
        return resp.read().decode("utf-8")


def _lapis_fetch(q: dict, limit: int = 5000) -> list:
    """Fetch from LAPIS details; return list of dicts (raw rows)."""
    q["limit"] = limit
    params = urllib.parse.urlencode(q)
    try:
        data = _fetch_url(LAPIS_MPOX_DETAILS, params)
        obj = json.loads(data)
        return obj.get("data") or []
    except Exception as e:
        print(f"Pathoplexus LAPIS failed: {e}", file=sys.stderr)
        return []


def _row_to_tuple(r: dict, date_key_priority=None) -> tuple | None:
    """(accession_version, date_sort_key, length, insdc) or None. date_key_priority: list of keys for date."""
    acc_ver = (r.get("accessionVersion") or r.get("accession") or "").strip()
    length = r.get("length")
    if not acc_ver or length is None or length < MIN_LENGTH_BP:
        return None
    insdc = (r.get("insdcAccessionFull") or "").strip() or None
    date = None
    for key in date_key_priority or ["sampleCollectionDate", "sampleCollectionDateRangeLower", "sampleCollectionDateRangeUpper"]:
        v = r.get(key)
        if v and isinstance(v, str) and v.strip():
            date = v.strip()
            break
    if not date:
        date = "9999-99-99"
    return (acc_ver, date, int(length), insdc)


def fetch_ib_kinshasa(limit: int = 5000) -> list:
    """Earliest 5 Ib Kinshasa (length >= 190k)."""
    rows = _lapis_fetch({
        "geoLocCountry": "Democratic Republic of the Congo",
        "geoLocAdmin1": "Kinshasa",
        "clade": "Ib",
        "lengthFrom": MIN_LENGTH_BP,
    }, limit=limit)
    out = []
    for r in rows:
        t = _row_to_tuple(r, ["sampleCollectionDate"])
        if t and t[1] != "9999-99-99":
            out.append(t)
    out.sort(key=lambda x: (x[1], x[0]))
    return out[:PER_GROUP]


def fetch_ia_kinshasa_sh2024(limit: int = 5000) -> list:
    """Earliest 5 Ia Kinshasa, outbreak sh2024, collection date >= 2024-08-19."""
    rows = _lapis_fetch({
        "geoLocCountry": "Democratic Republic of the Congo",
        "geoLocAdmin1": "Kinshasa",
        "clade": "Ia",
        "outbreak": "sh2024",
        "lengthFrom": MIN_LENGTH_BP,
    }, limit=limit)
    out = []
    for r in rows:
        t = _row_to_tuple(r, ["sampleCollectionDate"])
        if t and t[1] != "9999-99-99" and t[1] >= IA_SH2024_MIN_DATE:
            out.append(t)
    out.sort(key=lambda x: (x[1], x[0]))
    return out[:PER_GROUP]


def fetch_sh2017(limit: int = 5000) -> list:
    """Earliest 5 outbreak sh2017, length >= 190k."""
    rows = _lapis_fetch({
        "outbreak": "sh2017",
        "lengthFrom": MIN_LENGTH_BP,
    }, limit=limit)
    out = []
    for r in rows:
        t = _row_to_tuple(r, ["sampleCollectionDate"])
        if t and t[1] != "9999-99-99":
            out.append(t)
    out.sort(key=lambda x: (x[1], x[0]))
    return out[:PER_GROUP]


def fetch_iia_earliest(limit: int = 5000) -> list:
    """Earliest 5 IIa by date (use range fields if sampleCollectionDate missing)."""
    rows = _lapis_fetch({
        "clade": "IIa",
        "lengthFrom": MIN_LENGTH_BP,
    }, limit=limit)
    seen_base = set()
    out = []
    for r in rows:
        if r.get("versionStatus") != "LATEST_VERSION":
            continue
        t = _row_to_tuple(r, ["sampleCollectionDate", "sampleCollectionDateRangeLower", "sampleCollectionDateRangeUpper"])
        if not t:
            continue
        base = (t[0].split(".")[0], t[1])
        if base in seen_base:
            continue
        seen_base.add(base)
        out.append(t)
    out.sort(key=lambda x: (x[1], x[0]))
    return out[:PER_GROUP]


def fetch_fasta_pathoplexus(accession_version: str, out_path: Path) -> bool:
    url = f"{PATHOPLEXUS_FASTA}/{accession_version}.fa"
    try:
        data = _fetch_url(url, "")
        if not data.strip() or "not found" in data.lower()[:200]:
            return False
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(data)
        return True
    except Exception:
        return False


def fetch_fasta_ncbi(accession: str, out_path: Path) -> bool:
    params = f"db=nucleotide&id={urllib.parse.quote(accession.strip())}&rettype=fasta&retmode=text"
    data = _fetch_url(NCBI_EFETCH, params)
    if not data.strip() or "Error" in data[:200]:
        return False
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(data)
    return True


def align_group(combined_fa: Path, out_dir: Path, stem: str, use_clade_ii: bool) -> Path | None:
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


def build_consensus_from_aln(aln_path: Path, consensus_stem: str) -> tuple[str, int] | None:
    """Return (consensus_content_with_header, len_ungapped) or None."""
    seqs = {}
    with open(aln_path) as f:
        current_id = None
        current_seq = []
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
    consensus = []
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
    n_non_n = sum(1 for b in consensus_ungapped if b.upper() in "ACGT")
    total = len(consensus_ungapped)
    body = "\n".join(consensus_ungapped[i : i + 80] for i in range(0, len(consensus_ungapped), 80)) + "\n"
    header = f">{consensus_stem}\n"
    return (header + body, len(consensus_ungapped))


def main() -> int:
    ap = argparse.ArgumentParser(description="Earliest 5 per selected clade → one consensus FASTA per clade.")
    ap.add_argument("--clades", nargs="+", default=["Ia", "Ib", "IIa", "IIb"], metavar="CLADE", help="Clades: Ia, Ib, IIa, IIb (comma or space separated)")
    ap.add_argument("--out-dir", type=Path, default=None, help="If set, write consensus FASTA files only here (e.g. for RecMpox); otherwise script directory")
    args = ap.parse_args()
    valid = {"Ia", "Ib", "IIa", "IIb"}
    clades = []
    for c in args.clades:
        clades.extend(x.strip() for x in c.split(",") if x.strip())
    args.clades = [c for c in clades if c in valid]
    if not args.clades and clades:
        print("Invalid --clades. Choose from: Ia, Ib, IIa, IIb", file=sys.stderr)
        return 1
    if not args.clades:
        args.clades = ["Ia", "Ib", "IIa", "IIb"]
    script_dir = Path(__file__).resolve().parent
    consensus_output_dir = args.out_dir.resolve() if (getattr(args, "out_dir", None) and args.out_dir) else script_dir
    out_dir = Path(tempfile.mkdtemp(prefix="earliest_consensus_"))
    try:
        fasta_dir = out_dir / "fasta"
        fasta_dir.mkdir(parents=True, exist_ok=True)

        all_groups = [
            ("Ib_Kinshasa", fetch_ib_kinshasa, "ib_kinshasa", False, "Ib", "sh2023Ib"),
            ("Ia_Kinshasa_sh2024", fetch_ia_kinshasa_sh2024, "ia_kinshasa", False, "Ia", "sh2024Ia"),
            ("sh2017", fetch_sh2017, "sh2017", True, "IIb", "sh2017IIb"),
            ("IIa", fetch_iia_earliest, "iia", True, "IIa", "iia"),
        ]
        groups = [g for g in all_groups if g[4] in args.clades]
        if not groups:
            print("No clades selected.", file=sys.stderr)
            return 1

        n_written = 0
        for name, fetch_fn, consensus_stem, use_clade_ii, _, header_stem in groups:
            print(f"Fetching {name}...")
            rows = fetch_fn()
            if len(rows) < 2:
                print(f"  Skip {name}: need at least 2 samples (got {len(rows)}).", file=sys.stderr)
                continue

            group_fastas = []
            for acc_ver, date, length, insdc in rows:
                safe_id = acc_ver.replace(".", "_").replace("/", "_")
                path = fasta_dir / f"{name}_{safe_id}.fa"
                if not fetch_fasta_pathoplexus(acc_ver, path) and insdc:
                    fetch_fasta_ncbi(insdc.split(".")[0], path)
                if path.exists():
                    group_fastas.append(path)
            if len(group_fastas) < 2:
                print(f"  Skip {name}: could not download enough FASTAs.", file=sys.stderr)
                continue

            combined_fa = out_dir / f"samples_combined_{consensus_stem}.fa"
            with open(combined_fa, "w") as out:
                for p in sorted(group_fastas):
                    text = p.read_text()
                    out.write(text)
                    if text and not text.endswith("\n"):
                        out.write("\n")

            aln_path = align_group(combined_fa, out_dir, consensus_stem, use_clade_ii)
            if not aln_path or not aln_path.exists():
                print(f"  Skip {name}: alignment failed (install squirrel and/or mafft).", file=sys.stderr)
                continue

            result = build_consensus_from_aln(aln_path, header_stem)
            if not result:
                print(f"  Skip {name}: consensus build failed.", file=sys.stderr)
                continue
            content, length_bp = result
            lines = content.split("\n")
            if lines and lines[0].startswith(">"):
                lines[0] = f">{header_stem}"
            content = "\n".join(lines) + ("\n" if content.endswith("\n") else "")
            out_fa = consensus_output_dir / f"{header_stem}.fa"
            out_fa.parent.mkdir(parents=True, exist_ok=True)
            out_fa.write_text(content)
            print(f"  Wrote {out_fa.name} (length {length_bp} bp)")
            n_written += 1

        if n_written == 0:
            print("No consensus files produced.", file=sys.stderr)
            return 1
        return 0
    finally:
        shutil.rmtree(out_dir, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
