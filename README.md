# RecMpox

Current release: **v0.0.3**

RecMpox is a command-line tool that **flags potential recombination events** in monkeypox viruses. It does not confirm recombination, but highlights genomes that may be recombinant and warrant further investigation. RecMpox works by detecting regions within a genome that appear to originate from two different parental viruses. Such patterns are not conclusive evidence of recombination, as similar signals can also arise from shared ancestral variation, convergent mutations, mixed populations (e.g., co-infections or laboratory contamination), or sequencing and assembly errors.


### How RecMpox Works?
1. **References are required**: RecMpox compares your genomes against two reference sequences (for example, Clade Ia vs. Ib, or Ib vs. IIb), because recombination can only occur between two distinct lineages.
2. **Alignment and diagnostic SNPs**: The two reference genomes are aligned using [Squirrel](https://github.com/aineniamh/squirrel), so that the same genomic positions correspond across all sequences. RecMpox then identifies positions where the two references differ at the same coordinates. These positions are defined as diagnostic SNPs, because they distinguish between the reference lineages. Positions where the references are identical are ignored, as they do not provide information for detecting recombination.
3. **Consensus genome classification**: Your consensus genomes are aligned to the same references. At each diagnostic SNP, the base is classified as matching reference 1, reference 2, or other (e.g., gaps or ambiguous bases).
4. **Flagging potential recombinants**: If both references contribute at least 10% of the diagnostic positions in a genome, RecMpox flags it as a potential recombinant, since no single lineage clearly dominates.
5. **Recombination tracts and breakpoints**: By examining the pattern of reference matches along the genome, RecMpox infers recombination tracts and identifies their breakpoints (start and end positions). By default, no consecutive-SNP filtering is applied (minimum run length = 1), but you can ignore single-SNP runs by adding `-breakpoint-snp` (or `-b`), which sets the minimum run length to 2.
6. **Outputs**:
   - TSV file: or each genome, reports the number and proportion of diagnostic SNPs matching each reference, the resulting recombinant flag, and summary statistics used for tract inference.
   - Interactive HTML report: Provides sortable tables, summary plots, per-sample visualisations, and genome-wide displays of inferred recombination tracts and breakpoints.
   - Aligned FASTA: Contains the aligned reference and query sequences used for analysis.

⚠️ **Note**: RecMpox is primarily designed to investigate potential recombination between viruses circulating in sustained human outbreaks (for example, SH2017, SH2023b, and SH2024a). The reference genomes provided by default in the tool correspond to these sustained outbreak lineages. When applying RecMpox outside this context, it is crucial to select reference genomes that are genetically close to your consensus sequences. Using distant or poorly matched references can reduce the interpretability of diagnostic SNPs and may lead to misleading recombinant signals.

🔴 **Caution — Intra-clade comparisons**: When comparing within a clade, set the threshold higher than the default. Intra-clade comparisons yield far fewer diagnostic SNPs (e.g. as low as ~120 between Ia and Ib), meaning each individual SNP carries more weight and small percentages can arise from convergent evolution rather than true recombination. We recommend a minimum threshold of **20%** (`-m 20`) for intra-clade recombinant screening.

## Installation

### Prerequisites
First, install conda if you haven't already:
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Then, ensure you have the required channels:
```bash
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels defaults
conda config --set channel_priority strict
```

### Option 1: Using Conda (Recommended)
Standard install [RecMpox via Conda](https://anaconda.org/bioconda/recmpox):
```bash
conda create -n recmpox -c conda-forge -c bioconda recmpox -y
conda activate recmpox
```

OR, if the above fails:
```bash
conda create -y -n recmpox -c conda-forge -c bioconda recmpox python=3.11 --solver libmamba --strict-channel-priority
conda activate recmpox
```

### Option 2: From Source Code
1. Create conda environment with required tools and install RecMpox
   ```bash
   git clone https://github.com/DaanJansen94/RecMpox.git
   cd RecMpox
   conda env create -f environment-recmpox.yml
   conda activate recmpox
   pip install .
   ```

2. Re-installation (when updates are available):
   ```bash
   conda activate RecMpox  # Make sure you're in the right environment
   cd RecMpox
   git pull  # Get the latest updates from GitHub
   pip uninstall RecMpox
   pip install .
   ```
   
## Usage

### Basic usage

```bash
# Use built-in references: consensus of 5 earliest genomes obtained from sustained outbreaks
recmpox -i fasta/ -o output -ref Ia,Ib -t 4
recmpox -i fasta/ -o output -ref Ib,IIb -t 4

# Input can be: FASTA file, directory of .fa/.fasta/.fna, or NCBI accession(s)
recmpox -i consensus.fa -o output -ref Ia,Ib
recmpox -i OZ375330.1 -o output -ref Ib,IIb   # UK recombinant case example
recmpox -i accessions.txt -o output -ref Ia,Ib   # one accession per line or comma-separated
```

**Note**: Either `-ref` (e.g. `Ia,Ib`, `IIa,IIb`, or `Ib,IIb`) or both `-ref1` and `-ref2` are required. With `-ref`, default references are used (Ia=OZ254474.1, Ib=PP601219.1, IIa=OZ287284.1, IIb=NC_063383.1).

### Command-line options

#### Required
- `-i, --input`: Input: FASTA file, directory of `.fa`/`.fasta`/`.fna`, `.txt` file of accessions (one per line or comma-separated), or NCBI accession(s). Accessions are downloaded and used as queries.

#### Reference (use one of)
- `-ref`: Reference pair: `Ia,Ib` or `IIa,IIb`. Uses built-in defaults.
- `-ref1`, `-ref2`: Custom references (path or NCBI accession). Use with `-ref1_g`/`-ref2_g` for labels (e.g. `-ref1_g Ia -ref2_g Ib`).

#### Optional
- `-o, --output`: Output directory (default: `output/`)
- `-ref1_g`, `-ref2_g`: Genotype labels for TSV/HTML (default from `-ref` or accession)
- `-include-indels`: Include diagnostic indels (default: SNPs only)
- `-min-indel-size`: Min indel length (bp) when using `-include-indels` (default: 100)
- `-m, -minor-ref-pct`: Minor reference % threshold for calling "potential recombinant" (default: 10). Increase to be more conservative (e.g. 15, 20).
- `-t, --threads`: Number of threads
- `-q, --quiet`: Log to file only

### Examples

```bash
# Custom references
recmpox -i fasta/ -o output -ref1 NC_003310.1 -ref2 PP601219.1 -ref1_g Ia -ref2_g Ib -t 4

# Mixed clades (e.g. Ia vs IIb)
recmpox -i fasta/ -o output -ref1 ACC1 -ref2 ACC2 -ref1_g Ia -ref2_g IIb

# Include diagnostic indels
recmpox -i fasta/ -o output -ref Ia,Ib -include-indels
```

## Output files

- **recmpox_results.tsv**: Per-genome counts (n_ref1, n_ref2, n_other), percentages (pct_ref1, pct_ref2, pct_other), and recombinant call (no recombinant / potential recombinant).
- **recmpox_results.html**: Interactive report (summary, sortable table, stacked bar chart, diagnostic SNP positions, diagnostic sites per sample, recombination tracts and breakpoints per sample). Split into multiple files + index when >100 genomes.
- **potential_recombinants_diagnostic_sites.tsv**: Diagnostic site classification per potential recombinant (when any exist).
- **diagnostic_snps.txt**: List of diagnostic SNP positions (ref1 vs ref2 alleles).
- **.recmpox.log**: Log file (in output directory).
- With **-extract-tracts**: **tracts/** — per-sample FASTA with only Ia tract positions (rest N) and only Ib tract positions (rest N): `Ia_recombinant_ancestral_tract.fa`, `Ib_recombinant_ancestral_tract.fa` (clade names depend on -ref).
- With **-phylogeny**: **phylogeny/** folder containing **phylogeny_alignment.fasta**, **phylogeny_tree.treefile** (midpoint-rooted), **phylogeny_tree.pdf**, and **phylogeny_tree.svg**. The pipeline runs a bundled R script (requires **ape**, **phytools**, **ggtree**, **ggplot2**). To test the R script on an existing tree (e.g. after a run where the R step failed), from your **output directory** run:
  ```bash
  Rscript /path/to/recmpox/references/root_tree_figure.R phylogeny/phylogeny_tree.treefile
  ```
  (Replace `/path/to/recmpox` with the RecMpox install or source path; the script writes **rooted.tree**, **tree_figure.pdf**, and **tree_figure.svg** into **phylogeny/**.)

## Interpretation

- **No recombinant**: One ref dominates (minor ref &lt; 10% of diagnostic sites).
- **Potential recombinant**: Both refs contribute ≥10% (minor ref % ≥ 10%). The HTML report shows recombination tracts (beginning/end of each tract) and breakpoints between tracts. A single tract means the genome is entirely one clade (no recombination).
- **High pct_other**: Many Ns, gaps, or non-ref bases at diagnostic sites (poor coverage or alignment).

## HTML Output example

Example HTML output for one sample:

![HTML Output Example](html_example.png)

Also on [**Zenodo**](https://doi.org/10.5281/zenodo.18495962) and [**Docker**](https://hub.docker.com/repository/docker/daanjansen94/recmpox/general).

## Citation

If you use RecMpox in your research, please cite:

```
Jansen, D., & Vercauteren, K. RecMpox: A Command-Line Tool for Flagging Potential Recombination Events in Monkeypox Viruses (v0.0.3). Zenodo. https://doi.org/10.5281/zenodo.18495962
```

## Acknowledgements

RecMpox integrates several external bioinformatics tools. Please also cite these tools as appropriate when using RecMpox:

- [Squirrel](https://github.com/aineniamh/squirrel) 
- [minimap2](https://github.com/lh3/minimap2) 
- [samtools](https://www.htslib.org/)

## License

This project is licensed under the GNU General Public License v3.0 (GPL-3.0) - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Support

If you encounter any problems or have questions, please open an issue on GitHub.