# monq_adaptive_genomics

This repository contains the bioinformatic and statistical scripts used in a study of adaptive genomics in Montezuma quail (*Cyrtonyx montezumae*), combining pangenomics, landscape genomics, ecological niche modeling, and forward-time simulation. All analyses were run on the Purdue University computing clusters (Bell and Negishi) using SLURM.

---

## Repository Structure

```
monq_adaptive_genomics/
├── analysis/
│   ├── external_scripts/
│   │   ├── adaptive_index.R
│   │   ├── genomic_offset.R
│   │   └── rdadapt.R
│   ├── monq_admixture.sh
│   ├── monq_gea-go.R
│   ├── monq_genomic_load_VEP.sh
│   ├── monq_heterozygosity.sh
│   ├── monq_kraken2.sh
│   ├── monq_mapping.sh
│   ├── monq_qc.sh
│   ├── monq_ref_formating.sh
│   ├── monq_snp.sh
│   ├── monq_snpcall_byID.sh
│   ├── monq_sv.sh
│   └── monq_vg.sh
├── enm/
│   ├── MONQ_ENM_workflow.txt
│   ├── monq_enm.R
│   ├── monq_env_processing.ipynb
│   └── monq_occ_processing.R
├── pangenome/
│   ├── monq_longread_assembly.sh
│   ├── monq_longread_kraken2.sh
│   └── monq_pangenome_assembly.sh
├── plotting/
│   ├── MONQ_load_plots.R
│   ├── MONQ_slim_plots_refined.R
│   └── MONQ_slim.gnx_plotting.R
└── simulation/
    ├── utilities/
    │   ├── add_T_at_segment_ends.py
    │   ├── concat_fasta_fromVCF.py
    │   ├── make_dummy_fasta_fromVCF_forSLiM.py
    │   ├── MONQ_langen_cc_prelocation.R
    │   └── recode_gt_ge2_to1.py
    ├── MONQ_langen_cc_postburnin_resume2.slim
    ├── MONQ_langen_cc_preburnin.slim
    ├── MONQ_langen_cc_preburnin_resume.slim
    ├── MONQ_rescue.slim
    ├── MONQ_rescue_burnin.slim
    ├── MONQ_VEP_forSLiM.sh
    ├── monq_geonomics.py
    ├── monq_gnx.py
    └── monq_gnx_params.py
```

---

## Folder Descriptions

### `pangenome/`

Contains scripts for assembling individual-level genomes and constructing the Montezuma quail pangenome.

Long-read (PacBio HiFi) assemblies are generated for ten individuals using `hifiasm` by `monq_longread_assembly.sh`, duplicate haplotypes are removed using `purge_dups`, and assembled contigs are screened for contamination with Kraken2 (`monq_longread_kraken2.sh`). The resulting assemblies, after renaming and repeat-masking (described in a separate pipeline), are used to build the pangenome using the Minigraph-Cactus workflow (`monq_pangenome_assembly.sh`). The final pangenome graph (`.gfa` / `.vg`) serves as the reference for short-read variant calling in the `analysis/` folder (`monq_vg.sh`).

---

### `analysis/`

Contains the core genomic analysis pipeline, covering everything from raw read processing through population genetics and genotype–environment association (GEA) analysis.

The workflow begins with quality control and contamination filtering of short-read sequencing data (`monq_qc.sh`, `monq_kraken2.sh`), followed by reference genome preparation (`monq_ref_formating.sh`). Reads are then mapped to the pangenome using `vg giraffe` and processed per-chromosome (`monq_vg.sh`). SNP calling is performed genome-wide using GATK HaplotypeCaller, with a parallelized per-individual, per-chromosome approach orchestrated by `monq_snpcall_byID.sh`. Structural variants (SVs) are called using the `vg call` framework and merged, filtered, and profiled in `monq_sv.sh`. Population-level SNP filtering, LD pruning, and PCA are handled in `monq_snp.sh` using ANGSD, ngsLD, and PCAngsd.

Population structure is further assessed via admixture analysis with NGSadmix (`monq_admixture.sh`), and individual heterozygosity is estimated from site allele frequency spectra (`monq_heterozygosity.sh`).

Genomic load is quantified using Ensembl VEP with a combined Japanese quail and chicken W-chromosome annotation (`monq_genomic_load_VEP.sh`), classifying variants as HIGH, MODERATE, LOW, or Modifier impact.

GEA and genomic offset analyses are conducted in `monq_gea-go.R`, which integrates genotype data (outlier SVs identified via partial RDA), population structure (PCAngsd), and bioclimatic + topographic variables (USGS/CHELSA). This script calls three external R functions (from Capblancq & Forester, 2021) sourced from `external_scripts/`:

- **`rdadapt.R`** — performs the RDA-based genome scan (Mahalanobis distance test with q-value correction)
- **`adaptive_index.R`** — projects adaptive genetic composition onto the landscape under current climate
- **`genomic_offset.R`** — computes per-pixel genomic offset between current and future climate scenarios (SSP1-2.6 and SSP5-8.5, periods 2011–2040, 2041–2070, and 2071–2100)

Genotype-phenotype regression and reverse offset analyses are included in `monq_gea-go.R` script, too.
Gene Ontology (GO) enrichment and KEGG pathway analysis of GEA outlier SVs are also performed in `monq_gea-go.R` using `goseq`.

`monq_vg.sh` additionally handles all preprocessing steps needed for downstream SLiM simulations: concatenating SNP and SV VCFs, recoding multi-allelic genotypes to biallelic (calling `simulation/utilities/recode_gt_ge2_to1.py`), generating population-specific VCF subsets, constructing a dummy FASTA for SLiM input (calling `simulation/utilities/make_dummy_fasta_fromVCF_forSLiM.py` and `simulation/utilities/add_T_at_segment_ends.py`), running VEP, and estimating contemporary effective population sizes with currentNe2.

---

### `simulation/`

Contains forward-time simulation scripts modeling landscape-level population dynamics and genetic rescue under climate change.
SLiM simulations (.slim scripts) consist of two separate but related simulation sets.

Climate change simulations track the long-term demographic and genetic consequences of habitat loss driven by projected climate change, using ENM-derived habitat suitability rasters as inputs:

**Climate change simulations**

1. **First / Resumed burnin** (`MONQ_langen_cc_preburnin.slim`, `MONQ_langen_cc_preburnin_resume.slim`) — initializes populations on the ENM raster and runs to genomic equilibrium under current habitat conditions
2. **Initial / Resumed main run with climate change** (`MONQ_langen_cc_postburnin2.slim`, `MONQ_langen_cc_postburnin_resume2.slim`) — simulates population declines under projected habitat loss from ENM outputs, tracking genetic load (potential and realized), ROH, heterozygosity, mean fitness, and other population genetic statistics over time

Rescue simulations model the genetic and demographic effects of translocating individuals from a large, genetically diverse donor population into a small, isolated recipient population. Three populations are represented: Central Texas (CTX; recipient), New Mexico (NM; donor), and southeastern Mexico (SEMX; donor):

1. **Burnin** (`MONQ_rescue_burnin.slim`) — initializes the three populations from empirical VCF data and runs to demographic equilibrium
2. **Genetic rescue** (`MONQ_rescue.slim`) — introduces migrants from a donor population (NM or SEMX) into CTX at configurable intervals, tracking migrant ancestry (genome-wide, at neutral A-sites, and at adaptive G-sites marked by GEA outlier SVs), genetic load (potential and realized, by impact class), ROH, heterozygosity, and mean fitness over 100 post-rescue generations

`MONQ_VEP_forSLiM.sh` prepares VEP-annotated variant effect files (SIFT scores) as input for parameterizing selection coefficients in the SLiM models.

**Geonomics simulations** (`monq_geonomics.py`, `monq_gnx.py`, `monq_gnx_params.py`) use the Geonomics framework to simulate spatially explicit landscape genomics scenarios, with parameters defined in `monq_gnx_params.py`.

**Utility scripts** in `simulation/utilities/` are called upstream (primarily from `analysis/monq_vg.sh`) to prepare SLiM-compatible input files:

- **`make_dummy_fasta_fromVCF_forSLiM.py`** — builds a concatenated dummy FASTA from a `.fai` index and VCF, placing REF alleles at variant sites and `C` elsewhere; outputs a segment coordinate file
- **`add_T_at_segment_ends.py`** — appends `T` nucleotides at chromosome segment boundaries in the dummy FASTA to ensure proper recombination breakpoints in SLiM
- **`recode_gt_ge2_to1.py`** — recodes multi-allelic genotype calls (allele index ≥ 2) to 1 in a VCF, making it compatible with SLiM's biallelic expectation
- **`concat_fasta_fromVCF.py`** — an alternative FASTA builder used in earlier pipeline versions
- **`MONQ_langen_cc_prelocation.R`** — prepares spatial location and habitat suitability inputs for the landscape change simulation

---

### `enm/`

Contains scripts for ecological niche modeling (ENM) of Montezuma quail under current and future climate scenarios.

Occurrence data are downloaded and cleaned from GBIF and iNaturalist using `monq_occ_processing.R` (via the `rgbif`, `rinat`, and `CoordinateCleaner` packages). Environmental variable processing and spatial operations (e.g., raster clipping, correlation filtering, bias layer generation, and resampling) are described in `MONQ_ENM_workflow.txt` and partly handled in a python notebook (`monq_env_processing.ipynb`) using ArcGIS Pro. Model fitting and projection are performed in `monq_enm.R` using `biomod2` and `ENMeval`, benchmarking multiple algorithms (MaxEnt, Random Forest, GBM, GLM, ANN, and others) and selecting a Random Forest ensemble model for projections across present and future (SSP126 and SSP585; 2040/2070/2100). The resulting habitat suitability rasters feed into the `enm/` outputs used by the simulation scripts.

---

### `plotting/`

Contains R scripts for generating publication-quality figures from downstream analysis results.

`MONQ_load_plots.R` visualizes genomic load across SNPs, all SVs, and GEA outlier SVs, producing stacked bar plots of variant impact composition and boxplots of deleterious variant proportions by variant type and population (MX, WTX, CTX). Statistical comparisons use Kruskal-Wallis and Wilcoxon tests.

`MONQ_slim_plots_refined.R` and `MONQ_slim.gnx_plotting.R` visualize outputs from the SLiM genetic rescue simulations and Geonomics landscape simulations, respectively, including trajectories of population size, fitness, genetic load, heterozygosity, ROH, and migrant ancestry proportion over time.

---

## Notes

- All SLURM scripts were developed for Purdue's Bell and Negishi HPC clusters. Paths (e.g., `/scratch/bell/jeon96/`, `/scratch/negishi/allen715/`) and module names will need to be updated for use on other systems.
- The `analysis/external_scripts/` folder contains R functions sourced from published workflows for RDA-based landscape genomics. Please cite the original source (Capblancq & Forester, 2021) if you use them.
- Effective population size estimates used in simulations (CTX: ~132, NM: ~36,417, SEMX: ~434) were derived using currentNe2 and prior literature (Mathur & DeWoody 2021).
