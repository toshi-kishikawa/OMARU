<p align="center">
  <img width="500" src="https://user-images.githubusercontent.com/78953113/122951074-f5b9c300-d3b7-11eb-852c-7da153bc4ee0.jpg" />
</p>

# OMARU
OMARU (**O**mnibus **m**etagenome-wide **a**ssociation study with **r**ob**u**stness) is a new end-to-end analysis workflow that covers a wide range of microbiome analysis from phylogenetic and functional profiling to case-control metagenome-wide association studies (MWAS).

- implement rigorous quality control (QC) of shotgun sequence reads, samples, clades, and genes.
- construct phylogenetic and functional profiling of the metagenome. 
- conduct three major components of the case-control association tests of MWAS (i.e., phylogenetic, gene, and biological pathway analyses) with rigorous handling of false positives in statistical analysis.
- visualize attractive figures which enables comprehensive summary of the association test results.
- evaluate pathway-level links between the metagenome and the germline genome-wide association studies (GWAS) of the host genome, 
- link taxa and genes in the metagenome. 

## Overview
![Graphical_abstract](https://user-images.githubusercontent.com/78953113/122957772-39fb9200-d3bd-11eb-86d0-f2421ab638a7.png)

## Publication/Citation
Currently under development


## Requirements
OMARU requires [**Conda**](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) (version 4.8 or greater) for installation and operation.


## Installation of OMARU
To install OMARU via Conda, create a new environment using the following command:

```bash
    # Download conda-pack of OMARU.
    conda create -n metalaffa metalaffa -c bioconda -c borenstein-lab
   
   # Unpack environment into directory <OMARU_dir>
    mkdir -p OMARU_dir
    git clone https://github.com/toshi-kishikawa/OMARU
    tar -xzf OMARU.tar.gz -C OMARU_dir

    # Activate the environment. This adds `OMARU_dir/bin` to your path
    source OMARU_dir/bin/activate
 ```
 
## Download of reference databases 

 Users can flexibly customize the reference data. Default reference databases can be downloaded and prepared for OMARU using scripts in `OMARU_dir/OMARU_databases` as follows. These databases will be installed in `OMARU_dir/OMARU_databases`. Activate the `OMARU` environment and then run:
```bash
    # Download and prepare reference databases for read QC such as phix, adapters (in Trimmomatic), and human genome (hg38)
    Prepare_reference_read_QC.sh OMARU_dir/OMARU_databases
```
```bash
    # Download and prepare reference databases of phylogenetic analyses (based on ChocoPhlAn) in a FASTA format.
    Prepare_reference_ChocoPhlAn.sh OMARU_dir/OMARU_databases
   ```
if you adopt your original reference data, the additional following data should be prepared in `OMARU_dir/OMARU_databases`.
  1 `NCBI_species_scaffold_<phylogenetic_reference>.txt` (Refer to `NCBI_species_scaffold_EXAMPLE.txt`)
  **Row**  One scaffold in FASTA files per row
  **Column** 1.NCBI_Accession ID 2.Species 3.Scaffold
  
  2 `NCBI_lineage_<phylogenetic_reference>.txt` (Refer to `NCBI_lineage_EXAMPLE.txt`)
  **Row**  One FASTA file per row
  **Column** 1.NCBI_Accession ID 2~8.Kingdom ~ Species

  3 `eachL_lineage_<phylogenetic_reference>.txt` (Refer to `eachL_lineage_EXAMPLE.txt`) 
  **Row**  One clade per row
  **Column** 1.Clade 2~8.Kingdom ~ Species


    
    # Download and prepare reference databases of functional analyses (based on UniRef90 and GO term)
    Prepare_reference_UniRef90.sh OMARU_dir/OMARU_databases

**Note**: This process can be time and resource intensive, taking several hours, ~108GB of free disk space, and ~40GB of RAM.

By default, MetaLAFFA is able to interface with Sun Grid Engine (SGE) and HTCondor clusters. This is achieved via the use of Python job submission wrapper scripts, included in the `$CONDA_PREFIX/MetaLAFFA/src/` directory (`$CONDA_PREFIX/MetaLAFFA/src/sge_submission_wrapper.py` and `$CONDA_PREFIX/src/condor_submission_wrapper.py` respectively). If your cluster uses a different cluster management system, then you will need to create your own job submission wrapper by following these steps:

## Create project directory <OMARU_project_dir> and set your data
Users can create a new project directory as follows:

    prepare_project_dir.sh OMARU_project_dir OMARU_dir/OMARU_scripts

Put your input data of metagenomic shotgun sequencing (FASTQ format) to predetermined folder (`OMARU_project_dir/data/original_fastq`) according to the following format:

**Name** `<Sample_ID>_R1.fastq.gz` `<Sample_ID>_R2.fastq.gz`

Put your sample list with metadata to predetermined folder (`OMARU_project_dir/data`) according to the following format:

**Name** `original_sample_list.txt`

**Row**  One sample per row

**Column** The first three columns are sample ID, gender, age, and other metadata from the fourth column onwards.

Arrange some parameters of `OMARU_project_dir/config.yaml` that you may want to change. 


## Usage

Snakemakeに準じた様々なオプションが使える。
このページを参照にしてください。 http://
下記にいくつか例を挙げる
If you are running MetaLAFFA locally, use:

    ./MetaLAFFA.py

If you are running MetaLAFFA on a cluster, use:

    ./MetaLAFFA.py --use-cluster

-n
-j


This option will use the submission wrapper and jobscript specified in your project's `config/cluster.py` configuration submodule to submit jobs to your cluster system.

### Step 1: Read QC

```bash
    cd OMARU_project_dir
    snakemake -s OMARU_read_QC.sm 
```

#### output
You can check tables and figures of the statistical summary in the QC process at the output directory, `OMARU_project_dir/result/<Phenotype>_summary`.

Select the samples that has passed QC, and update the sample list with the name `QCed1_sample_list.txt` at `OMARU_project_dir/data`. 


### Step 2: Construct phylogenetic and functional profiling

```bash
    cd OMARU_project_dir
    snakemake -s OMARU_QCed1.sm 
```

#### output
As in the previous step, you can check tables and figures of the statistical summary in the profiling process at the output directory, `OMARU_project_dir/result/<Phenotype>_summary`.

Select the samples with sufficient quality in the profiling process, and update the sample list with the name `QCed2_sample_list.txt` at `OMARU_project_dir/data`. 

### Step 3: Check the phylogenetic and gene abundance data

```bash
    cd OMARU_project_dir
    snakemake -s OMARU_QCed2.sm 
```

#### output
You can check tables and figures of the statistical summary of the abundance data at the output directory, `OMARU_project_dir/result/PHYL_QCed2/<Phylogenetic_reference>_graph_basic` and `OMARU_project_dir/result/PHYL_QCed2/<Gene_reference>_graph_basic`.

Select the samples with appropriate profiling data for analysis, and update the sample list with the name `QCed3_sample_list.txt` at `OMARU_project_dir/data`. 

### Step 4-1: Case-control association test for phylogenetic abundance data

```bash
    cd OMARU_project_dir
    snakemake -s OMARU_Phyl_AS.sm 
```

#### output
You can check tables and figures of the statistical summary in the QC process at the output directory, `OMARU_project_dir/result/<Phenotype>_summary`.

 

### Step 4-2: Phenotype permutation and visualization of the result in phylogetic association test


```bash
    cd OMARU_project_dir
    snakemake -s OMARU_read_QC.sm 
```

#### output
You can check tables and figures of the statistical summary in the QC process at the output directory, `OMARU_project_dir/result/<Phenotype>_summary`.

Select the samples that has passed QC, and update the sample list with the name `QCed1_sample_list.txt` at `OMARU_project_dir/data`. 

### Step 5-1: Case-control association test for gene abundance data

```bash
    cd OMARU_project_dir
    snakemake -s OMARU_read_QC.sm 
```

#### output
You can check tables and figures of the statistical summary in the QC process at the output directory, `OMARU_project_dir/result/<Phenotype>_summary`.

Select the samples that has passed QC, and update the sample list with the name `QCed1_sample_list.txt` at `OMARU_project_dir/data`. 

### Step 5-2: Phenotype permutation in gene association test

```bash
    cd OMARU_project_dir
    snakemake -s OMARU_read_QC.sm 
```

#### output
You can check tables and figures of the statistical summary in the QC process at the output directory, `OMARU_project_dir/result/<Phenotype>_summary`.

Select the samples that has passed QC, and update the sample list with the name `QCed1_sample_list.txt` at `OMARU_project_dir/data`. 

### Step 5-2: Phenotype permutation in gene analysis

```bash
    cd OMARU_project_dir
    snakemake -s OMARU_read_QC.sm 
```

#### output
You can check tables and figures of the statistical summary in the QC process at the output directory, `OMARU_project_dir/result/<Phenotype>_summary`.

Select the samples that has passed QC, and update the sample list with the name `QCed1_sample_list.txt` at `OMARU_project_dir/data`. 

### Step 5-3: Gene set enrichment analysis

```bash
    cd OMARU_project_dir
    snakemake -s OMARU_read_QC.sm 
```

#### output
You can check tables and figures of the statistical summary in the QC process at the output directory, `OMARU_project_dir/result/<Phenotype>_summary`.

Select the samples that has passed QC, and update the sample list with the name `QCed1_sample_list.txt` at `OMARU_project_dir/data`. 

### Step 5-4: Comparison of the pathway analysis results between the MWAS and the GWAS

```bash
    cd OMARU_project_dir
    snakemake -s OMARU_read_QC.sm 
```

#### output
You can check tables and figures of the statistical summary in the QC process at the output directory, `OMARU_project_dir/result/<Phenotype>_summary`.

Select the samples that has passed QC, and update the sample list with the name `QCed1_sample_list.txt` at `OMARU_project_dir/data`. 

### Step 6: Link of phylogenetic and gene abundance data

```bash
    cd OMARU_project_dir
    snakemake -s OMARU_read_QC.sm 
```

#### output
You can check tables and figures of the statistical summary in the QC process at the output directory, `OMARU_project_dir/result/<Phenotype>_summary`.

Select the samples that has passed QC, and update the sample list with the name `QCed1_sample_list.txt` at `OMARU_project_dir/data`. 









Note: Any configuration changes made in the configuration module located at $CONDA_PREFIX/MetaLAFFA/config will be the default configurations for any newly created projects. Thus, if you have custom settings that you think should be preset in any new projects, you should make those changes to this base configuration module.

1) If you input ICD-10 code (for example, F20 for Schizophrenia as below), you will also get gold-standard drug (approved drugs for ICD-10 F20 in ChEMBL and TTD [Therapeutic Target Database]) in an output Q-Q plot data.
ICD-10 codes which are not listed in ChEMBL and TTD are not applicable.
The example command is as follows;
```bash
cd ./script
./Trans-Phar.sh F20
```

or 2) If you need not get gold-standard Q-Q plot, you only enter the example command as follows;
```bash
cd ./script
./Trans-Phar.sh
```

### Step 2: Put your input data to predetermined folder (named as Input_GWASsummary)
## Output

1) The example TWAS result outputs are as follows (if you use tutorial GWAS data);

```bash
#TWAS results according to each 29 GTEx (v7) tissue and combined files from all 29 tissues at Output/Schizo/TWASresults.

#For Example
cd ../Output/Schizo/TWASresults/ALLTISSUE
less GTEx_Adipose_Subcutaneous.chr_all.focus_shaped.tsv #TWAS result file (shaped), file format is described in https://github.com/bogdanlab/focus/wiki/Fine-mapping-TWAS-associations
#TWAS result png files are also in Output/Schizo/TWASresults/ALLTISSUE
```

2) The example Spearman result outputs are as follows (if you use tutorial GWAS data);

```bash
#Output p-values for Negative Spearmans's correlation tests according to total 308,872 pairs of TWAS tissue - CMap cell - Compunds
#Data of TWAS tissue - CMap cell - Compunds whose P-value < 0.0001 are in Output/Schizo/Spearmanresults/spearman_eachpair_results and Output/Schizo/Spearmanresults/spearman_eachpair_coplots
#For Example
cd ../../Spearmanresults/spearman_totalresults
less ALLpairs_spearmanresults.txt
#Q-Q plot for distribution of these P-value is also in Output/Schizo/Spearmanresults/spearman_totalresults
```

## Acknowledgements
* The original [FOCUS](https://github.com/bogdanlab/focus) was written by Nicholas Mancuso et al.

## Licence
This software is freely available for academic users. Usage for commercial purposes is not allowed.
Please refer to the [LICENCE](https://github.com/konumat/Trans-Phar/blob/master/LICENSE.md/LICENSE.md) page.

