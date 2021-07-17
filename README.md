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
    # Download conda-pack of OMARU into directory <OMARU_dir>
    mkdir -p OMARU_dir
    git clone https://github.com/toshi-kishikawa/OMARU
    tar -xzf OMARU.tar.gz -C OMARU_dir

    # Activate the environment. This adds `OMARU_dir/bin` to your path
    source OMARU_dir/bin/activate
 ```
 
## Download of reference databases 

 Users can flexibly customize the reference data.
 
 Default reference databases can be downloaded and prepared for OMARU in `<OMARU_dir>/OMARU_databases`. 
 
 Activate the `OMARU` environment and then run as follows.

- **Download and prepare reference databases for read QC such as phix, adapters (in Trimmomatic), and human genome (hg38)**
```bash
    $ Prepare_reference_read_QC.sh OMARU_dir/OMARU_databases
```
- **Download and prepare reference databases of phylogenetic analyses (based on ChocoPhlAn) in a FASTA format.**
```bash
    $ Prepare_reference_ChocoPhlAn.sh OMARU_dir/OMARU_databases
   ```
if you adopt your original phylogenetic reference data, FASTA file should be converted to the format of bowtie reference, and the following data should be prepared in `<OMARU_dir>/OMARU_databases`.

&nbsp; 1 `NCBI_species_scaffold_<phylogenetic_reference>.txt` (refer to `NCBI_species_scaffold_EXAMPLE.txt`)

&emsp; **Row** One scaffold in FASTA files per row
  
&emsp; **Column** 1.NCBI Accession ID 2.Species 3.Scaffold
  
&nbsp; 2 `NCBI_lineage_<phylogenetic_reference>.txt` (refer to `NCBI_lineage_EXAMPLE.txt`)　　
  
&emsp; **Row**  One FASTA file per row
  
&emsp; **Column** 1.NCBI Accession ID 2~8.Kingdom ~ Species

&nbsp; 3 `eachL_lineage_<phylogenetic_reference>.txt` (refer to `eachL_lineage_EXAMPLE.txt`) 　　

&emsp; **Row**  One clade per row

&emsp; **Column** 1.Clade 2~8.Kingdom ~ Species

- **Download and prepare reference databases of functional analyses (based on UniRef90 and GO term)**
```bash
    $ Prepare_reference_UniRef90.sh OMARU_dir/OMARU_databases
```
if you adopt your original functional reference data, FASTA file should be converted to the format of bowtie reference, and the additional following data should be prepared in `<OMARU_dir>/OMARU_databases`.

&nbsp; 1 `<gene_reference>_annotatioin.txt.gz` (refer to `EXAMPLE_gene_annotation.txt`)

&emsp; **Row** One gene per row
  
&emsp; **Column** 1.Gene ID 2~. Metadata of gene
  
&nbsp; 2 `header_<gene_reference>_annotatioin.txt` (refer to `header_EXAMPLE_gene_annotation.txt`)
  
&nbsp; 3 `<pathway_reference>_annotatioin.txt.gz` (refer to `EXAMPLE_pathway_annotation.txt`)

&emsp; **Row** One pathway per row
  
&emsp; **Column** 1.pathway ID 2~. Metadata of pathway
  
&nbsp; 2 `header_<pathway_reference>_annotatioin.txt` (refer to `header_EXAMPLE_pathway_annotation.txt`)


**Note**: This process can be time and resource intensive, taking several hours, almost 200GB of free disk space.

## Create project directory <OMARU_project_dir> and set your data
Users can create a new project directory as follows:

    prepare_project_dir.sh OMARU_project_dir OMARU_dir/OMARU_scripts

Put your input data of metagenomic shotgun sequencing (FASTQ format) to predetermined folder (`<OMARU_project_dir>/data/original_fastq`) according to the following format:

**Name** `<Sample_ID>_R1.fastq.gz` `<Sample_ID>_R2.fastq.gz`

Put your sample list with metadata to predetermined folder (`<OMARU_project_dir>/data`) according to the following format:

**Name** `original_sample_list.txt`

**Row**  One sample per row

**Column** The first three columns are sample ID, gender, age, and other metadata from the fourth column onwards.

Arrange some parameters of `<OMARU_project_dir>/config.yaml` that you may want to change. 
| Parameter | Description | Default |
|:-----------:|:-----------|:------------|
|PHENOTYPE|one word of your project such as phenotype|project_phenotype|
|DB_DIR|absolute pathname of your database directory|/OAMRU/OMARU_databases|
|THREAD|number of threads used by each shell script|4|
|REF: HUMAN|reference name of human genome|Homo_sapiens_assembly38|
|REF: PHYL|phylogenetic reference name|CHOCO|
|REF: GENE|reference name of gene|UniRef90|
|REF: PATH|reference name of pathway|GO|
|PHYL_THRESHOLD|cutoff value for relative abundance rate of clades <br> 5 indicates 1 × 10<sup>-5</sup>|5|
|SUFFIX_COV|suffixes for list file of covariates in association tests <br> See details in `Step 2`|\_wBMI|
|N_PCs|list of the number of PCs evaluated in association test|[0,1,2,3]|
|N_PC_PHYL|number of PC finally adopted in the phylogenetic association test|2|
|N_PC_GENE|number of PC finally adopted in the gene association test|2|
|PQ|threshold of p-value and false discovery rate|0.05|
|N_SIG_CLADE|number of clades with significant differences in the phylogenetic association test|2|
|TARGETS|genes to be evaluated for link with phylogenetic data|[\"XXX\",\"YYY\"]|

## Usage

various options can be used according to Snakemake such as:

    # run OMARU locally:
    snakemake -s script/OMARU.sm --jobs 10
    
    # perform a dry-run locally
    snakemake -s script/OMARU.sm --jobs 10 -n

    # run OMARU on a cluster:
    snakemake -s script/OMARU.sm --cluster qsub --jobs 20


For more details, see the "Executing Snakemake" section of the
[Snakemake docs](https://snakemake.readthedocs.io/en/v5.1.4/index.html).

### Step 1: Read QC

```bash
    cd OMARU_project_dir
    snakemake -s OMARU_read_QC.sm 
```

#### output
You can find QCed FASTQ files in the output directory, `<OMARU_project_dir>/result/QC/QCed_fastq`

You can also check tables and figures of the statistical summary in the QC process at the output directory, `<OMARU_project_dir>/result/<Phenotype>_summary`.

For the next step, <ins>select the samples that has passed QC</ins>, and <ins>update the sample list</ins> with the name `QCed1_sample_list.txt` at `<OMARU_project_dir>/data` . 


### Step 2: Construct phylogenetic and functional profiling

```bash
    cd OMARU_project_dir
    snakemake -s <OMARU>.sm 
```

#### output

You can find QCed FASTQ files in the output directory, `<OMARU_project_dir>/result/QC/QCed_fastq`
As in the previous step, you can check tables and figures of the statistical summary in the profiling process at the output directory, `<OMARU_project_dir>/result/<Phenotype>_summary`.

For the next step, <ins>select the samples with sufficient quality in the profiling process</ins>, and <ins>update the sample list</ins> with the name `QCed2_sample_list.txt` at `<OMARU_project_dir>/data`. 

### Step 3: Check the phylogenetic and gene abundance data

```bash
    cd OMARU_project_dir
    snakemake -s OMARU_QCed2.sm 
```

#### output
You can check tables and figures of the statistical summary of the abundance data at the output directory, `<OMARU_project_dir>/result/PHYL_QCed2/<Phylogenetic_reference>_graph_basic` and `<OMARU_project_dir>/result/PHYL_QCed2/<Gene_reference>_graph_basic`.

For the next step, <ins>select the samples with appropriate profiling data for analysis</ins>, and <ins>update the sample list</ins> with the name `QCed3_sample_list.txt` at `<OMARU_project_dir>/data`. 

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

