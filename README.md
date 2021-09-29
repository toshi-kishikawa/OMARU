<p align="center">
  <img width="500" src="https://user-images.githubusercontent.com/78953113/122951074-f5b9c300-d3b7-11eb-852c-7da153bc4ee0.jpg" />
</p>

# OMARU
OMARU (**O**mnibus **m**etagenome-wide **a**ssociation study with **r**ob**u**stness) is a new end-to-end analysis workflow that covers a wide range of microbiome analysis from phylogenetic and functional profiling to case-control metagenome-wide association studies (MWAS).

- implement rigorous quality control (QC) of shotgun sequence reads, samples, clades, and genes.
- construct phylogenetic and functional profiling of the metagenome. 
- conduct three major components of the case-control association tests of MWAS (i.e., phylogenetic, gene, and biological pathway analyses) with rigorous handling of false positives in statistical analysis.
- visualize attractive figures which enables comprehensive summary of the association test results.
- evaluate pathway-level links between the metagenome and the germline genome-wide association studies (GWAS) of the host genome. 
- link taxa and genes in the metagenome. 

## Overview
![Graphical_abstract](https://user-images.githubusercontent.com/78953113/126056587-d6cfe8d8-7164-4a2a-90b2-6e9f40e17455.png)

## Publication/Citation
Currently under development


## Requirements
OMARU requires [**Conda**](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) (version 4.8 or greater) for installation and operation.


## Installation of OMARU
To install OMARU via Conda, create a new environment using the following command:

```bash
    # Download conda-pack of OMARU into directory <OMARU_dir>
    $ mkdir -p OMARU_dir
    $ wget https://zenodo.org/record/5201539/files/OMARU.tar.gz
    
    # Activate the usage environment of downloaded conda-pack. This adds `OMARU_dir/bin` to your path
    $ tar -xzf OMARU.tar.gz -C OMARU_dir
    $ cd OMARU_dir
    $ source ./bin/activate
    $ conda-unpack
 ```
 
## Download of reference databases 

 Users can flexibly customize the reference data.
 
 Default reference databases can be downloaded and prepared for OMARU in `<OMARU_dir>/OMARU_databases`. 
 
 Activate the `OMARU` environment and then run as follows.

- **Download and prepare reference databases for read QC such as phix, adapters (in Trimmomatic), and the human genome (hg38)**
```bash
    $ Prepare_reference_read_QC.sh OMARU_dir/OMARU_databases
```
- **Download and prepare reference databases of phylogenetic analyses (based on ChocoPhlAn) in a FASTA format.**
```bash
    $ Prepare_reference_ChocoPhlAn.sh OMARU_dir/OMARU_databases
   ```
If you adopt your original phylogenetic reference data, FASTA file should be converted to the format of Bowtie reference, and the following data should be prepared in `<OMARU_dir>/OMARU_databases`.

&nbsp; 1 `NCBI_species_scaffold_<phylogenetic_reference>.txt` (refer to `NCBI_species_scaffold_EXAMPLE.txt`)  
&emsp; **Row** One scaffold in FASTA files per row   
&emsp; **Column** 1.NCBI Accession ID 2.Species 3.Scaffold
  
&nbsp; 2 `NCBI_lineage_<phylogenetic_reference>.txt` (refer to `NCBI_lineage_EXAMPLE.txt`)    
&emsp; **Row**  One FASTA file per row  
&emsp; **Column** The first column is NCBI Accession ID. The second and subsequent columns are lineages (Kingdom ~ Species).

&nbsp; 3 `eachL_lineage_<phylogenetic_reference>.txt` (refer to `eachL_lineage_EXAMPLE.txt`)  
&emsp; **Row**  One clade per row  
&emsp; **Column** The first column is the clade. The second and subsequent columns are lineages (Kingdom ~ Species).

- **Download and prepare reference databases of functional analyses (based on UniRef90 and GO term)**
```bash
    $ Prepare_reference_UniRef90.sh OMARU_dir/OMARU_databases
```
If you adopt your original functional reference data, FASTA file should be converted to the format of Bowtie reference, and the additional following data should be prepared in `<OMARU_dir>/OMARU_databases`.

&nbsp; 1 `<gene_reference>_annotatioin.txt.gz` (refer to `EXAMPLE_gene_annotation.txt`)  
&emsp; **Row** One gene per row    
&emsp; **Column** The first column is gene ID. The second and subsequent columns are metadata of genes.
  
&nbsp; 2 `header_<gene_reference>_annotatioin.txt` (refer to `header_EXAMPLE_gene_annotation.txt`)
  
&nbsp; 3 `<pathway_reference>_annotatioin.txt.gz` (refer to `EXAMPLE_pathway_annotation.txt`)  
&emsp; **Row** One pathway per row    
&emsp; **Column** The first column is pathway ID. The second and subsequent columns are metadata of pathways.
  
&nbsp; 4 `header_<pathway_reference>_annotatioin.txt` (refer to `header_EXAMPLE_pathway_annotation.txt`)


**Note**: This process can be time and resource-intensive, taking several hours, almost 200GB of free disk space.

## Create project directory <OMARU_project_dir> and set your data
Users can create a new project directory as follows:
```bash
    $ Prepare_project_dir.sh OMARU_project_dir OMARU_dir
```
Put your input data of metagenomic shotgun sequencing (FASTQ format) to predetermined folder (`<OMARU_project_dir>/data/original_fastq`) according to the following format:  
&emsp; **Name** `<Sample_ID>_R1.fastq.gz` `<Sample_ID>_R2.fastq.gz`

Put your sample list with metadata to predetermined folder (`<OMARU_project_dir>/data`) according to the following format:  
&emsp; **Name** `original_sample_list.txt`  
&emsp; **Row**  The first row is the header. The first four terms should be "Sample", "Phenotype",	"Sex", and "Age".  
&emsp;&emsp;&emsp;   One sample per row from the second row onwards.  
&emsp; **Column** The first four columns are sample ID, phenotype, gender, and age, in that order.  
&emsp;&emsp;&emsp;&emsp;&emsp; The fifth column and subsequent columns are the other metadata.  
&emsp;&emsp;&emsp;&emsp;&emsp; At phenotype column, positive and negative samples should be 1 and 0, respectively. 

Put covariate list for phylogenetic and gene association tests at `<OMARU_project_dir>/data` according to the following format:

**Name** `covariates.txt`

**Row**  One covariate per row. (Each word should be the same as the word of the header of `original_sample_list.txt`)

Customize parameters of `<OMARU_project_dir>/config.yaml` that you may want to change. 
| Parameter | Description | Default |
|:-----------:|:-----------|:------------|
|PHENOTYPE|one word of your project such as phenotype|project_phenotype|
|DB_DIR|absolute pathname of your database directory|<OAMRU_dir>/OMARU_databases|
|THREAD|number of threads used by each shell script|4|
|REF: HUMAN|reference name of the human genome|Homo_sapiens_assembly38|
|REF: PHYL|phylogenetic reference name|CHOCO|
|REF: GENE|reference name of gene|UniRef90|
|REF: PATH|reference name of pathway|GO|
|PHYL_THRESHOLD|cutoff value for relative abundance rate of clades <br> 5 indicates 1 × 10<sup>-5</sup>|5|
|SUFFIX_COV|suffix for a list file of covariates in association tests <br> See details in `Step 4-1`|\_wBMI|
|N_PCs|list of numbers of PCs as covariates to be tried in association test <br> See details in `Step 4-1`|[0,1,2,3]|
|N_PC_PHYL|number of PC finally adopted in the phylogenetic association test|2|
|N_PC_GENE|number of PC finally adopted in the gene association test|2|
|PQ|threshold of p-value and false discovery rate|0.05|
|N_SIG_CLADE|number of clades with significant differences in the phylogenetic association test|2|
|TARGETS|genes to be evaluated for links with phylogenetic data|[\"XXX\",\"YYY\"]|

## Usage
Various options are available d according to the function of `Snakemake` such as:
```bash
# run OMARU locally:
$ snakemake -s script/OMARU.sm --jobs 10
    
# perform a dry-run locally:
$ snakemake -s script/OMARU.sm --jobs 10 -n

# run OMARU on a cluster:
$ snakemake -s script/OMARU.sm --cluster qsub --jobs 20
```

For more details, see the "Executing Snakemake" section of the
[Snakemake docs](https://snakemake.readthedocs.io/en/v5.1.4/index.html).

### Step 1: Read QC

```bash
   $ cd OMARU_project_dir
   $ snakemake -s Snakefiles/OMARU_read_QC.sm 
```

#### Output
You can find QCed FASTQ files in the output directory, `<OMARU_project_dir>/result/QC/QCed_fastq`.

You can also check tables and figures of the statistical summary in the QC process at the output directory, `<OMARU_project_dir>/result/<Phenotype>_summary`.

For the next step, <ins>select the samples that have passed QC</ins>, and <ins>update the sample list</ins> with the name `QCed1_sample_list.txt` at `<OMARU_project_dir>/data` . 


### Step 2: Construct phylogenetic and functional profiling

```bash
   $ cd OMARU_project_dir
   $ snakemake -s Snakefiles/OMARU_QCed1.sm 
```

#### Output
You can find QCed FASTQ files in the output directory, `<OMARU_project_dir>/result/QC/QCed_fastq`
As in the previous step, you can check tables and figures of the statistical summary in the profiling process at the output directory, `<OMARU_project_dir>/result/<Phenotype>_summary`.

For the next step, <ins>select the samples with sufficient quality in the profiling process</ins>, and <ins>update the sample list</ins> with the name `QCed2_sample_list.txt` at `<OMARU_project_dir>/data`. 

### Step 3: Check the phylogenetic and gene abundance data

```bash
   $ cd OMARU_project_dir
   $ snakemake -s Snakefiles/OMARU_QCed2.sm 
```

#### Output

You can check tables and figures of the statistical summary of the abundance data at the output directory, `<OMARU_project_dir>/result/PHYL_QCed2/<Phylogenetic_reference>_graph_basic` and `<OMARU_project_dir>/result/PHYL_QCed2/<Gene_reference>_graph_basic`.

For the next step, <ins>select the samples with appropriate profiling data for analysis</ins>, and <ins>update the sample list</ins> with the name `QCed3_sample_list.txt` at `<OMARU_project_dir>/data`. 

### Step 4-1: Case-control association test for phylogenetic abundance data

The covariate list can be changed from the default of `<OMARU_project_dir>/data/covariates.txt`. (The default is only sex and age).  
Make the new list of covariates with the name `covariates<SUFFIX_COV>.txt` and change the parameter of `SUFFIX_COV` in `<OMARU_project_dir>/config.yaml`.

You can change the list of numbers of PCs used as covariates. For example, to try a range of four PCs to seven PCs, change the parameter of `N_PCs` in `<OMARU_project_dir>/config.yaml` to [4,5,6,7].


```bash
   $ cd OMARU_project_dir
   $ snakemake -s Snakefiles/OMARU_Phyl_AS.sm 
```

#### Output
You can find summary statistics of result of phylogenetic association tests with the name `sumstats_*.txt` at  `<OMARU_project_dir>/result/PHYL_QCed3/<Phylogenetic_reference>_association_test/<covariates>`.  
Also, figures at `<OMARU_project_dir>/result/PHYL_QCed3/<Phylogenetic_reference>_association_graph`.

For the next step, <ins>select the number of PCs to be adopted as covariates</ins>, and <ins>change the parameter of `N_PC_PHYL` in `<OMARU_project_dir>/config.yaml`. 
  

### Step 4-2: Phenotype permutation and visualization of the result in the phylogetic association tests

```bash
   $ cd OMARU_project_dir
   $ snakemake -s Snakefiles/OMARU_Phyl_permutation_visualization.sm 
```

#### Output
You can find summary statistics that integrate permutation results and annotation information with the name `sumstats_*_annot.txt` at  `<OMARU_project_dir>/result/PHYL_QCed3/<Phylogenetic_reference>_association_test/<covariates>`.  
Also, figures at `<OMARU_project_dir>/result/PHYL_QCed3/<Phylogenetic_reference>_association_graph`.

A phylogenetic tree indicating the association results is at `result/PHYL_QCed3/<Phylogenetic_reference>_ggtree/<covariates>`.

### Step 5-1: Case-control association test for gene abundance data
 
You can customize the parameters of `SUFFIX_COV` and `N_PCs` in `<OMARU_project_dir>/config.yaml` as in `Step 4-1`.

```bash
   $ cd OMARU_project_dir
   $ snakemake -s Snakefiles/OMARU_Func_AS.sm 
```
#### Output
You can find summary statistics of results of gene association tests with the name `sumstats_*.txt` at  `<OMARU_project_dir>/result/FUNC_QCed3/<Gene_reference>_association_test/<covariates>`.  
Also, figures at `<OMARU_project_dir>/<OMARU_project_dir>/result/FUNC_QCed3/<Gene_reference>_association_graph`.

For the next step, <ins>select the number of PCs to be adopted as covariates</ins>, and <ins>change the parameter of `N_PC_FUNC` in `<OMARU_project_dir>/config.yaml`. 

### Step 5-2: Phenotype permutation in gene association tests

```bash
   $ cd OMARU_project_dir
   $ snakemake -s Snakefiles/OMARU_Func_permutation.sm 
```

#### Output
You can find summary statistics that integrate permutation results and annotation information with the name `sumstats_*_annot.txt` at  `result/FUNC_QCed3/<Gene_reference>_association_test/<covariates>`.  
Also, figures at `<OMARU_project_dir>/result/FUNC_QCed3/<Gene_reference>_association_graph`.
  
### Step 5-3: Gene set enrichment analysis (GSEA) using the ranking of the genes

```bash
   $ cd OMARU_project_dir
   $ snakemake -s Snakefiles/OMARU_FUNC_GSEA.sm 
```

#### Output
You can find summary statistics of result of GSEA with the name `<Phenotype>_result_GSEA_<Gene_reference>_<Pathway_reference>_annot.txt` at  `<OMARU_project_dir>/result/FUNC_QCed3/<Gene_reference>_association_test/<covariates>/<Gene_reference>_<Pathway_reference>`.  
Also, figures at `<OMARU_project_dir>/result/FUNC_QCed3/<Gene_reference>_association_graph/<covariates>/GSEA_<covariates>`.

### Step 5-4: Links between the microbe MWAS and the germline GWAS of host
Use a tool for pathway analysis with summary statistics from GWAS (e.g., [PASKAL](https://www2.unil.ch/cbg/index.php?title=Pascal)) in order to determine pathway enrichment of the human genome in your target phenotype.
  
Put the pathway enrichment data of the human genome to predetermined folder (`<OMARU_project_dir>/data`) according to the following format:

**Name** `<Phenotype>_GWAS_<Pathway_reference>.txt`  
**Row**  The first row is header. One pathway per row from the second row onwards.  
**Column** The first column is pathway ID and the second column is p-value.

```bash
   $ cd OMARU_project_dir
   $ snakemake -s Snakefiles/OMARU_Func_MWAS_GWAS.sm 
```

#### Output
You can check tables and figures of the links between MWAS and GWAS at the output directory, `<OMARU_project_dir>/result/FUNC_QCed3/MWAS_GWAS/<Pathway_reference>`.

### Step 6: Links between taxa and genes in the metagenome.
Set the genes to be evaluated for links with phylogenetic data as the parameter of `TARGETS` in `<OMARU_project_dir>/config.yaml`.

```bash
   $ cd OMARU_project_dir
   $ snakemake -s Snakefiles/OMARU_Phy_Fun_link.sm 
```

#### Output
You can check tables and figures of the statistical summary in the QC process at the output directory, `<OMARU_project_dir>/result/PHYL_FUNC_link/<TARGETS>/<TARGETS>_result`.

## Licence
This software is freely available for academic users. Usage for commercial purposes is not allowed.  
Please refer to the [LICENCE](https://github.com/toshi-kishikawa/OMARU/blob/master/LICENSE.md) page.

## Contact
Toshihiro Kishikawa ([tkishikawa@ent.med.osaka-u.ac.jp](mailto:tkishikawa@ent.med.osaka-u.ac.jp))

