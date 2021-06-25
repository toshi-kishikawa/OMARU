<p align="center">
  <img width="500" src="https://user-images.githubusercontent.com/78953113/122951074-f5b9c300-d3b7-11eb-852c-7da153bc4ee0.jpg" />
</p>

# OMARU
OMARU (**O**mnibus **m**etagenome-wide **a**ssociation study with **r**ob**u**stness) is a new end-to-end analysis workflow that covers a wide range of microbiome analysis from phylogenetic and functional profiling to case-control metagenome-wide association studies (MWAS).

- implementation of rigorous quality control (QC) of shotgun sequence reads, samples, clades, and genes, OMARU constructs phylogenetic and functional profiling of the metagenome, the two main analytical pipelines. Three major components of the case-control association tests of MWAS (i.e., phylogenetic, gene, and biological pathway analyses) are subsequently conducted with rigorous handling of false positives in statistical analysis. OMARU visualizes attractive figures which enables comprehensive summary of the association test results. Furthermore, OMARU can evaluate pathway-level links between the metagenome and the germline genome-wide association studies (GWAS) of the host genome, as well as the links between taxa and genes in the metagenome. OMARU is a flexible and extensible workflow that can be customized, such as adding an up-to-date database. 

## Overview
![Graphical_abstract](https://user-images.githubusercontent.com/78953113/122957772-39fb9200-d3bd-11eb-86d0-f2421ab638a7.png)

## Publication/Citation



## Requirements
- R
- dichromat (R package) (installing by install.packages("dichromat"))
- python 3.X
- scipy
- numpy
- pandas
- math
- cycler
- kiwisolver
- FOCUS (Fine-mapping Of CaUsal gene Sets) as TWAS software


## Installation (Trans-Phar)
In order to get started with **Trans-Phar**, you can just clone this repo as follows;
```bash
git lfs clone https://github.com/konumat/Trans-Phar.git
cd ./Trans-Phar

#unzip QCed Cmap L1000 data
cd ./Cmap_QCeddata
for filename in $( ls *.gz ); do
echo ${filename}
gunzip ${filename}
done

cd ../
```

## Installation (FOCUS)
You have to install [FOCUS (Fine-mapping Of CaUsal gene Sets) soft ware] (https://github.com/bogdanlab/focus) as follows.
For detailed explanations, please visit [the original repository and installing tutorial] (https://github.com/bogdanlab/focus) and [wiki] (https://github.com/bogdanlab/focus/wiki).

When installing FOCUS, please make focus folder under the Trans-Phar folder.

```bash
git clone https://github.com/bogdanlab/focus.git
cd ./focus
python setup.py install
cd ../
```

## Usage
### Step 1: Prepare your input 
All you need is a text file with GWAS summary statistics. (A file extension is .sumstats)

| Column | Column name | Descriptions |
|:-----------:|:-----------:|:------------|
|1|CHR|Chromosome|
|2|SNP|rsID|
|3|BP|BP position|
|4|A1|Effect allele|
|5|A2|Other allele|
|6|MAF|Minor allele frequency (*optional*)|
|7|N|#Samples|
|8|BETA|Beta (effect allele)|
|9|P|P-value|

Please have a look at an example input at `./tutorial_input/Schizo.sumstats`.


### Step 2: Put your input data to predetermined folder (named as Input_GWASsummary)


```bash
mkdir ./Input_GWASsummary
mkdir ./Input_GWASsummary_done
mkdir ./Output

#if you use tutorial GWAS summary data;
gunzip ./tutorial_input/Schizo.sumstats.gz
cp ./tutorial_input/Schizo.sumstats ./Input_GWASsummary
```

### Step 3: Trans-Phar from GWAS summary to chemical compounds in all-in-one script

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

