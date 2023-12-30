 

# ACMGA

ACMGA is a reference-free multiple-genome alignment pipeline. The simplified schema of the pipeline is shown below.
![ACMGA](https://github.com/HFzzzzzzz/ACMGA/raw/master/workflow/image/schematic.jpg)





# Building Environment
ACMGA requires Python version 3.10 along with Biopython libraries. If you did not install have Biopython, please ensure to install it.
```
pip install biopython
```
## ACMGA supports building the environment either locally or using a docker image.
### Building the local environment:
- python3.10
- [Snakemake(>6.0)](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- [AnchorWave](https://github.com/baoxingsong/AnchorWave)
- [Cactus](https://github.com/ComparativeGenomicsToolkit/cactus)
- [SAMtools](http://www.htslib.org/)
- [Minimap2](https://github.com/lh3/minimap2)
- [bedtools](https://github.com/arq5x/bedtools2)
- [bedToGenePred](https://github.com/ENCODE-DCC/kentUtils/tree/master/src/hg/bedToGenePred)
- [genePredToGtf](https://github.com/ENCODE-DCC/kentUtils/tree/master/src/hg/genePredToGtf)
- [gffread](https://github.com/gpertea/gffread)
- [K8](https://github.com/attractivechaos/k8)
- [last](https://github.com/UCSantaCruzComputationalGenomicsLab/last/tree/master)

Using this approach, slight modifications to some of the paths within `command.sh` are necessary.
### Using docker image:
ACMGA currently relies on snakemake (>6.0.0), docker, and singularity. Please make sure these dependencies are installed before running ACMGA. We strongly recommend using this approach.



# Quickstart
For a quickstart with your own data, you can follow the instructions below. We recommend testing the pipeline with our test data first (see section [**Testing the pipeline**](#section2)), to ensure the pipeline will work correctly.

To get started, clone this repository.
```
git clone https://github.com/HFzzzzzzz/ACMGA.git
```
You can now prepare the run with the pipeline by doing the following:
 1.  Placing your FASTA sequences, Gff files (suffixed with  `.gff3`), and a  guide tree  into  `ACMGA/data`
 2.  Placing the CDS sequences set from all the input genomes into `ACMGA/data`
 3. For example
 
	 3.1 Copying  `ACMGA/config/config.yaml` to `ACMGA/config/myconfig.yaml` 
	 
	 3.2 Editing the `ACMGA/config/myconfig.yaml` to include :
	 
	-  Input FASTA sequences name (parameter  `fasta:` )
	-  Input GFF file name and ancestral GFF file name (parameter  `gff3:` )
	-  Path for the collection of CDS (parameter `nonDuplicateCDS:` ), using this [script](https://github.com/HFzzzzzzz/ACMGA/blob/master/workflow/scripts/CombineCDS.py) to merge CDS files and obtain `non_duplicate_CDS.fa`
	-  Path of the FASTA and the GFF files (parameter  `path:` )
	-  Species name (parameter  `species:` )
	-  Ancestor name (parameter `ancestor:` )
	-  Path of guide tree (parameter `Tree:` ), generated using recommended [steps](#section1)



The pipeline can then be executed from the  `ACMGA/`  directory in two steps. 
```
cd ACMGA
snakemake  -j 5 --configfile config/myconfig.yaml   --use-singularity  --singularity-args "-B $(pwd)"
```
1.The first step generates the `command.sh` script in the data file.

2.The second step is to enter the docker environment and run `command.sh`

```
docker login
docker run -v $(pwd):/data --rm -it mgatools/acmga:1.0
sh command.sh
```

# <a name="section2">Testing the pipeline</a>
## Building Environment using Docker

### 1、Create a conda environment named "acmga" with Python 3.10 and Snakemake
```
conda install -n base -c conda-forge mamba
conda activate base
mamba create -c conda-forge -c bioconda -n acmga python=3.10 snakemake
```
### 2、Install docker and singularity following the documentation instructions

 - [singularity installation guide](https://github.com/sylabs/singularity/blob/master/INSTALL.md)
 - [docker installation guide](https://docs.docker.com/engine/install/ubuntu/)

## Running ACMGA
To test the pipeline before running on your own data, you can align some Arabidopsis sequences. 
### 1、Download the code and activat the environment
```
git clone https://github.com/HFzzzzzzz/ACMGA.git
conda activate acmga
```

### 2、Download Arabidopsis data
```
cd ACMGA/data
sh download.sh
cd ..
```
### 3、Generate `command.sh`
```
snakemake  -j 5 --configfile config/config.yaml   --use-singularity  --singularity-args "-B  $(pwd)"
```
### 4、Run `command.sh`
```
docker login
docker run -v $(pwd):/data --rm -it mgatools/acmga:1.0
sh command.sh
```
 # <a name="section1">Generating a Guide Tree</a>
 ## 1、Use the [GEAN](https://github.com/baoxingsong/GEAN) tool to generate protein sequences by inputting FASTA, GFF files, and corresponding CDS and gene sequence files
```
./gean gff2seq -i /media/zhf/ext1/Downloads/gff/gff_chr1-5/Cvi.protein-coding.genes.v2.5.2019-10-09.gff3 -r /media/zhf/ext1/Downloads/fasta/fasta_chr1-5/Cvi.chr.all.v2.0.fasta -p /media/zhf/ext1/Downloads/protein/Cvi.protein.fa -c /media/zhf/ext1/Downloads/protein/Cvi.cds.fa -g /media/zhf/ext1/Downloads/protein/Cvi.gene.fa

./gean gff2seq -i /media/zhf/ext1/Downloads/gff/gff_chr1-5/An-1.protein-coding.genes.v2.5.2019-10-09.gff3 -r /media/zhf/ext1/Downloads/fasta/fasta_chr1-5/An-1.chr.all.v2.0.fasta -p /media/zhf/ext1/Downloads/protein/An-1.protein.fa -c /media/zhf/ext1/Downloads/protein/An-1.cds.fa -g /media/zhf/ext1/Downloads/protein/An-1.gene.fa

./gean gff2seq -i /media/zhf/ext1/Downloads/gff/gff_chr1-5/Ler.protein-coding.genes.v2.5.2019-10-09.gff3 -r /media/zhf/ext1/Downloads/fasta/fasta_chr1-5/Ler.chr.all.v2.0.fasta -p /media/zhf/ext1/Downloads/protein/Ler.protein.fa -c /media/zhf/ext1/Downloads/protein/Ler.cds.fa -g /media/zhf/ext1/Downloads/protein/Ler.gene.fa

./gean gff2seq -i /media/zhf/ext1/Downloads/gff/gff_chr1-5/Arabidopsis_thaliana.TAIR10.56.gff3 -r /media/zhf/ext1/Downloads/fasta/fasta_chr1-5/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -p /media/zhf/ext1/Downloads/protein/Arabidopsis_thaliana.protein.fa -c /media/zhf/ext1/Downloads/protein/Arabidopsis_thaliana.cds.fa -g /media/zhf/ext1/Downloads/protein/Arabidopsis_thaliana.gene.fa
```
## 2、Use the [OrthoFinder](https://github.com/davidemms/OrthoFinder) tool to generate a guide tree
 Create the ExampleData folder in the OrthoFinder directory. Put the generated protein sequence into the ExampleData folder. Use the following command to generate a guide tree.
```
OrthoFinder/orthofinder -f OrthoFinder/ExampleData
```
`OrthoFinder/ExampleData/OrthoFinder/Results_xxx/Species_Tree/SpeciesTree_rooted_node_labels.txt` is the generated guide tree file

# Explanation of Output files
Intermediate results of ACMGA are written to the `data` directory or subdirectories with outputs from different steps of the pipeline. The final output, the multiple genome alignment result, is saved as `result/evolverPlants.hal`.

# HAL Calling Vairants 
Refer to this [document](https://github.com/HFzzzzzzz/ACMGA/blob/master/result/README.md) for calling SNPs and INDELs in the HAL file. Additionally, you have the flexibility to perform further analysis using alternative methods. The Cactus software provides support for various extensions.

# Troubleshooting

## Common errors
#### Input file errors
When the snakemake run terminates with an error despite snakemake (version > 6.0.0) being correctly installed, there are several common causes related to input files:

-   Input FASTA files and GFF files in the data/ directory do not matching samples listed in the config file parameters  `species`.
-   Input FASTA files and GFF files having chromosomes/scaffolds with special characters; ideally, use names consisting of alphanumeric characters only,such as `chr1`.
-   The config.yaml ancestor parameters not being sufficient. Set the number to your ancestor nodesl; if unsure of the exact count of ancestor nodes, then set as many as possible within the maximum range (N0-N(2^(k-1)-1), where
    is the depth of the guide tree), to avoid the error of insufficient ancestor nodes.
- If the test case fails, please check for incomplete data downloads due to network problems.


# How to cite

If you use ACMGA in your work, please cite:
