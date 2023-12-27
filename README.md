 

# ACMGA

ACMGA is a reference-free Multiple-genome alignment pipeline. A simplified schema of the pipeline is shown below.
![ACMGA](https://github.com/HFzzzzzzz/ACMGA/raw/master/workflow/image/schematic.jpg)





# Building environment
ACMGA requires Python = 3.10 along with Biopython libraries. If you don't have biopython, you need to install it
```
pip install biopython
```
## ACMGA supports building a local environment and using docker image
### Building a local environment
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

Using this approach, you need to make slight modifications to some of the paths within `command.sh`
### Using docker image
The use of ACMGA currently requires the support of snakemake (>6.0.0), docker and singularity. Please ensure that they have been installed before running ACMGA. We recommend that you use this approach.



# Quickstart
For a quickstart with your own data, you can follow the instructions below. We recommend testing the pipeline with our test data first (see section [**Testing the pipeline**](#section2)), to ensure the pipeline will work correctly.

To get started, clone this repository.

You can now prepare the run with the pipeline by doing the following:
 1.  Place your FASTA sequences, gff files suffixed with  `.gff3` and a  guide tree  into  `ACMGA/data`
 2.  Place the CDS sequences set from all the input genomes into `ACMGA/data`
 3. For example
 
	 3.1 Copy  `ACMGA/config/config.yaml` to `ACMGA/config/myconfig.yaml` 
	 
	 3.2 Edit the `ACMGA/config/myconfig.yaml` to include :
	 
	-  Input FASTA sequences name ( parameter  `fasta:` )
	-  Input GFF file name and ancestral GFF file name ( parameter  `gff3:` )
	-  Path for the collection of CDS ( parameter `nonDuplicateCDS:` ), you can use this [script](https://github.com/HFzzzzzzz/ACMGA/blob/master/workflow/scripts/CombineCDS.py) to merge CDS files and obtain the `non_duplicate_CDS.fa`
	-  Path of the FASTA and the GFF files (parameter  `path:` )
	-  Species name ( parameter  `species:` )
	-  Ancestor name ( parameter `ancestor:` )
	-  Path of guide tree ( parameter `Tree:` ), you can use these [steps](#section1) to generate a guide tree



The pipeline can then be executed from the  `ACMGA/`  directory in two steps as shown below. 
```
cd ACMGA
snakemake  -j 5 --configfile config/myconfig.yaml   --use-singularity  --singularity-args "-B $(pwd)"
```
First step will generate the `command.sh` script in the data file.


Second step is to enter the docker environment and run `command.sh`

```
docker login
docker run -v $(pwd):/data --rm -it mgatools/acmga:1.0
sh command.sh
```

# <a name="section2">Testing the pipeline</a>
## Building Environment using docker

### 1、Creating a conda environment named "acmga" with Python 3.10 and Snakemake
```
conda install -n base -c conda-forge mamba
conda activate base
mamba create -c conda-forge -c bioconda -n acmga python=3.10 snakemake
```
### 2、Installing docker and singularity following the documentation instructions

 - [singularity installation guide](https://github.com/sylabs/singularity/blob/master/INSTALL.md)
 - [docker installation guide](https://docs.docker.com/engine/install/ubuntu/)

## Running ACMGA
To test the pipeline before running on your own data, you can align some Arabidopsis sequences. 
### 1、Downloading the code and activating the environment
```
git clone https://github.com/HFzzzzzzz/ACMGA.git
conda activate acmga
```

### 2、Downloading Arabidopsis data
```
cd ACMGA/data
sh download.sh
cd ..
```
### 3、Generating `command.sh`
```
snakemake  -j 5 --configfile config/config.yaml   --use-singularity  --singularity-args "-B  $(pwd)"`
```
### 4、Running `command.sh`
```
docker login
docker run -v $(pwd):/data --rm -it mgatools/acmga:1.0
sh command.sh
```
 # <a name="section1">Generating a guide tree</a>
 ## 1、Using the [GEAN](https://github.com/baoxingsong/GEAN) tool to generate protein sequences by inputting FASTA, GFF files and corresponding CDS, gene sequence files
```
./gean gff2seq -i /media/zhf/ext1/Downloads/gff/gff_chr1-5/Cvi.protein-coding.genes.v2.5.2019-10-09.gff3 -r /media/zhf/ext1/Downloads/fasta/fasta_chr1-5/Cvi.chr.all.v2.0.fasta -p /media/zhf/ext1/Downloads/protein/Cvi.protein.fa -c /media/zhf/ext1/Downloads/protein/Cvi.cds.fa -g /media/zhf/ext1/Downloads/protein/Cvi.gene.fa

./gean gff2seq -i /media/zhf/ext1/Downloads/gff/gff_chr1-5/An-1.protein-coding.genes.v2.5.2019-10-09.gff3 -r /media/zhf/ext1/Downloads/fasta/fasta_chr1-5/An-1.chr.all.v2.0.fasta -p /media/zhf/ext1/Downloads/protein/An-1.protein.fa -c /media/zhf/ext1/Downloads/protein/An-1.cds.fa -g /media/zhf/ext1/Downloads/protein/An-1.gene.fa

./gean gff2seq -i /media/zhf/ext1/Downloads/gff/gff_chr1-5/Ler.protein-coding.genes.v2.5.2019-10-09.gff3 -r /media/zhf/ext1/Downloads/fasta/fasta_chr1-5/Ler.chr.all.v2.0.fasta -p /media/zhf/ext1/Downloads/protein/Ler.protein.fa -c /media/zhf/ext1/Downloads/protein/Ler.cds.fa -g /media/zhf/ext1/Downloads/protein/Ler.gene.fa

./gean gff2seq -i /media/zhf/ext1/Downloads/gff/gff_chr1-5/Arabidopsis_thaliana.TAIR10.56.gff3 -r /media/zhf/ext1/Downloads/fasta/fasta_chr1-5/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -p /media/zhf/ext1/Downloads/protein/Arabidopsis_thaliana.protein.fa -c /media/zhf/ext1/Downloads/protein/Arabidopsis_thaliana.cds.fa -g /media/zhf/ext1/Downloads/protein/Arabidopsis_thaliana.gene.fa
```
## 2、Using the [OrthoFinder](https://github.com/davidemms/OrthoFinder) tool to generate guide tree
 Create the ExampleData folder in the OrthoFinder directory, put the generated protein sequence into the ExampleData folder, and use the following command to generate a guide tree
```
OrthoFinder/orthofinder -f OrthoFinder/ExampleData
```
`OrthoFinder/ExampleData/OrthoFinder/Results_xxx/Species_Tree/SpeciesTree_rooted_node_labels.txt` is the generated guide tree file

# Explanation of output files
The intermediate results of ACMGA are written to the `data` directory or subdirectories with outputs from different steps of the pipeline. The final output that multiple genome alignment result  as `result/evolverPlants.hal`

# HAL calling vairants 
You can refer to this [document](https://github.com/HFzzzzzzz/ACMGA/blob/master/result/README.md) for calling SNPs and INDELs in the HAL file. Additionally, you have the flexibility to perform further analysis using alternative methods. The Cactus software provides support for various extensions.

# Troubleshooting

## Common errors
#### Input file errors
When the snakemake run terminates with an error despite snakemake (version > 6.0.0) being correctly installed, there are several common causes related to input files:

-   Input FASTA files and GFF files in the data/ directory do not match samples listed in the config file parameters  `species`.
-   Input FASTA files and GFF files have chromosomes/scaffolds with special characters; ideally, we want names consisting of only alphanumeric characters,such as chr1.
-   The config.yaml ancestor parameters is not enough. It should set the number according to your ancestor nodes, If you are unsure of the exact count of ancestor nodes, then set as many as possible. The maximum range is N0-N(2^(k-1)-1)
 (k is depth of the guide tree), to avoid the error of insufficient ancestor nodes.
- If the test case is fails, please check whether the data download is incomplete due to network problems.


# How to cite

If you use ACMGA in your work, please cite:
