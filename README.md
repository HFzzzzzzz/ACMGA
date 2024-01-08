

 

# ACMGA

AnchorWave-Cactus Multiple Genome Alignment (ACMGA) is a reference-free multiple-genome alignment pipeline. It leverages the power of AnchorWave, a pairwise genome alignment tool that utilizes collinearity and global alignment algorithms. This enables ACMGA to effectively align repetitive sequence regions and accurately identify long INDELs (>50bp). Moreover, ACMGA incorporates the Progressive Cactus algorithm to generate ancestor sequences and implement progressive strategies. This combination of techniques makes ACMGA particularly well-suited for aligning plant genomes that are enriched with repetitive sequences. The simplified schema of the pipeline is depicted below.
![ACMGA](https://github.com/HFzzzzzzz/ACMGA/raw/master/workflow/image/schematic.jpg)


## Table of Contents
- [Building Environment](#section1)
- [Testing the pipeline](#section2)
- [Quickstart](#section3)
- [Please note](#section4)
- [Generat a guide tree](#section5)
- [Explanation of Output files](#section6)
- [HAL Calling Vairants](#section7)
- [Troubleshoot](#section8)
- [How to cite](#section9)


## <a name="section1">Building Environment</a>
### ACMGA supports building the environment using the Docker image or locally.

#### Using the Docker image:
ACMGA currently relies on Snakemake (>6.0.0), Docker, and Singularity. Please make sure these dependencies are installed before running ACMGA. We recommend using this approach.

#### Building the local environment:
- Python3.10
- Biopython
- [Snakemake(>6.0)](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- [AnchorWave](https://github.com/baoxingsong/AnchorWave)
- [Cactus](https://github.com/ComparativeGenomicsToolkit/cactus)
- [SAMtools](http://www.htslib.org/)
- [Minimap2](https://github.com/lh3/minimap2)
- [bedtools](https://github.com/arq5x/bedtools2)
- [bedToGenePred](https://github.com/ENCODE-DCC/kentUtils/tree/master/src/hg/bedToGenePred)
- [genePredToGtf](https://github.com/ENCODE-DCC/kentUtils/tree/master/src/hg/genePredToGtf)
- [Gffread](https://github.com/gpertea/gffread)
- [K8](https://github.com/attractivechaos/k8)
- [maf-convert](https://gitlab.com/mcfrith/last/-/blob/main/bin/maf-convert)

Using this approach, slight modifications to some of the paths (the path of the tools used in the pipeline and the path of the data used should be the same as the local) within `command.sh` are necessary, give the software executable permissions, and add the path of the executable program to the `PATH`.

## <a name="section2">Testing the pipeline</a>
### Building Environment using Docker

#### 1、Create a conda environment named "acmga" with Python 3.10 and Snakemake(>6.0)
```
conda install -n base -c conda-forge mamba
conda activate base
mamba create -c conda-forge -c bioconda -n acmga python=3.10 snakemake
```
#### 2、Install Docker and Singularity following the documentation instructions

 - [Singularity installation guide](https://github.com/sylabs/singularity/blob/master/INSTALL.md)
 - [Docker installation guide](https://docs.docker.com/engine/install/ubuntu/)

### Running ACMGA
To test the pipeline, you can align some Arabidopsis genomes. 
#### 1、Download the code, activate the environment and install biopython 
```
git clone https://github.com/HFzzzzzzz/ACMGA.git
conda activate acmga
conda install -c conda-forge biopython
```

#### 2、Download Arabidopsis data
```
cd ACMGA/data
sh download.sh
cd ..
```
#### 3、Generate `command.sh`
```
snakemake  -j 5 --configfile config/config.yaml   --use-singularity  --singularity-args "-B  $(pwd)"
```
#### 4、Run `command.sh`
```
docker login
docker run -v $(pwd):/data --rm -it mgatools/acmga:1.0
sh command.sh
```

## <a name="section3">Quickstart </a>
For a quickstart with your own data, you can follow the instructions below. We recommend testing the pipeline with our test data first to ensure the pipeline will work correctly.

After testing the pipeline, the environment has been build successfully. Now you just need to prepare your own data and modify the configuration `myconfig.yaml` file to run multiple genome alignment on your own data.

You can now prepare the run with the pipeline by doing the following:
 1.  Placing your FASTA sequences, Gff files (suffixed with  `.gff3`), and a  guide tree  into  `ACMGA/data/`.
 2.  Placing the CDS sequences set from all the input genomes into `ACMGA/data/`.
 3.  For example
 
	 3.1 Copying  `ACMGA/config/config.yaml` to `ACMGA/config/myconfig.yaml` 
	 
	 3.2 Editing the `ACMGA/config/myconfig.yaml` to include :
	 
		-  Input FASTA sequences name (parameter  `fasta:` ).
		-  Input GFF files name and ancestral GFF files name (parameter  `gff3:` ).
		-  Path for the collection of CDS (parameter `nonDuplicateCDS:` ), using this [script](https://github.com/HFzzzzzzz/ACMGA/blob/master/workflow/scripts/CombineCDS.py) to merge CDS files and obtain `non_duplicate_CDS.fa`.
		-  Path of the FASTA and the GFF files (parameter  `path:` ).
		-  Species name (parameter  `species:` ).
		-  The name of the ancestor sequence (parameter `ancestor:` ).
		-  Path of guide tree (parameter `Tree:` ), generated using recommended [steps](#section1).



The pipeline can then be executed from the  `ACMGA/`  directory in two steps.

1.The first step generates the `command.sh` script in the `ACMGA/` 
```
cd ACMGA
snakemake  -j 5 --configfile config/myconfig.yaml   --use-singularity  --singularity-args "-B $(pwd)"
```

2.The second step is to enter the docker environment and run `command.sh`

```
docker login
docker run -v $(pwd):/data --rm -it mgatools/acmga:1.0
sh command.sh
```
## <a name="section4">Please note</a>
- The suffix of the Gff file of the input genome and the generated ancestral genome is `.gff3` in parameter `gff:`. There is no requirement for the DNA sequence file name(`.fasta` or `.fa`).
- The species name (parameter  `species:`) should be a substring of the fasta and gff names.
- The names of the ancestor genome comes from the tree structure. (such as `(Ler:0.00493038,(Cvi:0.0145906,(Arabidopsis_thaliana:0.0117518,An-1:0.00833626)N2:0.0070939)N1:0.00493038)N0;`.The name of the ancestor genome are `N2`,`N1`,`N0`).
- The number of threads for AnchorWave and Cactus (parameter `proaliParamters`, `genoaliParamters`, `cactus_threads`)should be set according to your memory configuration to prevent memory overflow.
- The chromosome names in fasta and gff of the species need to be consistent. Special characters such as spaces cannot be carried. The best example is `chr{chromosome number}`.
- Sequences that are not at the chromosome level can be removed, such as `scaf`.


 ## <a name="section5">Generate a guide tree</a>
 ### 1、Use the [GEAN](https://github.com/baoxingsong/GEAN) tool to generate protein sequences by inputting FASTA, GFF files, and corresponding CDS and gene sequence files
```
./gean gff2seq -i /media/zhf/ext1/Downloads/gff/gff_chr1-5/Cvi.protein-coding.genes.v2.5.2019-10-09.gff3 -r /media/zhf/ext1/Downloads/fasta/fasta_chr1-5/Cvi.chr.all.v2.0.fasta -p /media/zhf/ext1/Downloads/protein/Cvi.protein.fa -c /media/zhf/ext1/Downloads/protein/Cvi.cds.fa -g /media/zhf/ext1/Downloads/protein/Cvi.gene.fa

./gean gff2seq -i /media/zhf/ext1/Downloads/gff/gff_chr1-5/An-1.protein-coding.genes.v2.5.2019-10-09.gff3 -r /media/zhf/ext1/Downloads/fasta/fasta_chr1-5/An-1.chr.all.v2.0.fasta -p /media/zhf/ext1/Downloads/protein/An-1.protein.fa -c /media/zhf/ext1/Downloads/protein/An-1.cds.fa -g /media/zhf/ext1/Downloads/protein/An-1.gene.fa

./gean gff2seq -i /media/zhf/ext1/Downloads/gff/gff_chr1-5/Ler.protein-coding.genes.v2.5.2019-10-09.gff3 -r /media/zhf/ext1/Downloads/fasta/fasta_chr1-5/Ler.chr.all.v2.0.fasta -p /media/zhf/ext1/Downloads/protein/Ler.protein.fa -c /media/zhf/ext1/Downloads/protein/Ler.cds.fa -g /media/zhf/ext1/Downloads/protein/Ler.gene.fa

./gean gff2seq -i /media/zhf/ext1/Downloads/gff/gff_chr1-5/Arabidopsis_thaliana.TAIR10.56.gff3 -r /media/zhf/ext1/Downloads/fasta/fasta_chr1-5/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -p /media/zhf/ext1/Downloads/protein/Arabidopsis_thaliana.protein.fa -c /media/zhf/ext1/Downloads/protein/Arabidopsis_thaliana.cds.fa -g /media/zhf/ext1/Downloads/protein/Arabidopsis_thaliana.gene.fa
```
### 2、Use the [OrthoFinder](https://github.com/davidemms/OrthoFinder) tool to generate a guide tree
 Create the ExampleData folder in the OrthoFinder directory. Put the generated protein sequence into the ExampleData folder. Use the following command to generate a guide tree.
```
OrthoFinder/orthofinder -f OrthoFinder/ExampleData
```
`OrthoFinder/ExampleData/OrthoFinder/Results_xxx/Species_Tree/SpeciesTree_rooted_node_labels.txt` is the generated guide tree file.

## <a name="section6">Explanation of Output files</a>
Intermediate results of ACMGA are written to the `data` directory or subdirectories with outputs from different steps of the pipeline. The final output, the multiple genome alignment result, is saved as `result/evolverPlants.hal`.

## <a name="section7">HAL Calling Vairants </a>
Refer to this [document](https://github.com/HFzzzzzzz/ACMGA/blob/master/result/README.md) for calling SNPs and INDELs in the HAL file. Additionally, you have the flexibility to perform further analysis using alternative methods. The Cactus toolkit provides support for various extensions.

## <a name="section8">Troubleshoot </a>

### Common errors
When the Snakemake run terminates with an error despite Snakemake (version > 6.0.0) being correctly installed, there are several common causes related to input files:

-   Input FASTA files and GFF files in the `data/` directory do not matching samples listed in the config file parameters  `species`.
-   Input FASTA files and GFF files having chromosomes/scaffolds with special characters; ideally, use names consisting of alphanumeric characters only, such as `chr1`.
-   The config.yaml ancestor parameters not being right. Set ancestor nodes according to your tree file. Ancestor nodes include all non-leaf nodes in the tree.
- If the test case fails, please check for incomplete data downloads due to network problems.


## <a name="section9">How to cite</a>

If you use ACMGA in your work, please cite:
