
# ACMGA

ACMGA is a reference-free Multiple-genome alignment pipeline. A simplified schema of the pipeline is shown below.
![ACMGA](https://github.com/HFzzzzzzz/ACMGA/raw/master/workflow/image/schematic.jpg)





# Build environment
ACMGA requires Python = 3.10 along with Biopython libraries. If you don't have biopython, you need to install it
```
pip install bioython
```
The use of ACMGA currently requires the support of snakemake (>6.0.0), docker and singularity. Please ensure that they have been installed before running ACMGA.

# Quickstart
For a quickstart with your own data, you can follow the instructions below. We recommend testing the pipeline with our test data first (see section  **Testing the pipeline**), to ensure the pipeline will work correctly.

To get started, clone this repository.

You can now prepare the run with the pipeline by doing the following:
1.  Place your FASTA sequences and gff file suffixed with  `.gff` and a  guide tree  into  `ACMGA/data`
2.  Place the CDS sequences set from all the input genomes into `ACMGA/data`
3.  Copy  `ACMGA/config/config.yaml`  to, for example,  `ACMGA/config/myconfig.yaml`  and edit the new config file to include your FASTA sequences (parameter  `fasta:`) and gff  (parameter  `gff:`). The Fasta and gff  must match filenames in  `ACMGA/data`.



The pipeline can then be executed from the  `/ACMGA`  directory in two steps as shown below. 
```
cd ACMGA
snakemake  -j 5 --configfile config/config.yaml   --use-singularity  --singularity-args "-B $(pwd)/data"
```
First step will generate the result.sh script in the data file.


The second step is to enter the docker environment and run the script

```
cd ACMGA/data
docker run -v $(pwd):/data --rm -it mgatools/acmga:1.0
sh command.sh
```

Use docker image from  [the latest release](https://hub.docker.com/repository/docker/mgatools/acmga/general) 
# Testing the pipeline

Assume that you already have a conda environment named **testPipeline** with python 3.10. Snakemake has been successfully installed in this environment. Docker and Singularity have been installed successfully.

To test the pipeline before running on your own data, you can align some Arabidopsis sequences. 
```
git clone https://github.com/HFzzzzzzz/ACMGA.git
conda activate testPipeline
#1、Prepare data
cd ACMGA/data
sh download.sh
cd ..
#2、Generate command.sh script
snakemake  -j 5 --configfile config/config.yaml   --use-singularity  --singularity-args "-B  $(pwd)/data" 
#3、Run command.sh
docker login
cd data
docker run -v $(pwd):/data --rm -it mgatools/acmga:1.0
sh command.sh
```
# Software requirements
Snakemake is required to run this pipeline and we recommend snakemake version 6.0.0 or higher. The recommended installation is shown below. For more details see  [snakemake installation guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

```
conda install -n base -c conda-forge mamba
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

Note that the installation should use the exact commands above, including the exact channel priority, otherwise snakemake may be improperly installed.

To run the pipeline using a prebuilt singularity container, you must have singularity installed on your system. Singularity is installed on many HPC systems and can be used without root privileges. However, note that installation of singularity does require root privileges. If you want to install singularity and have these privileges, you can find up-to-date instructions on how to do so in the official  [singularity installation guide](https://github.com/sylabs/singularity/blob/master/INSTALL.md).
# Explanation of output files
The intermediate results of ACMGA are written to the data directory or subdirectories with outputs from different steps of the pipeline. The final output that multiple genome alignment result  as `result/evolverPlants.hal`

# Troubleshooting

## Common errors
#### Input file errors
When the snakemake run terminates with an error despite snakemake (version > 6.0.0) being correctly installed, there are several common causes related to input files:

-   Input FASTA files and gff files in the /data directory do not match samples listed in the config file parameters  `species`.
-   Input FASTA files and gff files have chromosomes/scaffolds with special characters; ideally, we want names consisting of only alphanumeric characters.
-   The config.yaml ancestor parameters is not enough. It should set the number according to your ancestor nodes, If you are unsure of the exact count of ancestor nodes, then set as many as possible. The maximum range is N0-N(2^(k-1)^-1) (k is depth of the guide tree), to avoid the error of insufficient ancestor nodes.
- If the test case is fails, please check whether the data download is incomplete due to network problems.


# How to cite

If you use ACMGA in your work, please cite:
