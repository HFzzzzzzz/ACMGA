

# Result

The final results of ACMGA generate multiple-genome alignment results in HAL format, which are placed in the 
`ACMGA/result` directory.You can proceed with further analysis according to your specific needs.![ACMGA](https://github.com/HFzzzzzzz/ACMGA/raw/master/workflow/image/schematic.jpg)

# HAL file to VCF file
## HAL file to MAF file
The file in HAL format is converted into a chain file for pairwise genome alignments through [the document](https://github.com/ComparativeGenomicsToolkit/hal/blob/chaining-doc/doc/chaining-mapping.md). We convert it into a maf file through [chain2paf](https://github.com/AndreaGuarracino/chain2paf) and [paf2maf](https://github.com/wjwei-handsome/wgatools).
```
hal2fasta /home/zhf/bs674/14t2/huafeng/threeData/SlurmArabidopsis2/ancestor.hal An-1.protein       | faToTwoBit stdin An-1.protein.AnchorWave.2bit
hal2fasta  /home/zhf/bs674/14t2/huafeng/threeData/SlurmArabidopsis2/ancestor.hal Arabidopsis_thaliana.protein  | faToTwoBit stdin Arabidopsis_thaliana.protein.AnchorWave.2bit
halStats --bedSequences An-1.protein /home/zhf/bs674/14t2/huafeng/threeData/SlurmArabidopsis2/ancestor.hal > An-1.protein.AnchorWave.bed
halLiftover --outPSL /home/zhf/bs674/14t2/huafeng/threeData/SlurmArabidopsis2/ancestor.hal An-1.protein An-1.protein.AnchorWave.bed  Arabidopsis_thaliana.protein  /dev/stdout | pslPosTarget stdin An-1.protein_Arab.AnchorWave.psl
axtChain -psl -linearGap=loose An-1.protein_Arab.AnchorWave.psl Arabidopsis_thaliana.protein.AnchorWave.2bit  An-1.protein.AnchorWave.2bit     An-1.protein-Arab.AnchorWave.chain
chain2paf -i An-1.protein-Arab.AnchorWave.chain > /media/zhf/ext1/Arabidopsis/AnchorWave/paf/An-1.protein_Arab.AnchorWave.paf
/media/zhf/ext1/software/wgatools/target/release/wgatools paf2maf   --target /media/zhf/ext1/Arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa --query /media/zhf/ext1/gene2_copy/SlurmArabidopsis/An-1.protein.chr.all.v2.0.fasta /media/zhf/ext1/Arabidopsis/AnchorWave/paf/An-1.protein_Arab.AnchorWave.paf  -o /media/zhf/ext1/Arabidopsis/AnchorWave/maf/An-1.protein_Arab.AnchorWave.maf   -r

hal2fasta /home/zhf/bs674/14t2/huafeng/threeData/SlurmArabidopsis2/ancestor.hal Ler.protein       | faToTwoBit stdin Ler.protein.AnchorWave.2bit
hal2fasta  /home/zhf/bs674/14t2/huafeng/threeData/SlurmArabidopsis2/ancestor.hal Arabidopsis_thaliana.protein  | faToTwoBit stdin Arabidopsis_thaliana.protein.AnchorWave.2bit
halStats --bedSequences Ler.protein /home/zhf/bs674/14t2/huafeng/threeData/SlurmArabidopsis2/ancestor.hal > Ler.protein.AnchorWave.bed
halLiftover --outPSL /home/zhf/bs674/14t2/huafeng/threeData/SlurmArabidopsis2/ancestor.hal Ler.protein Ler.protein.AnchorWave.bed  Arabidopsis_thaliana.protein  /dev/stdout | pslPosTarget stdin Ler.protein_Arab.AnchorWave.psl
axtChain -psl -linearGap=loose Ler.protein_Arab.AnchorWave.psl Arabidopsis_thaliana.protein.AnchorWave.2bit  Ler.protein.AnchorWave.2bit     Ler.protein-Arab.AnchorWave.chain
chain2paf -i Ler.protein-Arab.AnchorWave.chain > /media/zhf/ext1/Arabidopsis/AnchorWave/paf/Ler.protein_Arab.AnchorWave.paf
/media/zhf/ext1/software/wgatools/target/release/wgatools paf2maf   --target /media/zhf/ext1/Arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa --query /media/zhf/ext1/gene2_copy/SlurmArabidopsis/Ler.protein.chr.all.v2.0.fasta /media/zhf/ext1/Arabidopsis/AnchorWave/paf/Ler.protein_Arab.AnchorWave.paf  -o /media/zhf/ext1/Arabidopsis/AnchorWave/maf/Ler.protein_Arab.AnchorWave.maf   -r

hal2fasta /home/zhf/bs674/14t2/huafeng/threeData/SlurmArabidopsis2/ancestor.hal Cvi.protein       | faToTwoBit stdin Cvi.protein.AnchorWave.2bit
hal2fasta  /home/zhf/bs674/14t2/huafeng/threeData/SlurmArabidopsis2/ancestor.hal Arabidopsis_thaliana.protein  | faToTwoBit stdin Arabidopsis_thaliana.protein.AnchorWave.2bit
halStats --bedSequences Cvi.protein /home/zhf/bs674/14t2/huafeng/threeData/SlurmArabidopsis2/ancestor.hal > Cvi.protein.AnchorWave.bed
halLiftover --outPSL /home/zhf/bs674/14t2/huafeng/threeData/SlurmArabidopsis2/ancestor.hal Cvi.protein Cvi.protein.AnchorWave.bed  Arabidopsis_thaliana.protein  /dev/stdout | pslPosTarget stdin Cvi.protein_Arab.AnchorWave.psl
axtChain -psl -linearGap=loose Cvi.protein_Arab.AnchorWave.psl Arabidopsis_thaliana.protein.AnchorWave.2bit  Cvi.protein.AnchorWave.2bit     Cvi.protein-Arab.AnchorWave.chain
chain2paf -i Cvi.protein-Arab.AnchorWave.chain > /media/zhf/ext1/Arabidopsis/AnchorWave/paf/Cvi.protein_Arab.AnchorWave.paf
/media/zhf/ext1/software/wgatools/target/release/wgatools paf2maf   --target /media/zhf/ext1/Arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa --query /media/zhf/ext1/gene2_copy/SlurmArabidopsis/Cvi.protein.chr.all.v2.0.fasta /media/zhf/ext1/Arabidopsis/AnchorWave/paf/Cvi.protein_Arab.AnchorWave.paf  -o /media/zhf/ext1/Arabidopsis/AnchorWave/maf/Cvi.protein_Arab.AnchorWave.maf   -r
```

## MAF file to VCF file
MAF files can be converted into GVCF files through [TASSEL](https://www.maizegenetics.net/tassel), and then the VCF file can be obtained through [GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4).
 ```
awk '{$7=toupper($7)} 1' /home/zhouhf/my_data/Arabidopsis/AnchorWave/maf/An-1.protein_Arab.AnchorWave21.maf  > /home/zhouhf/my_data/Arabidopsis/AnchorWave/maf/An-1.protein_Arab.ANCHORWAVE.maf
awk '{$7=toupper($7)} 1' /home/zhouhf/my_data/Arabidopsis/AnchorWave/maf/Ler.protein_Arab.AnchorWave21.maf  > /home/zhouhf/my_data/Arabidopsis/AnchorWave/maf/Ler.protein_Arab.ANCHORWAVE.maf
awk '{$7=toupper($7)} 1' /home/zhouhf/my_data/Arabidopsis/AnchorWave/maf/Cvi.protein_Arab.AnchorWave21.maf  > /home/zhouhf/my_data/Arabidopsis/AnchorWave/maf/Cvi.protein_Arab.ANCHORWAVE.maf

sed  's/[DMWKSYR]/N/g' /home/zhouhf/my_data/Arabidopsis/AnchorWave/maf/An-1.protein_Arab.ANCHORWAVE.maf > /home/zhouhf/my_data/Arabidopsis/AnchorWave/maf/An-1.Arab.ANCHORWAVE.maf
sed  's/[DMWKSYR]/N/g' /home/zhouhf/my_data/Arabidopsis/AnchorWave/maf/Ler.protein_Arab.ANCHORWAVE.maf > /home/zhouhf/my_data/Arabidopsis/AnchorWave/maf/Ler.Arab.ANCHORWAVE.maf
sed  's/[DMWKSYR]/N/g' /home/zhouhf/my_data/Arabidopsis/AnchorWave/maf/Cvi.protein_Arab.ANCHORWAVE.maf > /home/zhouhf/my_data/Arabidopsis/AnchorWave/maf/Cvi.Arab.ANCHORWAVE.maf

/home/zhouhf/tassel-5-standalone/run_pipeline.pl -Xmx200G -debug -MAFToGVCFPlugin -referenceFasta /home/zhouhf/SlurmArabidopsis/Arabidopsis_thaliana.TAIR10.dna.toplevel.fasta  -mafFile /home/zhouhf/my_data/Arabidopsis/AnchorWave/maf/An-1.Arab.ANCHORWAVE.maf -sampleName An-1_anchorwave -gvcfOutput /home/zhouhf/my_data/Arabidopsis/AnchorWave/gvcf/An-1_ANCHORWAVE.gvcf -fillGaps false -bgzipAndIndex false >/home/zhouhf/my_data/Arabidopsis/AnchorWave/gvcflog/An-1_ANCHORWAVE.log
/home/zhouhf/tassel-5-standalone/run_pipeline.pl -Xmx200G -debug -MAFToGVCFPlugin -referenceFasta /home/zhouhf/SlurmArabidopsis/Arabidopsis_thaliana.TAIR10.dna.toplevel.fasta  -mafFile /home/zhouhf/my_data/Arabidopsis/AnchorWave/maf/Ler.Arab.ANCHORWAVE.maf -sampleName Ler_anchorwave -gvcfOutput /home/zhouhf/my_data/Arabidopsis/AnchorWave/gvcf/Ler_ANCHORWAVE.gvcf -fillGaps false -bgzipAndIndex false >/home/zhouhf/my_data/Arabidopsis/AnchorWave/gvcflog/Ler_ANCHORWAVE.log
/home/zhouhf/tassel-5-standalone/run_pipeline.pl -Xmx200G -debug -MAFToGVCFPlugin -referenceFasta /home/zhouhf/SlurmArabidopsis/Arabidopsis_thaliana.TAIR10.dna.toplevel.fasta  -mafFile /home/zhouhf/my_data/Arabidopsis/AnchorWave/maf/Cvi.Arab.ANCHORWAVE.maf -sampleName Cvi_anchorwave -gvcfOutput /home/zhouhf/my_data/Arabidopsis/AnchorWave/gvcf/Cvi_ANCHORWAVE.gvcf -fillGaps false -bgzipAndIndex false >/home/zhouhf/my_data/Arabidopsis/AnchorWave/gvcflog/Cvi_ANCHORWAVE.log

gatk --java-options "-Xmx200g" GenotypeGVCFs -R /home/zhouhf/SlurmArabidopsis/Arabidopsis_thaliana.TAIR10.dna.toplevel.delchr.fasta  -stand-call-conf 0 -ploidy 1 -V /home/zhouhf/my_data/Arabidopsis/AnchorWave/Filter/An-1_AnchorWave.gvcf -O /home/zhouhf/my_data/Arabidopsis/AnchorWave/vcf/An-1_anchorwave.vcf --cloud-prefetch-buffer 10000 --cloud-index-prefetch-buffer 10000 --genomicsdb-max-alternate-alleles 110 --max-alternate-alleles 100 --gcs-max-retries 1000
gatk --java-options "-Xmx200g" GenotypeGVCFs -R /home/zhouhf/SlurmArabidopsis/Arabidopsis_thaliana.TAIR10.dna.toplevel.delchr.fasta  -stand-call-conf 0 -ploidy 1 -V /home/zhouhf/my_data/Arabidopsis/AnchorWave/Filter/Ler_AnchorWave.gvcf -O /home/zhouhf/my_data/Arabidopsis/AnchorWave/vcf/Ler_anchorwave.vcf --cloud-prefetch-buffer 10000 --cloud-index-prefetch-buffer 10000 --genomicsdb-max-alternate-alleles 110 --max-alternate-alleles 100 --gcs-max-retries 1000
gatk --java-options "-Xmx200g" GenotypeGVCFs -R /home/zhouhf/SlurmArabidopsis/Arabidopsis_thaliana.TAIR10.dna.toplevel.delchr.fasta  -stand-call-conf 0 -ploidy 1 -V /home/zhouhf/my_data/Arabidopsis/AnchorWave/Filter/Cvi_AnchorWave.gvcf -O /home/zhouhf/my_data/Arabidopsis/AnchorWave/vcf/Cvi_anchorwave.vcf --cloud-prefetch-buffer 10000 --cloud-index-prefetch-buffer 10000 --genomicsdb-max-alternate-alleles 110 --max-alternate-alleles 100 --gcs-max-retries 1000

```
## Calling SNPs and INDELs from VCF files
```
gatk --java-options -XX:ParallelGCThreads=120 SelectVariants -V  /home/zhouhf/my_data/Arabidopsis/AnchorWave/vcf/An-1_anchorwave.vcf --select-type-to-include SNP    -O  /home/zhouhf/my_data/Arabidopsis/AnchorWave/SNP/An-1_anchorwave.SNP.vcf
gatk --java-options -XX:ParallelGCThreads=120 SelectVariants -V  /home/zhouhf/my_data/Arabidopsis/AnchorWave/vcf/Ler_anchorwave.vcf --select-type-to-include SNP    -O  /home/zhouhf/my_data/Arabidopsis/AnchorWave/SNP/Ler_anchorwave.SNP.vcf
gatk --java-options -XX:ParallelGCThreads=120 SelectVariants -V  /home/zhouhf/my_data/Arabidopsis/AnchorWave/vcf/Cvi_anchorwave.vcf --select-type-to-include SNP    -O  /home/zhouhf/my_data/Arabidopsis/AnchorWave/SNP/Cvi_anchorwave.SNP.vcf

gatk --java-options -XX:ParallelGCThreads=120 SelectVariants -V  /home/zhouhf/my_data/Arabidopsis/AnchorWave/vcf/An-1_anchorwave.vcf --select-type-to-include INDEL    -O /home/zhouhf/my_data/Arabidopsis/AnchorWave/INDEL/An-1_anchorwave.INDEL.vcf
gatk --java-options -XX:ParallelGCThreads=120 SelectVariants -V  /home/zhouhf/my_data/Arabidopsis/AnchorWave/vcf/Ler_anchorwave.vcf --select-type-to-include INDEL    -O /home/zhouhf/my_data/Arabidopsis/AnchorWave/INDEL/Ler_anchorwave.INDEL.vcf
gatk --java-options -XX:ParallelGCThreads=120 SelectVariants -V  /home/zhouhf/my_data/Arabidopsis/AnchorWave/vcf/Cvi_anchorwave.vcf --select-type-to-include INDEL    -O /home/zhouhf/my_data/Arabidopsis/AnchorWave/INDEL/Cvi_anchorwave.INDEL.vcf
```
