cactus_sanitizeFastaHeaders  /data/Arabidopsis/Ler.chr.all.v2.0.fasta  Ler.protein > /data/Arabidopsis/Ler.protein.sanitized.fa
cactus_sanitizeFastaHeaders  /data/Arabidopsis/Cvi.chr.all.v2.0.fasta  Cvi.protein > /data/Arabidopsis/Cvi.protein.sanitized.fa
cactus_sanitizeFastaHeaders  /data/Arabidopsis/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa  Arabidopsis_thaliana.protein > /data/Arabidopsis/Arabidopsis_thaliana.protein.sanitized.fa
cactus_sanitizeFastaHeaders  /data/Arabidopsis/An-1.chr.all.v2.0.fasta  An-1.protein > /data/Arabidopsis/An-1.protein.sanitized.fa
mkdir /data/Arabidopsis/subtree1
#align Arabidopsis_thaliana.TAIR10.dna.toplevel.fa	against	An-1.chr.all.v2.0.fasta

#m=============================================================================================================================================================================================================================1

anchorwave gff2seq -i /data/Arabidopsis/Arabidopsis_thaliana.TAIR10.56.gff3   -r /data/Arabidopsis/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa  -o Arabidopsis_thaliana_cds.fa
minimap2  -x splice -t 35 -k 12 -a -p 0.4 -N 20 /data/Arabidopsis/An-1.chr.all.v2.0.fasta  Arabidopsis_thaliana_cds.fa > An-1_Arabidopsis_thaliana_cds.sam
minimap2  -x splice -t 35 -k 12 -a -p 0.4 -N 20 /data/Arabidopsis/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa  Arabidopsis_thaliana_cds.fa > Arabidopsis_thaliana_Arabidopsis_thaliana_ref.sam
/usr/bin/time anchorwave genoAli -i /data/Arabidopsis/Arabidopsis_thaliana.TAIR10.56.gff3 -as Arabidopsis_thaliana_cds.fa   -r /data/Arabidopsis/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -a An-1_Arabidopsis_thaliana_cds.sam  -ar Arabidopsis_thaliana_Arabidopsis_thaliana_ref.sam -s /data/Arabidopsis/An-1.chr.all.v2.0.fasta -n /data/Arabidopsis/Arabidopsis_thaliana_An-1.anchors -o /data/Arabidopsis/Arabidopsis_thaliana_An-1.anchorwave.maf -f /data/Arabidopsis/Arabidopsis_thaliana_An-1.anchorwave.f.maf -t 5 -IV  > Arabidopsis_thaliana_An-1.log 2>&1
maf-convert sam /data/Arabidopsis/Arabidopsis_thaliana_An-1.anchorwave.maf > /data/Arabidopsis/Arabidopsis_thaliana_An-1.anchorwave.sam
samtools view -H --reference /data/Arabidopsis/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa  /data/Arabidopsis/Arabidopsis_thaliana_An-1.anchorwave.sam > /data/Arabidopsis/Arabidopsis_thaliana_An-1.2.anchorwave.sam
grep -v "^@" /data/Arabidopsis/Arabidopsis_thaliana_An-1.anchorwave.sam >> /data/Arabidopsis/Arabidopsis_thaliana_An-1.2.anchorwave.sam
k8 /data/paftools.js  sam2paf /data/Arabidopsis/Arabidopsis_thaliana_An-1.2.anchorwave.sam > /data/Arabidopsis/Arabidopsis_thaliana_An-1.tmp.anchorwave.paf
python /data/replace_ref_que.py -i /data/Arabidopsis/Arabidopsis_thaliana_An-1.tmp.anchorwave.paf -o /data/Arabidopsis/Arabidopsis_thaliana_An-1.anchorwave.paf -r Arabidopsis_thaliana.protein -q An-1.protein
cat /data/Arabidopsis/Arabidopsis_thaliana_An-1.anchorwave.paf > /data/Arabidopsis/Arabidopsis_thaliana_An-1.anchorwave_invert.paf
paf_invert -i /data/Arabidopsis/Arabidopsis_thaliana_An-1.anchorwave.paf >> /data/Arabidopsis/Arabidopsis_thaliana_An-1.anchorwave_invert.paf
paf_chain -i /data/Arabidopsis/Arabidopsis_thaliana_An-1.anchorwave_invert.paf --maxGapLength 1000000 --chainGapOpen 5000 --chainGapExtend 1 --trimFraction 0.02 --logLevel DEBUG  >> /data/Arabidopsis/subtree1/anchorwave_chain.paf

#m=============================================================================================================================================================================================================================1

#genome_a   Ler.chr.all.v2.0.fasta   genome_b   Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
anchorwave gff2seq -i /data/Arabidopsis/Ler.protein-coding.genes.v2.5.2019-10-09.gff3   -r /data/Arabidopsis/Ler.chr.all.v2.0.fasta  -o Ler_cds.fa
minimap2  -x splice -t 35 -k 12 -a -p 0.4 -N 20 /data/Arabidopsis/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa  Ler_cds.fa > Arabidopsis_thaliana_Ler_cds.sam
minimap2  -x splice -t 35 -k 12 -a -p 0.4 -N 20 /data/Arabidopsis/Ler.chr.all.v2.0.fasta  Ler_cds.fa > Ler_Ler_ref.sam
/usr/bin/time anchorwave genoAli -i /data/Arabidopsis/Ler.protein-coding.genes.v2.5.2019-10-09.gff3 -as Ler_cds.fa   -r /data/Arabidopsis/Ler.chr.all.v2.0.fasta -a Arabidopsis_thaliana_Ler_cds.sam  -ar Ler_Ler_ref.sam -s /data/Arabidopsis/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -n /data/Arabidopsis/Ler_Arabidopsis_thaliana.anchors -o /data/Arabidopsis/Ler_Arabidopsis_thaliana.anchorwave.maf -f /data/Arabidopsis/Ler_Arabidopsis_thaliana.anchorwave.f.maf -t 5 -IV  > Ler_Arabidopsis_thaliana.avx2.log 2>&1
maf-convert sam /data/Arabidopsis/Ler_Arabidopsis_thaliana.anchorwave.maf  > /data/Arabidopsis/Ler_Arabidopsis_thaliana.2.anchorwave.sam
samtools view -H --reference /data/Arabidopsis/Ler.chr.all.v2.0.fasta  /data/Arabidopsis/Ler_Arabidopsis_thaliana.2.anchorwave.sam > /data/Arabidopsis/Ler_Arabidopsis_thaliana.3.anchorwave.sam
grep -v "^@" /data/Arabidopsis/Ler_Arabidopsis_thaliana.2.anchorwave.sam >> /data/Arabidopsis/Ler_Arabidopsis_thaliana.3.anchorwave.sam
k8 /data/paftools.js   sam2paf /data/Arabidopsis/Ler_Arabidopsis_thaliana.3.anchorwave.sam > /data/Arabidopsis/Ler_Arabidopsis_thaliana.tmp.anchorwave.paf
python /data/replace_ref_que.py -i /data/Arabidopsis/Ler_Arabidopsis_thaliana.tmp.anchorwave.paf -o /data/Arabidopsis/Ler_Arabidopsis_thaliana.anchorwave.paf -r Ler.protein -q Arabidopsis_thaliana.protein

#m=============================================================================================================================================================================================================================1

#genome_a   Cvi.chr.all.v2.0.fasta   genome_b   Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
anchorwave gff2seq -i /data/Arabidopsis/Cvi.protein-coding.genes.v2.5.2019-10-09.gff3   -r /data/Arabidopsis/Cvi.chr.all.v2.0.fasta  -o Cvi_cds.fa
minimap2  -x splice -t 35 -k 12 -a -p 0.4 -N 20 /data/Arabidopsis/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa  Cvi_cds.fa > Arabidopsis_thaliana_Cvi_cds.sam
minimap2  -x splice -t 35 -k 12 -a -p 0.4 -N 20 /data/Arabidopsis/Cvi.chr.all.v2.0.fasta  Cvi_cds.fa > Cvi_Cvi_ref.sam
/usr/bin/time anchorwave genoAli -i /data/Arabidopsis/Cvi.protein-coding.genes.v2.5.2019-10-09.gff3 -as Cvi_cds.fa   -r /data/Arabidopsis/Cvi.chr.all.v2.0.fasta -a Arabidopsis_thaliana_Cvi_cds.sam  -ar Cvi_Cvi_ref.sam -s /data/Arabidopsis/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -n /data/Arabidopsis/Cvi_Arabidopsis_thaliana.anchors -o /data/Arabidopsis/Cvi_Arabidopsis_thaliana.anchorwave.maf -f /data/Arabidopsis/Cvi_Arabidopsis_thaliana.anchorwave.f.maf -t 5 -IV  > Cvi_Arabidopsis_thaliana.avx2.log 2>&1
maf-convert sam /data/Arabidopsis/Cvi_Arabidopsis_thaliana.anchorwave.maf  > /data/Arabidopsis/Cvi_Arabidopsis_thaliana.2.anchorwave.sam
samtools view -H --reference /data/Arabidopsis/Cvi.chr.all.v2.0.fasta  /data/Arabidopsis/Cvi_Arabidopsis_thaliana.2.anchorwave.sam > /data/Arabidopsis/Cvi_Arabidopsis_thaliana.3.anchorwave.sam
grep -v "^@" /data/Arabidopsis/Cvi_Arabidopsis_thaliana.2.anchorwave.sam >> /data/Arabidopsis/Cvi_Arabidopsis_thaliana.3.anchorwave.sam
k8 /data/paftools.js   sam2paf /data/Arabidopsis/Cvi_Arabidopsis_thaliana.3.anchorwave.sam > /data/Arabidopsis/Cvi_Arabidopsis_thaliana.tmp.anchorwave.paf
python /data/replace_ref_que.py -i /data/Arabidopsis/Cvi_Arabidopsis_thaliana.tmp.anchorwave.paf -o /data/Arabidopsis/Cvi_Arabidopsis_thaliana.anchorwave.paf -r Cvi.protein -q Arabidopsis_thaliana.protein
cat /data/Arabidopsis/Cvi_Arabidopsis_thaliana.anchorwave.paf     /data/Arabidopsis/Ler_Arabidopsis_thaliana.anchorwave.paf > /data/Arabidopsis/cat_Ler_Arabidopsis_thaliana_Cvi_Arabidopsis_thaliana.2.anchorwave.paf
cat /data/Arabidopsis/cat_Ler_Arabidopsis_thaliana_Cvi_Arabidopsis_thaliana.2.anchorwave.paf > /data/Arabidopsis/Cvi_Arabidopsis_thaliana.2.anchorwave_invert.paf
paf_invert -i /data/Arabidopsis/cat_Ler_Arabidopsis_thaliana_Cvi_Arabidopsis_thaliana.2.anchorwave.paf >> /data/Arabidopsis/Cvi_Arabidopsis_thaliana.2.anchorwave_invert.paf
paf_chain -i /data/Arabidopsis/Cvi_Arabidopsis_thaliana.2.anchorwave_invert.paf --maxGapLength 1000000 --chainGapOpen 5000 --chainGapExtend 1 --trimFraction 0.02 --logLevel DEBUG  >> /data/Arabidopsis/subtree1/anchorwave_chain.paf

#m=============================================================================================================================================================================================================================1

#genome_a   Ler.chr.all.v2.0.fasta   genome_b   An-1.chr.all.v2.0.fasta
anchorwave gff2seq -i /data/Arabidopsis/Ler.protein-coding.genes.v2.5.2019-10-09.gff3   -r /data/Arabidopsis/Ler.chr.all.v2.0.fasta  -o Ler_cds.fa
minimap2  -x splice -t 35 -k 12 -a -p 0.4 -N 20 /data/Arabidopsis/An-1.chr.all.v2.0.fasta  Ler_cds.fa > An-1_Ler_cds.sam
minimap2  -x splice -t 35 -k 12 -a -p 0.4 -N 20 /data/Arabidopsis/Ler.chr.all.v2.0.fasta  Ler_cds.fa > Ler_Ler_ref.sam
/usr/bin/time anchorwave genoAli -i /data/Arabidopsis/Ler.protein-coding.genes.v2.5.2019-10-09.gff3 -as Ler_cds.fa   -r /data/Arabidopsis/Ler.chr.all.v2.0.fasta -a An-1_Ler_cds.sam  -ar Ler_Ler_ref.sam -s /data/Arabidopsis/An-1.chr.all.v2.0.fasta -n /data/Arabidopsis/Ler_An-1.anchors -o /data/Arabidopsis/Ler_An-1.anchorwave.maf -f /data/Arabidopsis/Ler_An-1.anchorwave.f.maf -t 5 -IV  > Ler_An-1.avx2.log 2>&1
maf-convert sam /data/Arabidopsis/Ler_An-1.anchorwave.maf  > /data/Arabidopsis/Ler_An-1.2.anchorwave.sam
samtools view -H --reference /data/Arabidopsis/Ler.chr.all.v2.0.fasta  /data/Arabidopsis/Ler_An-1.2.anchorwave.sam > /data/Arabidopsis/Ler_An-1.3.anchorwave.sam
grep -v "^@" /data/Arabidopsis/Ler_An-1.2.anchorwave.sam >> /data/Arabidopsis/Ler_An-1.3.anchorwave.sam
k8 /data/paftools.js   sam2paf /data/Arabidopsis/Ler_An-1.3.anchorwave.sam > /data/Arabidopsis/Ler_An-1.tmp.anchorwave.paf
python /data/replace_ref_que.py -i /data/Arabidopsis/Ler_An-1.tmp.anchorwave.paf -o /data/Arabidopsis/Ler_An-1.anchorwave.paf -r Ler.protein -q An-1.protein

#m=============================================================================================================================================================================================================================1

#genome_a   Cvi.chr.all.v2.0.fasta   genome_b   An-1.chr.all.v2.0.fasta
anchorwave gff2seq -i /data/Arabidopsis/Cvi.protein-coding.genes.v2.5.2019-10-09.gff3   -r /data/Arabidopsis/Cvi.chr.all.v2.0.fasta  -o Cvi_cds.fa
minimap2  -x splice -t 35 -k 12 -a -p 0.4 -N 20 /data/Arabidopsis/An-1.chr.all.v2.0.fasta  Cvi_cds.fa > An-1_Cvi_cds.sam
minimap2  -x splice -t 35 -k 12 -a -p 0.4 -N 20 /data/Arabidopsis/Cvi.chr.all.v2.0.fasta  Cvi_cds.fa > Cvi_Cvi_ref.sam
/usr/bin/time anchorwave genoAli -i /data/Arabidopsis/Cvi.protein-coding.genes.v2.5.2019-10-09.gff3 -as Cvi_cds.fa   -r /data/Arabidopsis/Cvi.chr.all.v2.0.fasta -a An-1_Cvi_cds.sam  -ar Cvi_Cvi_ref.sam -s /data/Arabidopsis/An-1.chr.all.v2.0.fasta -n /data/Arabidopsis/Cvi_An-1.anchors -o /data/Arabidopsis/Cvi_An-1.anchorwave.maf -f /data/Arabidopsis/Cvi_An-1.anchorwave.f.maf -t 5 -IV  > Cvi_An-1.avx2.log 2>&1
maf-convert sam /data/Arabidopsis/Cvi_An-1.anchorwave.maf  > /data/Arabidopsis/Cvi_An-1.2.anchorwave.sam
samtools view -H --reference /data/Arabidopsis/Cvi.chr.all.v2.0.fasta  /data/Arabidopsis/Cvi_An-1.2.anchorwave.sam > /data/Arabidopsis/Cvi_An-1.3.anchorwave.sam
grep -v "^@" /data/Arabidopsis/Cvi_An-1.2.anchorwave.sam >> /data/Arabidopsis/Cvi_An-1.3.anchorwave.sam
k8 /data/paftools.js   sam2paf /data/Arabidopsis/Cvi_An-1.3.anchorwave.sam > /data/Arabidopsis/Cvi_An-1.tmp.anchorwave.paf
python /data/replace_ref_que.py -i /data/Arabidopsis/Cvi_An-1.tmp.anchorwave.paf -o /data/Arabidopsis/Cvi_An-1.anchorwave.paf -r Cvi.protein -q An-1.protein
cat /data/Arabidopsis/Cvi_An-1.anchorwave.paf     /data/Arabidopsis/Ler_An-1.anchorwave.paf > /data/Arabidopsis/cat_Ler_An-1_Cvi_An-1.2.anchorwave.paf
cat /data/Arabidopsis/cat_Ler_An-1_Cvi_An-1.2.anchorwave.paf > /data/Arabidopsis/Cvi_An-1.2.anchorwave_invert.paf
paf_invert -i /data/Arabidopsis/cat_Ler_An-1_Cvi_An-1.2.anchorwave.paf >> /data/Arabidopsis/Cvi_An-1.2.anchorwave_invert.paf
paf_chain -i /data/Arabidopsis/Cvi_An-1.2.anchorwave_invert.paf --maxGapLength 1000000 --chainGapOpen 5000 --chainGapExtend 1 --trimFraction 0.02 --logLevel DEBUG  >> /data/Arabidopsis/subtree1/anchorwave_chain.paf
paf_tile -i /data/Arabidopsis/subtree1/anchorwave_chain.paf --logLevel DEBUG  > /data/Arabidopsis/subtree1/anchorwave_tile.paf
grep -v tp:A:S /data/Arabidopsis/subtree1/anchorwave_tile.paf > /data/Arabidopsis/N2_primary.paf || true
grep  tp:A:S /data/Arabidopsis/subtree1/anchorwave_tile.paf> /data/Arabidopsis/N2_secondary.paf || true
cactus_consolidated --sequences " Arabidopsis_thaliana.protein /data/Arabidopsis/Arabidopsis_thaliana.protein.sanitized.fa An-1.protein /data/Arabidopsis/An-1.protein.sanitized.fa Ler.protein /data/Arabidopsis/Ler.protein.sanitized.fa Cvi.protein /data/Arabidopsis/Cvi.protein.sanitized.fa "  --speciesTree  " (Ler.protein:0.00493038,(Cvi.protein:0.0145906,(Arabidopsis_thaliana.protein:0.0117518,An-1.protein:0.00833626)N2:0.0070939)N1:0.00493038)N0;" --logLevel DEBUG --alignments  /data/Arabidopsis/N2_primary.paf   --params  /workflow/envs/cactus_progressive_config.xml --outputFile /data/Arabidopsis/N2.c2h  --outputHalFastaFile /data/Arabidopsis/N2.c2h.fa --outputReferenceFile /data/Arabidopsis/N2.ref --outgroupEvents Ler.protein Cvi.protein --referenceEvent N2 --t 36 --secondaryAlignments /data/Arabidopsis/N2_secondary.paf
mkdir /data/Arabidopsis/subtree2
#align Cvi.chr.all.v2.0.fasta	against	N2

#m=============================================================================================================================================================================================================================2

anchorwave gff2seq -i /data/Arabidopsis/Cvi.protein-coding.genes.v2.5.2019-10-09.gff3   -r /data/Arabidopsis/Cvi.chr.all.v2.0.fasta  -o Cvi_cds.fa
minimap2 -x splice -t 35 -k 12 -a -p 0.4 -N 20 /data/Arabidopsis/N2.ref Cvi_cds.fa > N2_Cvi_cds.sam
minimap2 -x splice -t 35 -k 12 -a -p 0.4 -N 20 /data/Arabidopsis/Cvi.chr.all.v2.0.fasta  Cvi_cds.fa > Cvi_Cvi_ref.sam
/usr/bin/time anchorwave proali -i /data/Arabidopsis/Cvi.protein-coding.genes.v2.5.2019-10-09.gff3 -as Cvi_cds.fa   -r /data/Arabidopsis/Cvi.chr.all.v2.0.fasta -a N2_Cvi_cds.sam  -ar Cvi_Cvi_ref.sam -s /data/Arabidopsis/N2.ref -n /data/Arabidopsis/Cvi_N2.anchors -o /data/Arabidopsis/Cvi_N2.anchorwave.maf -f /data/Arabidopsis/Cvi_N2.anchorwave.f.maf -Q 1 -R 1 -t 10   > Cvi_N2.log 2>&1
maf-convert sam /data/Arabidopsis/Cvi_N2.anchorwave.maf > /data/Arabidopsis/Cvi_N2.anchorwave.sam
samtools view -H --reference /data/Arabidopsis/Cvi.chr.all.v2.0.fasta  /data/Arabidopsis/Cvi_N2.anchorwave.sam > /data/Arabidopsis/Cvi_N2.2.anchorwave.sam
grep -v "^@" /data/Arabidopsis/Cvi_N2.anchorwave.sam >> /data/Arabidopsis/Cvi_N2.2.anchorwave.sam
k8 /data/paftools.js  sam2paf /data/Arabidopsis/Cvi_N2.2.anchorwave.sam > /data/Arabidopsis/Cvi_N2.tmp.anchorwave.paf
python /data/replace_ref_que.py -i /data/Arabidopsis/Cvi_N2.tmp.anchorwave.paf -o /data/Arabidopsis/Cvi_N2.anchorwave.paf -r Cvi.protein -q N2
cat /data/Arabidopsis/Cvi_N2.anchorwave.paf > /data/Arabidopsis/Cvi_N2.anchorwave_invert.paf
paf_invert -i /data/Arabidopsis/Cvi_N2.anchorwave.paf >> /data/Arabidopsis/Cvi_N2.anchorwave_invert.paf
paf_chain -i /data/Arabidopsis/Cvi_N2.anchorwave_invert.paf --maxGapLength 1000000 --chainGapOpen 5000 --chainGapExtend 1 --trimFraction 0.02 --logLevel DEBUG  >> /data/Arabidopsis/subtree2/anchorwave_chain.paf

#m=============================================================================================================================================================================================================================2

#genome_a   Ler.chr.all.v2.0.fasta   genome_b   Cvi.chr.all.v2.0.fasta
anchorwave gff2seq -i /data/Arabidopsis/Ler.protein-coding.genes.v2.5.2019-10-09.gff3   -r /data/Arabidopsis/Ler.chr.all.v2.0.fasta  -o Ler_cds.fa
minimap2  -x splice -t 35 -k 12 -a -p 0.4 -N 20 /data/Arabidopsis/Cvi.chr.all.v2.0.fasta  Ler_cds.fa > Cvi_Ler_cds.sam
minimap2  -x splice -t 35 -k 12 -a -p 0.4 -N 20 /data/Arabidopsis/Ler.chr.all.v2.0.fasta  Ler_cds.fa > Ler_Ler_ref.sam
/usr/bin/time anchorwave genoAli -i /data/Arabidopsis/Ler.protein-coding.genes.v2.5.2019-10-09.gff3 -as Ler_cds.fa   -r /data/Arabidopsis/Ler.chr.all.v2.0.fasta -a Cvi_Ler_cds.sam  -ar Ler_Ler_ref.sam -s /data/Arabidopsis/Cvi.chr.all.v2.0.fasta -n /data/Arabidopsis/Ler_Cvi.anchors -o /data/Arabidopsis/Ler_Cvi.anchorwave.maf -f /data/Arabidopsis/Ler_Cvi.anchorwave.f.maf -t 5 -IV  > Ler_Cvi.avx2.log 2>&1
maf-convert sam /data/Arabidopsis/Ler_Cvi.anchorwave.maf  > /data/Arabidopsis/Ler_Cvi.2.anchorwave.sam
samtools view -H --reference /data/Arabidopsis/Ler.chr.all.v2.0.fasta  /data/Arabidopsis/Ler_Cvi.2.anchorwave.sam > /data/Arabidopsis/Ler_Cvi.3.anchorwave.sam
grep -v "^@" /data/Arabidopsis/Ler_Cvi.2.anchorwave.sam >> /data/Arabidopsis/Ler_Cvi.3.anchorwave.sam
k8 /data/paftools.js   sam2paf /data/Arabidopsis/Ler_Cvi.3.anchorwave.sam > /data/Arabidopsis/Ler_Cvi.tmp.anchorwave.paf
python /data/replace_ref_que.py -i /data/Arabidopsis/Ler_Cvi.tmp.anchorwave.paf -o /data/Arabidopsis/Ler_Cvi.anchorwave.paf -r Ler.protein -q Cvi.protein
cat /data/Arabidopsis/Ler_Cvi.anchorwave.paf > /data/Arabidopsis/cat_Ler_Cvi_Ler_Cvi.2.anchorwave.paf
cat /data/Arabidopsis/cat_Ler_Cvi_Ler_Cvi.2.anchorwave.paf > /data/Arabidopsis/Ler_Cvi.2.anchorwave_invert.paf
paf_invert -i /data/Arabidopsis/cat_Ler_Cvi_Ler_Cvi.2.anchorwave.paf >> /data/Arabidopsis/Ler_Cvi.2.anchorwave_invert.paf
paf_chain -i /data/Arabidopsis/Ler_Cvi.2.anchorwave_invert.paf --maxGapLength 1000000 --chainGapOpen 5000 --chainGapExtend 1 --trimFraction 0.02 --logLevel DEBUG  >> /data/Arabidopsis/subtree2/anchorwave_chain.paf

#m=============================================================================================================================================================================================================================2

#genome_a   Ler.chr.all.v2.0.fasta   genome_b   N2
anchorwave gff2seq -i /data/Arabidopsis/Ler.protein-coding.genes.v2.5.2019-10-09.gff3   -r /data/Arabidopsis/Ler.chr.all.v2.0.fasta  -o Ler_cds.fa
minimap2   -x splice -t 35 -k 12 -a -p 0.4 -N 20 /data/Arabidopsis/N2.ref Ler_cds.fa > N2_Ler_cds.sam
minimap2   -x splice -t 35 -k 12 -a -p 0.4 -N 20 /data/Arabidopsis/Ler.chr.all.v2.0.fasta  Ler_cds.fa > Ler_Ler_ref.sam
/usr/bin/time anchorwave proali -i /data/Arabidopsis/Ler.protein-coding.genes.v2.5.2019-10-09.gff3 -as Ler_cds.fa   -r /data/Arabidopsis/Ler.chr.all.v2.0.fasta -a N2_Ler_cds.sam  -ar Ler_Ler_ref.sam -s /data/Arabidopsis/N2.ref -n /data/Arabidopsis/Ler_N2.anchors -o /data/Arabidopsis/Ler_N2.anchorwave.maf -f /data/Arabidopsis/Ler_N2.anchorwave.f.maf -Q 1 -R 1 -t 10   > Ler_N2.log 2>&1
maf-convert sam /data/Arabidopsis/Ler_N2.anchorwave.maf  > /data/Arabidopsis/Ler_N2.2.anchorwave.sam
samtools view -H --reference /data/Arabidopsis/Ler.chr.all.v2.0.fasta  /data/Arabidopsis/Ler_N2.2.anchorwave.sam > /data/Arabidopsis/Ler_N2.3.anchorwave.sam
grep -v "^@" /data/Arabidopsis/Ler_N2.2.anchorwave.sam >> /data/Arabidopsis/Ler_N2.3.anchorwave.sam
k8 /data/paftools.js   sam2paf /data/Arabidopsis/Ler_N2.3.anchorwave.sam > /data/Arabidopsis/Ler_N2.tmp.anchorwave.paf
python /data/replace_ref_que.py -i /data/Arabidopsis/Ler_N2.tmp.anchorwave.paf -o /data/Arabidopsis/Ler_N2.anchorwave.paf -r Ler.protein -q N2
cat /data/Arabidopsis/Ler_N2.anchorwave.paf > /data/Arabidopsis/cat_Ler_N2_Ler_N2.2.anchorwave.paf
cat /data/Arabidopsis/cat_Ler_N2_Ler_N2.2.anchorwave.paf > /data/Arabidopsis/Ler_N2.2.anchorwave_invert.paf
paf_invert -i /data/Arabidopsis/cat_Ler_N2_Ler_N2.2.anchorwave.paf >> /data/Arabidopsis/Ler_N2.2.anchorwave_invert.paf
paf_chain -i /data/Arabidopsis/Ler_N2.2.anchorwave_invert.paf --maxGapLength 1000000 --chainGapOpen 5000 --chainGapExtend 1 --trimFraction 0.02 --logLevel DEBUG  >> /data/Arabidopsis/subtree2/anchorwave_chain.paf
paf_tile -i /data/Arabidopsis/subtree2/anchorwave_chain.paf --logLevel DEBUG  > /data/Arabidopsis/subtree2/anchorwave_tile.paf
grep -v tp:A:S /data/Arabidopsis/subtree2/anchorwave_tile.paf > /data/Arabidopsis/N1_primary.paf || true
grep  tp:A:S /data/Arabidopsis/subtree2/anchorwave_tile.paf> /data/Arabidopsis/N1_secondary.paf || true
cactus_consolidated --sequences " Cvi.protein /data/Arabidopsis/Cvi.protein.sanitized.fa N2 /data/Arabidopsis/N2.ref Ler.protein /data/Arabidopsis/Ler.protein.sanitized.fa "  --speciesTree  " (Ler.protein:0.00493038,(Cvi.protein:0.0145906,N2:0.0070939)N1:0.00493038)N0;" --logLevel DEBUG --alignments  /data/Arabidopsis/N1_primary.paf   --params  /workflow/envs/cactus_progressive_config.xml --outputFile /data/Arabidopsis/N1.c2h  --outputHalFastaFile /data/Arabidopsis/N1.c2h.fa --outputReferenceFile /data/Arabidopsis/N1.ref --outgroupEvents Ler.protein --referenceEvent N1 --t 36 --secondaryAlignments /data/Arabidopsis/N1_secondary.paf
mkdir /data/Arabidopsis/subtree3
#align Ler.chr.all.v2.0.fasta	against	N1

#m=============================================================================================================================================================================================================================3

anchorwave gff2seq -i /data/Arabidopsis/Ler.protein-coding.genes.v2.5.2019-10-09.gff3   -r /data/Arabidopsis/Ler.chr.all.v2.0.fasta  -o Ler_cds.fa
minimap2 -x splice -t 35 -k 12 -a -p 0.4 -N 20 /data/Arabidopsis/N1.ref Ler_cds.fa > N1_Ler_cds.sam
minimap2 -x splice -t 35 -k 12 -a -p 0.4 -N 20 /data/Arabidopsis/Ler.chr.all.v2.0.fasta  Ler_cds.fa > Ler_Ler_ref.sam
/usr/bin/time anchorwave proali -i /data/Arabidopsis/Ler.protein-coding.genes.v2.5.2019-10-09.gff3 -as Ler_cds.fa   -r /data/Arabidopsis/Ler.chr.all.v2.0.fasta -a N1_Ler_cds.sam  -ar Ler_Ler_ref.sam -s /data/Arabidopsis/N1.ref -n /data/Arabidopsis/Ler_N1.anchors -o /data/Arabidopsis/Ler_N1.anchorwave.maf -f /data/Arabidopsis/Ler_N1.anchorwave.f.maf -Q 1 -R 1 -t 10   > Ler_N1.log 2>&1
maf-convert sam /data/Arabidopsis/Ler_N1.anchorwave.maf > /data/Arabidopsis/Ler_N1.anchorwave.sam
samtools view -H --reference /data/Arabidopsis/Ler.chr.all.v2.0.fasta  /data/Arabidopsis/Ler_N1.anchorwave.sam > /data/Arabidopsis/Ler_N1.2.anchorwave.sam
grep -v "^@" /data/Arabidopsis/Ler_N1.anchorwave.sam >> /data/Arabidopsis/Ler_N1.2.anchorwave.sam
k8 /data/paftools.js  sam2paf /data/Arabidopsis/Ler_N1.2.anchorwave.sam > /data/Arabidopsis/Ler_N1.tmp.anchorwave.paf
python /data/replace_ref_que.py -i /data/Arabidopsis/Ler_N1.tmp.anchorwave.paf -o /data/Arabidopsis/Ler_N1.anchorwave.paf -r Ler.protein -q N1
cat /data/Arabidopsis/Ler_N1.anchorwave.paf > /data/Arabidopsis/Ler_N1.anchorwave_invert.paf
paf_invert -i /data/Arabidopsis/Ler_N1.anchorwave.paf >> /data/Arabidopsis/Ler_N1.anchorwave_invert.paf
paf_chain -i /data/Arabidopsis/Ler_N1.anchorwave_invert.paf --maxGapLength 1000000 --chainGapOpen 5000 --chainGapExtend 1 --trimFraction 0.02 --logLevel DEBUG  >> /data/Arabidopsis/subtree3/anchorwave_chain.paf
paf_tile -i /data/Arabidopsis/subtree3/anchorwave_chain.paf --logLevel DEBUG  > /data/Arabidopsis/subtree3/anchorwave_tile.paf
grep -v tp:A:S /data/Arabidopsis/subtree3/anchorwave_tile.paf > /data/Arabidopsis/N0_primary.paf || true
grep  tp:A:S /data/Arabidopsis/subtree3/anchorwave_tile.paf> /data/Arabidopsis/N0_secondary.paf || true
cactus_consolidated --sequences " Ler.protein /data/Arabidopsis/Ler.protein.sanitized.fa N1 /data/Arabidopsis/N1.ref "  --speciesTree  " (Ler.protein:0.00493038,N1:0.00493038)N0;" --logLevel DEBUG --alignments  /data/Arabidopsis/N0_primary.paf   --params  /workflow/envs/cactus_progressive_config.xml --outputFile /data/Arabidopsis/N0.c2h  --outputHalFastaFile /data/Arabidopsis/N0.c2h.fa --outputReferenceFile /data/Arabidopsis/N0.ref --referenceEvent N0 --t 36 --secondaryAlignments /data/Arabidopsis/N0_secondary.paf
halAppendCactusSubtree /data/Arabidopsis/N0.c2h  /data/Arabidopsis/N0.c2h.fa  "(Ler.protein:0.00493038,N1:0.00493038)N0;"    result/evolverPlants.hal
halAppendCactusSubtree /data/Arabidopsis/N1.c2h  /data/Arabidopsis/N1.c2h.fa  "(Cvi.protein:0.0145906,N2:0.0070939)N1;"    result/evolverPlants.hal --outgroups Ler.protein
halAppendCactusSubtree /data/Arabidopsis/N2.c2h  /data/Arabidopsis/N2.c2h.fa  "(Arabidopsis_thaliana.protein:0.0117518,An-1.protein:0.00833626)N2;"    result/evolverPlants.hal --outgroups Ler.protein,Cvi.protein
