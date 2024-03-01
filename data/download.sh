mkdir Arabidopsis
cd Arabidopsis

wget -c https://1001genomes.org/data/MPIPZ/MPIPZJiao2020/releases/current/strains/An-1/An-1.chr.all.v2.0.fasta.gz --no-check-certificate
wget -c https://1001genomes.org/data/MPIPZ/MPIPZJiao2020/releases/current/strains/An-1/An-1.protein-coding.genes.v2.5.2019-10-09.CDS.fasta.gz --no-check-certificate
wget -c https://1001genomes.org/data/MPIPZ/MPIPZJiao2020/releases/current/strains/An-1/An-1.protein-coding.genes.v2.5.2019-10-09.gff3.gz --no-check-certificate


wget -c https://1001genomes.org/data/MPIPZ/MPIPZJiao2020/releases/current/strains/Cvi/Cvi.chr.all.v2.0.fasta.gz --no-check-certificate
wget -c https://1001genomes.org/data/MPIPZ/MPIPZJiao2020/releases/current/strains/Cvi/Cvi.protein-coding.genes.v2.5.2019-10-09.CDS.fasta.gz --no-check-certificate--no-check-certificate
wget -c https://1001genomes.org/data/MPIPZ/MPIPZJiao2020/releases/current/strains/Cvi/Cvi.protein-coding.genes.v2.5.2019-10-09.gff3.gz --no-check-certificate


wget -c https://1001genomes.org/data/MPIPZ/MPIPZJiao2020/releases/current/strains/Ler/Ler.chr.all.v2.0.fasta.gz --no-check-certificate           
wget -c https://1001genomes.org/data/MPIPZ/MPIPZJiao2020/releases/current/strains/Ler/Ler.protein-coding.genes.v2.5.2019-10-09.CDS.fasta.gz --no-check-certificate
wget -c https://1001genomes.org/data/MPIPZ/MPIPZJiao2020/releases/current/strains/Ler/Ler.protein-coding.genes.v2.5.2019-10-09.gff3.gz --no-check-certificate

wget -c https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz --no-check-certificate
wget -c https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/arabidopsis_thaliana/cds/Arabidopsis_thaliana.TAIR10.cds.all.fa.gz	--no-check-certificate
wget -c https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-56/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.56.gff3.gz --no-check-certificate

find -maxdepth 1 -name "*.gz"|xargs -i gunzip {}

mkdir CDS

mv Arabidopsis_thaliana.TAIR10.cds.all.fa Arabidopsis_thaliana.TAIR10.CDS.fasta
sed -i 's/^>\([0-9]*\).*/>chr\1/g' Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
sed -i 's/^\([0-9]*\)/chr\1/g' Arabidopsis_thaliana.TAIR10.56.gff3
cp *.CDS.fasta  CDS
cp ../../workflow/scripts/CombineCDS.py CDS
python CDS/CombineCDS.py > ../non_duplicate_CDS.fa

awk -i inplace '/^>chr[1]/{p=1; print; next} p && /^>/{exit} p' An-1.chr.all.v2.0.fasta
awk -i inplace '/^>chr[1]/{p=1; print; next} p && /^>/{exit} p' Cvi.chr.all.v2.0.fasta
awk -i inplace '/^>chr[1]/{p=1; print; next} p && /^>/{exit} p' Ler.chr.all.v2.0.fasta
awk -i inplace '/^>chr[1]/{p=1; print; next} p && /^>/{exit} p' Arabidopsis_thaliana.TAIR10.dna.toplevel.fa



sed -i -E  '/^#|^chr[1]/!d' An-1.protein-coding.genes.v2.5.2019-10-09.gff3
sed -i -E  '/^#|^chr[1]/!d' Cvi.protein-coding.genes.v2.5.2019-10-09.gff3
sed -i -E  '/^#|^chr[1]/!d' Ler.protein-coding.genes.v2.5.2019-10-09.gff3
sed -i -E  '/^#|^chr[1]/!d' Arabidopsis_thaliana.TAIR10.56.gff3