cd ..
cd xiphophorus_maculatus/
ls
mv xiphophorus_maculatus_vep_79_Xipmac4.4.2.tar.gz ..
cd ..
ls
tar xvzf xiphophorus_maculatus_vep_79_Xipmac4.4.2.tar.gz 
ls
rm*.gz
rm *.gz
ls
cd /applications/local/ensembl-tools-release-79/scripts/variant_effect_predictor/
variant_effect_predictor.pl -i example.vcf --cache
ls
variant_effect_predictor.pl -i example_GRCh37.vcf --cache
cd
variant_effect_predictor.pl -i /applications/local/ensembl-tools-release-79/scripts/variant_effect_predictor/example_GRCh37.vcf --cache
ls
cat Lo
cat Log.*
ls
rm variant_effect_output.txt*
bcftools 
bgzip 
tabix 
locate varutils.pl
su brewmaster
which picard
vcftools --version
vcftools
aptitude search picard
su brewmaster
htslib
cd Desktop/freebayes-0.9.20/
ls
make
cd ..
git clone --recursive git://github.com/ekg/freebayes.git
cd freebayes
git checkout b04023680443f48fc1eceee6fb5e58377ccc6be9
ls
make
ls
cd vcflib/
ls
ls bin/
make
cd ..
ls
less README.md 
su brewmaster
make
ls bin
cd vcflib/
ls
less README.md 
su brewmaster
variant_effect_predictor.pl -i example.vcf --cache
freebayes 
df -h
ls
cd Desktop/
ls
tar xvzf RNA-seq.tgz 
tar xvzf ChIP-seq_0914.tar.gz 
aptitude search fastx
fastqc
cutadapt 
cutadapt --version
aptitude search cutadapt
aptitude search sdfsdf
aptitude search cut
aptitude search cutadapt
which cutadapt 
ll /usr/local/bin/cutadapt
apt-file search cutadapt
cutadapt --version
trim_galore 
trim_galore --version
python -v
python --version
bwa --versipon
bwa --version
bwa
bamtools --version
aptitude search bamtools
which bamtools
ll /usr/local/bin/bamtools
ll /usr/local/bin/bamtools-2.3.0 
apt-file search bamtools
samtools --version
samtools tview 
bwa mem 
igv
gatk
GATK 
cd Desktop/GenomeAnalysisTK-3.3-0/
java -jar GenomeAnalysisTK.jar 
GATK
su brewmaster
ipython notebook
fastqc 
which fastqc
ls /applications/local/FastQC/
cat /applications/local/FastQC/fastqc
head /applications/local/FastQC/fastqc
su brewmaster
pymol 
rasmol
locate modeller
mod9.14 
modeller 
aptitude search pymol
su brewmaster
pymol
su brewmaster
cd Desktop/vcftools_0.1.12b/
make
ls
ls bin/
ls perl/
java -v
java -version
which jav
which java
bwa
samtools
picard 
locate picard
locate picard.jar
locate picard-tools
ls /applications/local/picard-tools/
vim .bashrc
echo $PICARD
java -jar $PICARD
GATK 
GATK -h
GATK -h | less
vim .bashrc
$GATK -h
java -jar $GATK -h
java -jar $GATK
java -jar $PICARD
java -jar $GATK
echo $GATK
ls /applications/local/GATK/
java -jar $GATK -h
Rscript 
aptitude search java
locate java7
/usr/share/doc/oracle-java7-installer
ls /usr/share/doc/oracle-java7-installer/
java --version
java -version
aptitude search java
aptitude search java | grep &i
aptitude search java | grep $i
aptitude search java | grep ^i
java -version
ls
ls tools/
ls tools/picard/
ls cambridge_mda14/
ls cambridge_mda14/calling/
ll
ll -rt
ls -lrt
rm -r cambridge_mda14
rm -r tools/
ls
su brewmaster
mv Course_Materials/GATK_var .
mv Course_Materials/GATK_basic .
su brewmaster
ls
cd Desktop/
git pull https://github.com/bioinformatics-core-shared-training/basicr.git
git clone https://github.com/bioinformatics-core-shared-training/basicr.git
cd ..
ls
ls R_course/
rm -r R_course/
mv Desktop/basicr/ R_course
ls /mounts/Open_Share/
su brewmaster
R
cd Course_Materials/
ls
rm -r ./*
ls
git clone https://github.com/pycam/python-intro.git
ls
sudo mv python-intro/* .
mv python-intro/* .
ls
rm -r python-intro/
y
ls
cat TOC.md 
less TOC.md 
ls
cd /data
cd databases/
ls -l
cd databases/
ls -l
ls
ls -l
rm *~
ls -l
man rsync
ssh -X dpjudge@hancocks
su brewmaster 
cat > foo.fasta
blastp 
blastp -help
blastp -db UniRef90 foo.fasta 
blastp -db UniRef90 -in foo.fasta 
blastp -db UniRef90 -query foo.fasta 
rm foo.fasta 
ls
ls Course_Materials/
cd
ls
cat Log.progress.out 
ls -l
cat Log.
cat Log.out 
ls 
rm Log.*
ls Local_Disk/
ls -ld Local_Disk/
ls -l Local_Disk/
ls -lR Local_Disk/
rm -r Local_Disk/
mount -l
mount -l | grep Local
ls -l Local_Disk/
mount -l | grep Local
cd Local_Disk/
top
ls -l Local_Disk/
su brewmaster
du -sh /data/databases/UniRef/
ls /databases/
ls /databases/BLASTDB/
ls -l /databases/BLASTDB/
ls -l /databases/
cd /databases/
ln -s /data/databases/UniRef/ .
su brewmaster 
su brewmaster 
macs2
macs2 --version
ls /applications/local/picard-tools/
ls $PICARD
samtools --version
bedtools --version
which bedtools
aptitude search tophat
aptitude search tophat2
aptitude search bowtie2
which bedtools
ll /applications/local/bin/bedtools
which ngsplotdb.py 
ll /applications/local/bin/ngsplotdb.py
which cufflinks
aptitude search cufflinks
apt-show-versions cufflinks
bedtools --version
echo $PICARD_HOME
ls $PICARD_HOME
bedtools --version
tophat2 --version
tophat --version
which tophat
bowtie2 --version
bowtie2
which bedGraphToBigWig 
bedGraphToBigWig --version
bedGraphToBigWig
ll /applications/local/bin/bedGraphToBigWig
ls /applications/local/UCSC_Software/
genomeCoverageBed
which genomeCoverageBed
which genomeCoverageBed | ll
which genomeCoverageBed
ll `which genomeCoverageBed`
which htseq-count 
htseq-count --version
htseq-count
cd Desktop/
mkdir ngsplor
mv ngsplor/ ngsplot
cd ngsplot/
git clone https://github.com/shenlab-sinai/ngsplot.git
vim ~/.bashrc
ls
cd ngsplot/
ls
head README
ls bin
cd
vim ~/.bashrc
ngsplotdb.py 
echo $NGSPLOT
ls $NGSPLOT
su brewmaster
cd Desktop/phantompeakqualtools/
ls
Rscript run_spp.R 
vim run_spp.R 
fastqc 
aptitude search fastx
apt-show-versions fastx-toolkit
bowtie2 --version
which bowtie2
aptitude search bowtie
apt-show-versions bowtie2
ls /applications/local/
cd ..
which tophat2
tophat2 --version
su brewmaster
vi ~/Desktop/ChIP-seq/commands.sh
source ~/Desktop/ChIP-seq/commands.sh
ls -l
ls -l bowtie_index/
mv Mus_musculus.GRCm38.dna.chromosome.1.fa bowtie_index/
vi ~/Desktop/ChIP-seq/commands.sh
ls -l bowtie_index/
head bowtie_index/Mus_musculus.GRCm38.dna.chromosome.1.fa 
vi ~/Desktop/ChIP-seq/commands.sh
rm -R bowtie_index/
source ~/Desktop/ChIP-seq/commands.sh
ls -l bo
ls -l 
head mm10.fa 
vi ~/Desktop/ChIP-seq/commands.sh
rm -R bowtie_index/
cd ..
rm -R bowtie_index/
source ~/Desktop/ChIP-seq/commands.sh
head mm10.fa 
head Mus_musculus.GRCm38.dna.chromosome.1.fa 
tail mm10.fa 
tail Mus_musculus.GRCm38.dna.chromosome.1.fa 
vi ~/Desktop/ChIP-seq/commands.sh
cd ..
ls -l
mv Oct4_ESC_rep1_chr1.fastq Oct4.fastq 
vi ~/Desktop/ChIP-seq/commands.sh
source ~/Desktop/ChIP-seq/commands.sh
ls -l
ls -l bowtie_index/
rm bowtie_index/Mus_musculus.GRCm38.dna.chromosome.1.fa
igv.sh &
igv
vi ~/Desktop/ChIP-seq/commands.sh
source ~/Desktop/ChIP-seq/commands.sh
ls -l
vi ~/Desktop/ChIP-seq/commands.sh
source ~/Desktop/ChIP-seq/commands.sh
vi ~/Desktop/ChIP-seq/commands.sh
source ~/Desktop/ChIP-seq/commands.sh
ls -l
head Oct4_peaks.narrowPeak 
vi ~/Desktop/ChIP-seq/commands.sh
source ~/Desktop/ChIP-seq/commands.sh
ls -l
vi ~/Desktop/ChIP-seq/commands.sh
source ~/Desktop/ChIP-seq/commands.sh
rm H3K4me3.bw 
source ~/Desktop/ChIP-seq/commands.sh
ls -l
cp -R ~/Desktop/Course_Materials/PeakAnalyzer_1.4/ .
cd PeakAnalyzer_1.4/
java -jar PeakAnalyzer.jar &
cd ../
less -S Oct4_peaks.narrowPeak 
awk -OFS=FS=$'\t' 'print $1,$2,$3,$4,$5' Oct4_peaks.narrowPeak | head
awk -OFS=FS=$'\t' '{print $1,$2,$3,$4,$5}' Oct4_peaks.narrowPeak | head
awk --OFS=FS=$'\t' '{print $1,$2,$3,$4,$5}' Oct4_peaks.narrowPeak | head
awk --OFS=FS='\t' '{print $1,$2,$3,$4,$5}' Oct4_peaks.narrowPeak | head
awk -OFS=FS='\t' '{print $1,$2,$3,$4,$5}' Oct4_peaks.narrowPeak | head
awk 'BEGIN{FS=OFS="\t"};{print $1,$2,$3,$4,$5}' Oct4_peaks.narrowPeak | head
awk 'BEGIN{FS=OFS="\t"};{print $1,$2,$3,$4,$5}' Oct4_peaks.narrowPeak > Oct4_peaks.bed
cd PeakAnalyzer_1.4/
java -jar PeakAnalyzer.jar &
awk 'BEGIN{FS=OFS="\t"};{print $1,$2,$3}' Oct4_peaks.narrowPeak > Oct4_peaks.bed
awk 'BEGIN{FS=OFS="\t"};{print $1,$2,$3}' ../Oct4_peaks.narrowPeak > ../Oct4_peaks.bed
cd ../#]
cd ..
macs2 --help
macs2 callpeak --help
ls -l
bedGraphToBigWig Oct4_treat_pileup.bdg bowtie_index/mouse.mm10.genome > Oct4_treat_pileup.bw
less -S Oct4_treat_pileup.bdg
mv Oct4_treat_pileup.bdg Oct4_treat_pileup.bedGraph
bedGraphToBigWig Oct4_treat_pileup.bedGraph bowtie_index/mouse.mm10.genome > Oct4_treat_pileup.bw
bedGraphToBigWig Oct4_treat_pileup.bdg bowtie_index/mouse.mm10.genome > Oct4_treat_pileup.bw
bedGraphToBigWig Oct4_treat_pileup.sorted.bedGraph bowtie_index/mouse.mm10.genome  Oct4_treat_pileup.bw
ls -l
tail Oct4_treat_pileup.bedGraph 
head Oct4_treat_pileup.bedGraph 
rm Oct4_peaks.*
rm Oct4_control_lambda.bdg 
rm Oct4_treat_pileup.*
ls -l
rm Oct4_summits.bed 
rm Oct4_model.r 
ls -l
vi ~/Desktop/ChIP-seq/commands.sh
source ~/Desktop/ChIP-seq/commands.sh
ls -l
head Oct4_summits.bed
history | grep "sort"
sort -k5 -nr Oct4_summits.bed > Oct4_summits.sorted.bed 
head Oct4_summits.sorted.bed
awk 'BEGIN{FS="\t"}; NR < 301 { print $1, $2-31, $3+30 }' Oct4_summits.sorted.bed > Oct4_top300_summits.bed
head Oct4_top300_summits.bed
awk 'BEGIN{FS=OFS="\t"}; NR < 301 { print $1, $2-30, $3+29 }' Oct4_summits.sorted.bed > Oct4_top300_summits.bed
head Oct4_top300_summits.bed
i=10
while [[ $i -lt 100 ]]; do echo "Looking at ${i}% of the reads"; macs2 randsample -t Oct4.bam -p $i -o Oct4.perc${i}.bed --outdir; macs2_downsample -f BAM; macs2 randsample -t gfp.bam -p $i -o gfp.perc${i}.bed --outdir; macs2_downsample -f BAM; macs2 callpeak -t macs2_downsample/Oct4.perc${i}.bed -c; macs2_downsample/gfp.perc${i}.bed --gsize 138000000 --format BED; --name macs2_downsample/Oct4.perc${i} --pvalue 1e-3; ((i = i + 10))
while [[ $i -lt 100 ]]; do echo "Looking at ${i}% of the reads"; macs2 randsample -t Oct4.bam -p $i -o Oct4.perc${i}.bed --outdir macs2_downsample -f BAM; macs2 randsample -t Input.bam -p $i -o Input.perc${i}.bed --outdir macs2_downsample -f BAM; macs2 callpeak -t macs2_downsample/Oct4.perc${i}.bed -c macs2_downsample/Input.perc${i}.bed --gsize 138000000 --format BED --name macs2_downsample/Oct4.perc${i} --pvalue 1e-3 --call-summits; ((i = i + 10)); done
bedtools getfasta -fi bowtie_index/mm10_full.fa -bed Oct4_top300_summits.bed -fo Oct4_top300_summits.fa
samtools faidx bowtie_index/mm10_full.fa
bedtools getfasta -fi bowtie_index/mm10_full.fa -bed Oct4_top300_summits.bed -fo Oct4_top300_summits.fa
ls -l
ls -l bowtie_index/
samtools faidx bowtie_index/mm10_full.fa
ls -l bowtie_index/
bedtools getfasta -fi bowtie_index/mm10_full.fa -bed Oct4_top300_summits.bed -fo Oct4_top300_summits.fa
ls -l
head Oct4_top300_summits.fa 
ls -l
rm H3K4me3.bam H3K4me3.bedgraph H3K4me3.fastq H3K4me3.sam H3K4me3.sorted.bam* Input.b* Input.s* Oct4.b* Oct4.s* Oct4_* 
ls -l
rm -R macs2_downsample/
ls -l
ls -l bowtie_index/
rm bowtie_index/mm10*.ebwt
ls -l bowtie_index/
ls -l
ls -l PeakAnalyzer_1.4/
rm PeakAnalyzer_1.4/Oct4_peaks.bed 
ls -l PeakAnalyzer_1.4/
rm commands.sh 
ls -l
history | grep "sort"
history | grep "getfasta"
fetchChromSizes --help
fetchChromSizes
macs2 --HELP
macs2 callpeak --help
igv &
cd ../PeakAnalyzer_1.4/
java -jar PeakAnalyzer.jar 
macs2 callpeak --help
bedtools 
cd ..
ls -l
ls -l macs2_downsample/
su brewmaster
ls
fastx_trimmer -h
fastq_quality_trimmer -Q 33 -t 20 -l 50 -i bad_example.fastq -o bad_example_quality_trimmed.fastq
ls
fastq &
fastqc &
ls
rm bad_example_*
ls
cd ../ChIP-seq/
ls
ls bowtie_index/
wget http://www.ebi.ac.uk/~mxenoph/chip-seq-June15/buecker_2014/H3K27ac_ESC_rep1.sorted.bam
ls
ngs.plot.r
mkdir Histone_modifications
ls
mv H3K27ac_ESC_rep1.sorted.bam Histone_modifications/
cd Histone_modifications/
ls
ngs.plot.r -G mm10 -O H3K27ac.tss -FL 150 -R tss -C H3K27ac_ESC_rep1.sorted.bam -T H3K27ac
wget http://www.ebi.ac.uk/~mxenoph/chip-seq-June15/buecker_2014/ChIP_mE14_H3K4me3.sorted.bam
wget http://www.ebi.ac.uk/~mxenoph/chip-seq-June15/buecker_2014/Input_ESC_rep2.sorted.bam
ngs.plot.r -G mm10 -O H3K27ac.tss -FL 150 -R tss -C H3K27ac_ESC_rep1.sorted.bam:Input_ESC_rep2.sorted.bam -T H3K27ac:Control
wget http://www.ebi.ac.uk/~mxenoph/chip-seq-June15/buecker_2014/ChIP_mE14_H3K4me3.sorted.bam
ls
rm ChIP_mE14_H3K4me3.sorted.bam1
rm ChIP_mE14_H3K4me3.sorted.bam.1 
ls
mv H3K27ac_ESC_rep1.sorted.bam H3K27ac.bam
mv ChIP_mE14_H3K4me3.sorted.bam H3K4me3.bam
ls
wget http://www.ebi.ac.uk/~mxenoph/chip-seq-June15/buecker_2014/Input_ESC_rep1.sorted.bam
ls
mv Input_ESC_rep1.sorted.bam Input.bam
ls
rm Input_ESC_rep2.sorted.bam*
ls
ls H3K27ac_ESC_rep1.sorted.bam.*
rm H3K27ac_ESC_rep1.sorted.bam.*
ls
wget http://www.ebi.ac.uk/~mxenoph/chip-seq-June15/buecker_2014/H3K4me1_ESC_rep1.sorted.bam
ls
mv H3K4me1_ESC_rep1.sorted.bam H3K4me1.bam
ls
ngs.plot.r -G mm10 -O H3K27ac.tss -FL 150 -R tss -C H3K27ac.bam:Input.bam -T H3K27ac:Control
ngs.plot.r -G mm10 -O H3K4me3.tss -FL 150 -R tss -C H3K4me3.bam -T H3K4me3
ngs.plot.r -G mm10 -O H3K4me1.H3K4me3 -FL 150 -R tss -C multhist.txt -T H3K4me1.H3K4me3
ls
wget http://www.ebi.ac.uk/~mxenoph/chip-seq-June15/buecker_2014/Oct4_ESC_rep1.bam
rm Oct4_ESC_rep1.bam 
ls
ngs.plot.r -G mm10 -O Oct4.tss -FL 150 -R tss -C Oct4.bam -T Oct4
wget http://www.ebi.ac.uk/~mxenoph/chip-seq-June15/buecker_2014/Oct4_ESC_rep1.bam
mv Oct4_ESC_rep1.bam Oct4.bam
ngs.plot.r -G mm10 -O Oct4.tss -FL 150 -R tss -C Oct4.bam -T Oct4
ngs.plot.r -G mm10 -O Oct4.tss -FL 150 -R tss -C Oct4.bam:Input.bam -T Oct4:Input
ls
ls *.bam*
ls *.bam?
ls *bai
rm *bai
ls *cnt
rm *cnt
ls
rm *zip
rm *pdf
ls
rm multhist.txt 
ls
pwd
ls ..
ls
wget http://www.ebi.ac.uk/~mxenoph/chip-seq-June15/buecker_2014/ChIP_mE14_WCE.fastq
ls
rm ChIP_mE14_WCE.fastq 
exit
echo $NGSPLOT
ls /applications/local/ngsplot/ngsplot
ls /applications/local/ngsplot/ngsplot/bin/
ls /applications/local/ngsplot/
ls /applications/local/ngsplot/ngsplot/database/
ngsplotdb.py -y list
ls /applications/local/ngsplot/ngsplot/database/mm10/
ls /applications/local/ngsplot
ls /applications/local/ngsplot/ngsplot
ngsplotdb.py install ~/Desktop/ngsplotdb_mm10_75_3.00.tar.gz 
ngsplotdb.py -y list
wget http://www.ebi.ac.uk/~mxenoph/chip-seq-June15/QC/
ls ~/Desktop/
ls
cd ~/Desktop/
mkdir QC
ls
cd QC/
ls
wget http://www.ebi.ac.uk/~mxenoph/chip-seq-June15/QC/good_example.fastq
wget http://www.ebi.ac.uk/~mxenoph/chip-seq-June15/QC/bad_example.fastq
wget http://www.ebi.ac.uk/~mxenoph/chip-seq-June15/QC/adapter_contamination.fastq
ls
fastqc
fastx_trimmer -f 1 -l 80 bad_example.fastq -o bad_example_length_trimmed.fastq 
ls
fastx_trimmer -f 1 -l 80 bad_example.fastq -o bad_example_trimmed01.fastq 
fastx_trimmer -f 1 -l 80 -i bad_example.fastq -o bad_example_trimmed01.fastq 
ls
ngsplotdb.py -y list
ls
exit
vim .bashrc
ls
rm -r GATK_*
ls
su brewmaster
ls $IDR
echo $IDR
ls $IDR_HOME
ls $PHANTOMPEAK_HOME
ls $PEAKANALYZER_HOME
ls $PEAKANALYZER
java -jar $PEAKANALYZER
su brewmaster
ls
cd D
cd Desktop/
git clone https://github.com/cambiotraining/r-intro.git
cd r-intro/
ls
git fetch origin
git checkout gh-pages
ls
sudo vim index.html 
vim index.html 
git push origin gh-pages
git commit -a -m "added materials links"
man git commit
git commit -a -m "added materials links" --author="Paul Judge <pj237@cam.ac.uk>"
vim index.html 
su brewmaster
su brewmaster 
start_art 
which start_art 
cd Desktop/
ll
cat foo~ 
ls -l foo~ 
rm foo~ 
ls -l
rm report\(1\).csv~ 
ls -l
ls -la
rm .Rhistory 
ls -la
cd
ls -lart
Art
art
act
su brewmaster 
la
cd Desktop/
mkdir datSORT
mv Cambridge_2012_data/ datSORT/
mv workshop.tar.gz datSORT/
cd datSORT/
ls -l
mv Cambridge_2012_data/ 2012
mkdir 2013
mkdir 2014
mv workshop.tar.gz 2014/
cd 2014/
ls
tar xvzf workshop.tar.gz 
cd ..
ls *
rm 2014/workshop.tar.gz 
ls
ls *
ls ..
mv 2014/* 2013
ls *
diff -r 2012/Module_1 2013/Module_1
diff -r 2012/Module_2 2013/Module_2
diff -r 2012/Module_3 2013/Module_3
ls -l *
mv 2014/Course_Data/* 2014
rmdir 2014/Course_Data/
ls *
diff -r 2012/Module_2 2013/Module_2
diff -r 2012/Module_2 2014/Module_2
diff -r 2013/Module_3 2014/Module_3
diff -r 2012/Module_3 2014/Module_3
diff -r 2013/Module_3 2014/Module_3
diff -r 2012/Module_1 2014/Module_1
diff -r 2013/Module_1 2014/Module_1
diff -r 2012 2014
ls -l 
ls 2014
mv 2014 ../Course_Data
cd /mounts/hancocks/Administration/COURSE_RECORDS/POSTGRADUATE/
ls
mkdir -p Theoretical/Figure_Design/2015/05.27/Desktop
mv /mounts/Open_Share/Course_Survey.URL Theoretical/Figure_Design/2015/05.27/Desktop
su brewmaster
ipython notebook
which pip
ll /usr/bin/pip
pip --version
nano
ls
rm -r R_course/
ls
su brewmaster
cd Desktop/
ls
unzip alignment.zip 
mkdir alignment
mv *.gz alignment
mkdir prioritization
cd prioritization/
unzip ../prioritization .
unzip ../prioritization
ls
cd ..
mkdir annotation
cd annotation/
unzip ../annotation.zip 
cd ..
mkdir rna_seq
cd rna_seq/
unzip ../rna_seq
cd ..
mkdir visualization
cd visualization/
unzip ../visualization.zip 
cd ..
ls
cd Course_
cd Course_Materials
cd calling/
ls
man unzip
mkdir gatk java-6 picard example1 genome mutect example3
unzip gatk.zip -d gatk
unzip java-6.zip -d java-6
ls
unzip genome.zip -d genome
ls
rm -r genome
ls
unzip picard.zip picard
unzip picard.zip -d picard
ls
unzip example1.zip -d example1
cd example1/
ls
rm 000-dna_chr21_100_hq_pe.bam 
unzip ../example1.zip 
ls
ll
rm 000-dna_chr21_100_hq_pe.bam 
ls
cd ..
ls
ls gatk
ls java-6
ls picard
ls
ls mutect
unzip mutect.zip -d mutect
ls
unzip example3.zip -d example3
samtools --version
bwa --version
bwa
tophat
tophat --version
tophat2 --version
bowtie--version
bowtie2 --version
STAR 
STAR --version
cd
cd Desktop/
mkdir STAR
cd STAR/
git clone https://github.com/alexdobin/STAR.git
ls
cd STAR/
ls
cat RE
cat README 
ls
cd source/
ls
make STAR
cd ..
ls bin/
ls bin/Linux_x86_64
ll bin/Linux_x86_64
STARlo
which STAR
aptitude search star
aptitude search star | grep bio
aptitude search star | grep align
aptitude search STAR | grep align
ls
which STAR
su brewmaster
STAR --version
STARlong 
ls /applications/local/
ls /applications/local/IDR/
fastqc 
passwd 
ifconfig
ll
pwd
#ll Course_Materials/annotation/
htop
samtools
java -versin
java -version
which java
ll /usr/bin/
ll /etc/alternatives/
lsb_release -a
w
df
df -h
ll
ls
ll bin/
ll
ls
cd igv/
ls
cd ..
ls
bwa 
bwa mem
bowtie
bowtie2
bowtie2-align-
bowtie2-align
bowtie2-build
bowtie2-
STAR
samtools
bowtie2 -v
bowtie2 -version
bowtie2 --version
igv
cd Course_Materials/annotation/
ll
cd cellbase/
ll
ln -s cellbasemaster/bin/cellbase.sh cellbase.sh
ll
ll cellbase
ll cellbase.sh 
ll examples/
./cellbase.sh variant-annotation -h
cd cellbasemaster/
cd bin/
./cellbase.sh variant-annotation -h
cd ..
echo $PATH
#export PATH=$PATH:
pwd
export PATH=$PATH:/home/participant/Course_Materials/annotation/cellbase/cellbasemaster/bin/
cellbase.sh variant-annotation -h
cd ..
ll
cd annovar/
ll
cd ..
#ll visualization/exam
igv
hpg-aligner 
hpg-aligner -v
hpg-aligner -version
hpg-aligner --version
hpg-aligner -h
tophat2
tophat2 --version
tree 
htop
locate *bam
cp /usr/local/lib/R/site-library/Rsamtools/extdata/tiny.bam /tmp/
cd /tmp/
ll
samtools sort 
samtools sort tiny.bam -o tiny_sorted.bam
ll
samtools sort 
samtools sort -o tiny_sorted.bam tiny.bam 
samtools view 
samtools view tiny.bam 
samtools view 
samtools view tiny.bam -o pepe
ll
less pepe 
samtools index 
samtools index tiny.bam tiny.bam.bai
samtools index tiny.bam > tiny.bam.bai
ll
llll
ls
pwd
fastqc
fastqc -v
fastqc --verison
which fastqc
fastqc -v
cutadapt
cutadapt --version
fastqc -v
ll
ls 
cd Desktop/
ls
cd Course_Materials
ls
ls quality_control/
ls
cd ..
cd ~
mkdir david_tests
cd david_tests/
ll
git clone https://github.com/ngs-course 
git https://github.com/ngs-course/ngs-course.github.io
git clone https://github.com/ngs-course/ngs-course.github.io
ll
mkdir -p /tmp/cambridge_mda/
cp -r Course_Materials/annotation /tmp/cambridge_mda/
cd /tmp/cambridge_mda/annotation/annovar/
perl convert2annovar.pl -format vcf4 example/example1.vcf > example/example1.annovar
#perl convert2annovar.pl -format vcf4 example/example1.vcf > example/example1.annovar
cat example/example1.annovar
perl convert2annovar.pl -format vcf4 -allsample -withfreq example/example1.vcf > example/example1.annovar
cat example/example1.annovar
cat example/example1.vcf
java -version
which java
ll /usr/bin/java
perl annotate_variation.pl -buildver hg19 -downdb snp135 humandb/
perl annotate_variation.pl -buildver hg19 -downdb 1000g2012apr humandb/
perl annotate_variation.pl -buildver hg19 -downdb refGene humandb/
perl annotate_variation.pl -buildver hg19 -downdb phastConsElements46way humandb/
perl annotate_variation.pl -buildver hg19 -downdb genomicSuperDups humandb/
perl annotate_variation.pl -buildver hg19 -downdb ljb2_all humandb/
perl annotate_variation.pl -buildver hg19 -downdb genomicSuperDups humandb/
perl annotate_variation.pl -buildver hg19 -downdb ljb26_all humandb/
perl annotate_variation.pl -buildver hg19 -downdb esp6500si_all humandb/
perl annotate_variation.pl -buildver hg19 -downdb esp6500siv2_all humandb/
mkdir results
perl annotate_variation.pl --geneanno example/example1.annovar humandb/ -build hg19 --outfile results/0-geneanno
cat results/0-geneanno.*
perl annotate_variation.pl -regionanno -dbtype cytoBand example/example1.annovar humandb/ -build hg19 --outfile results/1-regionanno
perl annotate_variation.pl -filter -dbtype snp142 example/example1.annovar humandb/ -build hg19 --outfile results/2-filter
perl annotate_variation.pl -buildver hg19 -downdb cytoBand humandb/
cd /tmp/cambridge_mda/
cd annovar
cd annotation/annovar/
perl annotate_variation.pl -buildver hg19 -downdb cytoBand humandb/
perl annotate_variation.pl -buildver hg19 -downdb snp142 humandb/
perl annotate_variation.pl -filter -dbtype 1000g2012apr_all -maf 0.01 example/example1.annovar humandb/ -build hg19 --outfile results/2-filter
#perl annotate_variation.pl -filter -dbtype 1000g2012apr_all -maf 0.01 example/example1.annovar humandb/ -build hg19 --outfile results/2-filter
cat results/2-filter.*
cat results/1-regionanno.hg19_cytoBand 
cd
cd Course_Materials/annotation/
ll
cd cellbase/
ll
rm cellbase.sh 
ln -s cellbasedevelop/bin/cellbase.sh cellbase.sh
ll
./cellbase.sh -h
cd
#cp -r /home/participant/Desktop/Open_Share/annotation /home/participant/cambridge_mda/
cd /tmp/cambridge_mda/annotation/cellbase/
./cellbase.sh -h
ll
cd
cd Course_Materials/annotation/cellbase/
ll
rm cellbase.sh 
cd /tmp/cambridge_mda/annotation/cellbase/
ll
rm cellbase.sh 
export PATH=$PATH:/tmp/cambridge_mda/annotation/cellbase/cellbasedevelop/bin
cellbase.sh -h
cellbase.sh variant-annotation -h
mkdir cellbase_results    
cellbase.sh variant-annotation -i /home/participant/cambridge_mda/annotation/cellbase/examples/CEU.exon.2010_03.genotypes.vcf -o /home/participant/cambridge_mda/annotation/cellbase_results/CEU.exon.2010_03.annotated.vep -s hsapiens -u wwwdev.ebi.ac.uk
#cellbase.sh variant-annotation -i /tmp/cambridge_mda/annotation/cellbase/examples/CEU.exon.2010_03.genotypes.vcf -o /tmp/cambridge_mda/annotation/cellbase_results/CEU.exon.2010_03.annotated.vep -s hsapiens -u wwwdev.ebi.ac.uk
#ll cellbase_results/
cellbase.sh variant-annotation -i /tmp/cambridge_mda/annotation/cellbase/examples/CEU.exon.2010_03.genotypes.vcf -o /tmp/cambridge_mda/annotation/cellbase_results/CEU.exon.2010_03.annotated.vep -s hsapiens -u wwwdev.ebi.ac.uk
#cellbase.sh variant-annotation -i /tmp/cambridge_mda/annotation/cellbase/examples/CEU.exon.2010_03.genotypes.vcf -o /tmp/cambridge_mda/annotation/cellbase_results/CEU.exon.2010_03.annotated.vep -s hsapiens -u wwwdev.ebi.ac.uk
mkdir results
cellbase.sh variant-annotation -i /tmp/cambridge_mda/annotation/cellbase/examples/CEU.exon.2010_03.genotypes.vcf -o /tmp/cambridge_mda/annotation/cellbase/results/CEU.exon.2010_03.annotated.vep -s hsapiens -u wwwdev.ebi.ac.uk
sfasd
cd 
cd Desktop/Course_Materials/visualization/
ll
cd ..
ll
cd /home/participant/Desktop/Course_Materials/visualization
cd example_0
ll
samtools index igv1.bam
ll
rm *.bai
cd /home/participant/Desktop/Course_Materials/visualization/example_1
samtools index NA12878_child.bam
samtools index NA12891_dad.bam
samtools index NA12892_mom.bam
ll
rm *.bai
cd /home/participant/Desktop/Course_Materials/visualization/example_2
samtools index heart.bodyMap.bam
samtools index brain.bodyMap.bam
ll
rm *.bai
top
rm -rf /tmp/cambridge_mda
ll
ls
mkdir nacho
wget https://www.dropbox.com/sh/4qkqch7gyt888h7/AAD6uW2mF1yRN-C4j8COVYpGa/alignment/dna_chr21_100_hq_read1.fastq.gz?dl=0
wget https://www.dropbox.com/sh/4qkqch7gyt888h7/AACDss3p99EZDU-vi08qVK1na/alignment/dna_chr21_100_hq_read2.fastq.gz?dl=0
ll
cd nacho/
ll
cd ..
ls -ltr
mv dna_chr21_100_hq_read* nacho/
cd nacho/
ll
mv dna_chr21_100_hq_read1.fastq.gz\?dl\=0 dna_chr21_100_hq_read1.fastq.gz
mv dna_chr21_100_hq_read2.fastq.gz\?dl\=0 dna_chr21_100_hq_read2.fastq.gz
ll
gunzip dna_chr21_100_hq_read1.fastq.gz 
gunzip dna_chr21_100_hq_read2.fastq.gz 
ll
bwa 
bwa mem
bwa mem -h
bwa mem
wget https://www.dropbox.com/sh/4qkqch7gyt888h7/AAAZwGqS-KBt4RAN9l-r9UaAa/f000_chr21_ref_genome_sequence.fa?dl=0
ll
zless f000_chr21_ref_genome_sequence.fa\?dl\=0 
mv f000_chr21_ref_genome_sequence.fa\?dl\=0 chr21_ref_genome_sequence.fa
ll
bwa index 
bwa index chr21_ref_genome_sequence.fa 
ll
rm chr21_ref_genome_sequence.fa.*
ll
mkdir index
mv chr21_ref_genome_sequence.fa index/
bwa index index/chr21_ref_genome_sequence.fa 
cat /mounts/Course_Share/calling/java-6/LICENSE
cd
cat Desktop/Course_Share/calling/java-6/LICENSE 
cat Desktop/Course_Share/calling/java-6/README 
cat Desktop/Course_Share/calling/java-6/Welcome.html 
cat Desktop/Course_Share/calling/java-6/plugin/desktop/sun_java.desktop 
ls -l Desktop/Course_Share/calling/java-6/plugin/desktop/sun_java.desktop 
ls -l Desktop/Course_Share/calling/java-6/plugin/desktop/
ls -l Desktop/Course_Share/calling/java-6/plugin/
diff -r /mounts/Course_Share/ /home/participant/Desktop/Course_Materials
diff -r /mounts/Course_Share/ /home/participant/Desktop/Course_Materials | less
ls -ld /mounts/Course_Share/calling/java-6/man/man1/
git clone https://github.com/pycam/python-intro.git
mv python-intro/ Course_Materials/
cd Course_Materials/
ls
mv python-intro/* .
rm -r python-intro/
y
l
su brewmaster
python
pwd
git clone https://github.com/bioinformatics-core-shared-training/cruk-bioinf-sschool.git
ls cruk-bioinf-sschool/
ls
cd cruk-bioinf-sschool/
ls
ls -l
participant@exmoor:~/cruk-bioinf-sschool$ ls -l
wget https://www.dropbox.com/s/b8gix98mzlzdrqq/SRR576933.fastq.gz -P Day1
fastqc Day1/SRR576933.fastq.gz 
R
cd ..
rm -r cruk-bioinf-sschool/
rm -rf cruk-bioinf-sschool/
ls
git clone https://github.com/bioinformatics-core-shared-training/cruk-bioinf-sschool.git
ls
rm -rf cruk-bioinf-sschool/
git clone https://github.com/bioinformatics-core-shared-training/cruk-bioinf-sschool.git
cd cruk-bioinf-sschool/
cat softwareInstall.sh 
ls
ls -l
chmod 755 getData.sh 
./getData.sh 
samtools view -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/other_exome_alignments/HG00096/exome_alignment/HG00096.mapped.illumina.mosaik.GBR.exome.20111114.bam 22 | samtools view -bS - > Day2/HG00096.chr22.bam
samtools index Day2/HG00096.chr22.bam
./getData.sh 
pwd
ls
ls -l ref_data/
pwd
cd..
cd ..
rm -rf cruk-bioinf-sschool/
git clone https://github.com/bioinformatics-core-shared-training/cruk-bioinf-sschool.git
cd cruk-bioinf-sschool/
chmod 755 getData.sh 
./getData.sh 
cd ..
rm -rf cruk-bioinf-sschool/
git clone https://github.com/bioinformatics-core-shared-training/cruk-bioinf-sschool.git
chmod 755 getData.sh 
cd cruk-bioinf-sschool/
chmod 755 getData.sh 
./getData.sh 
ls -l 
ls -l Software/
cd ..
rm -f cruk-bioinf-sschool/
rm -rf cruk-bioinf-sschool/
git clone https://github.com/bioinformatics-core-shared-training/cruk-bioinf-sschool.git
cd cruk-bioinf-sschool/
chmod 755 getData.sh 
./getData.sh 
R -f installBiocPkgs.R 
R
ls
cd ref_data/
mkdir bwa
ln -s chr22.fa bwa/
ls -l bwa/
bwa index
pwd
`pwd
'pwd
'pwd'
ln -s 'pwd'/chr22.fa bwa/
rm bwa/chr22.fa 
ln -s 'pwd'/chr22.fa bwa/
ls -l bwa/
'pwd'
rm bwa/chr22.fa 
ln -s $(pwd)/chr22.fa bwa/
ls -l bwa/
bwa index bwa/chr22.fa 
ln -s $(pwd)/chr22.fa bwa/
mkdir bowtie
ln -s $(pwd)/chr22.fa bwa/
mkdir bowtie
ln -s $(pwd)/chr22.fa bowtie/
../Software/bowtie2-2.2.5/bowtie2-build chr22.fa chr22
cd ..
pwd
rm -rf cruk-bioinf-sschool/
git clone https://github.com/bioinformatics-core-shared-training/cruk-bioinf-sschool.git
rm -rf cruk-bioinf-sschool/
git clone https://github.com/bioinformatics-core-shared-training/cruk-bioinf-sschool.git
cd cruk-bioinf-sschool/
chmod 755 getData.sh 
./getData.sh 
pwd
ls -l Software/
cd ..
rm -rf cruk-bioinf-sschool/
git clone https://github.com/bioinformatics-core-shared-training/cruk-bioinf-sschool.git
cd cruk-bioinf-sschool/
chmod 755 getData.sh 
./getData.sh 
tophat2 
cd ref_data/
../Software/bowtie2-2.2.5/bowtie2-build bowtie/chr22.fa bowtoe/chr22
../Software/bowtie2-2.2.5/bowtie2-build bowtie/chr22.fa bowtie/chr22
cd ..
rm -rf cruk-bioinf-sschool/
git clone https://github.com/bioinformatics-core-shared-training/cruk-bioinf-sschool.git
cd cruk-bioinf-sschool/
ls
ls -l
chmod 755 getData.sh 
./getData.sh 
cd ..
rm -rf cruk-bioinf-sschool/
git clone https://github.com/bioinformatics-core-shared-training/cruk-bioinf-sschool.git
cd cruk-bioinf-sschool/
chmod 755 getData.sh 
./getData.sh 
cd Day1/
fastqc data/sample.fastq 
fastq_quality_filter -v -q 20 -p 75 -i data/sample.fastq -o data/sample.filtered.fastq
cd data/
ls
ls -l
fastq_quality_filter -v -q 20 -p 75 -i sample.fastq -o sample.filtered.fastq 
fastq_quality_filter 
fastx_trimmer -v -f 7 -l 36 -i sample.fastq 
fastx_trimmer -v -f 7 -l 36 -i sample.fastq -o sample.trimmted.fastq
fastq_quality_filter -h
fastq_quality_filter -v -q 20 -i sample.fastq 
fastq_quality_filter -v -q 20 -p 75 -i sample.fastq -o sample.filtered.fastq 
fastq_quality_filter -v -q 20 -p 75 -i sample.fastq -o sample_filtered.fastq
fastq_quality_filter -v -Q 33 -q 20 -p 75 -i sample.fastq -o sample_filtered.fastq
fastq_quality_filter -v -Q 64 -q 20 -p 75 -i sample.fastq -o sample_filtered.fastq
cd ..
ls
bwa aln test.reads_1.fq ../ref_data/chr22.fa 
bwa aln
bwa aln ../ref_data/chr22.fa test.reads_1.fq 
bwa aln ../ref_data/bwa/chr22.fa test.reads_1.fq 
bwa aln ../ref_data/bwa/chr22.fa test.reads_1.fq > 1.sai
bwa aln ../ref_data/bwa/chr22.fa test.reads_1.fq > test.reads_1.sai
bwa aln ../ref_data/bwa/chr22.fa test.reads_2.fq > test.reads_2.sai
bwa sampe ../ref_data/bwa/chr22.fa test.reads_1.sai test.reads_2.sai test.reads_1.fq test.reads_2.fq > test.sam
head test.sam 
samtools view -bS test.sam > test.bam
ls -lrt
cd ..
rm -rf cruk-bioinf-sschool/
git clone https://github.com/bioinformatics-core-shared-training/cruk-bioinf-sschool.git
cd cruk-bioinf-sschool/
chmod 755 getData.sh 
./getData.sh 
wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.5/bowtie2-2.2.5-linux-x86_64.zip
cd ..
rm -rf cruk-bioinf-sschool/
git clone https://github.com/bioinformatics-core-shared-training/cruk-bioinf-sschool.git
cd cruk-bioinf-sschool/
chmod 755 getData.sh 
./getData.sh 
cd Software/
wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.5/bowtie2-2.2.5-linux-x86_64.zip
cd ..
./getData.sh 
cd ..
rm -rf cruk-bioinf-sschool/
git clone https://github.com/bioinformatics-core-shared-training/cruk-bioinf-sschool.git
cd Software/
cd cruk-bioinf-sschool/
chmod 755 getData.sh 
./getData.sh 
ls
ls -l Day2/
ls -l
cd /home/participant/Desktop/DATA_FOR_DAY2/bam/
ls
samtools view -h 16N_aligned.bam chr22 | samtools view -bS - > 16N_chr22.bam
ls -lrt
samtools view 16N_chr22.bam | head
bamToFastq
bamToFastq -i 16N_chr22.bam -fq 16N_reads.fq
ls -lrt
rm 16N_chr22.bam
head 16N_reads.fq 
cp 16N_reads.fq /home/participant/cruk-bioinf-sschool/Day2/
cd /home/participant/cruk-bioinf-sschool/Day2/
ls
bowtie2
bowtie2 ../ref_data/bowtie/chr22 -U 16N_reads.fq 
bowtie2 -x ../ref_data/bowtie/chr22 -U 16N_reads.fq 
la
ls
docker
bwa
samtools
fast
fastqc 
igv
R
which fastqc
fastqc --version
cd Course_Materials/
git clone https://github.com/bioinformatics-core-shared-training/cruk-bioinf-sschool.git
ls
cp -r cruk-bioinf-sschool/* .
ls
rmdir cruk-bioinf-sschool/
rm -r cruk-bioinf-sschool/
ls
chmod +x getData.sh 
./getData.sh 
vim getData.sh 
su brewmaster
which bismark
ll /applications/local/bin/bismark
ls /applications/local/Bismark/
bismark --version
su brewmaster
bismark --version
samtools --version
bowtie2 --version
trim_galore --version
seqmonk --version
seqmonk 
ls /applications/local/SeqMonk/SeqMonk/Resources/monk.svg
fastqc
passwd
ls
ls -l
ls -l cruk-bioinf-sschool/
ls -l cruk-bioinf-sschool/Day2/
ls -l cruk-bioinf-sschool/Day2/bam/
ls -l cruk-bioinf-sschool/ref_data/
cutadapt 
ls -l
ls Course_Materials/
ls -l Course_Materials/
ls -l Course_Materials/Day1/
ls -l Course_Materials/Day1/alignment-demo/
ls -l */*
ls -l Course_Materials/*/*
ls -l Course_Materials/Day3/
ls -l Course_Materials/Day2
ls Course_Materials/ref_data/
R
ls -l Course_Materials/Day1
cd Course_Materials/Day1/alignment-demo/
ls
ls -l
cd ..
cd qa/
ls
ls -l
cd ..
ls
ls -l
cd ..
ls 
ls ref_data/
bwa
samtools
cd ..
cd Course_Materials/Day1/
ls
cd ..
ls Day2/bam/
cd Day2/bam/
ls
pwd
macs2
cd ..
cd ../Day4/
ls
cd data_for_practical/
ls
macs2 callpeak -t TF_1.bam -c Input.bam
ls -lrt
cd ..
cd Day1/
ls
cd qa/
ls
fastqc sample.fastq 
fastx_trimmer -v -Q 64 -q 20 -p 75 -i sample.fastq -o sample_filtered.fastq
fastq_quality_filter -v -Q 64 -q 20 -p 75 -i sample.fastq -o sample_filtered.fastq
cutadapt -m 20 -e 0.1 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACACA sample2.fastq -o sample2--cutadapt.fastq
cd ..
sl
ls
cd ..
ls 
ls -l
cd Day2/
LS
ls
sudo apt-get install sra-toolkit
fastq-dump
cd ../ref_data/
ls
mkdir whole_genome
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
ls
mv human_g1k_v37.fasta.gz whole_genome/
ls
cd whole_genome/
bwa index human_g1k_v37.fasta.gz 
ls
ls -l
cd ..
cd Day1/alignment-demo/
ls
fastq-dump -A SRR1186252
fastq-dump
fastq-dump SRR1186252
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX483%2FSRX483591/SRR1186252/SRR1186252.sra
ls
fastq-dump SRR1186252.sra 
ls -l
ls -lh SRR1186252.
ls -lh SRR1186252.fastq 
ls -l /home/participant/Course_Materials/Day1/alignment-demo/
bwa aln SRR1186252.fastq ../../ref_data/whole_genome/human_g1k_v37.fasta.gz
bwa aln ../../ref_data/whole_genome/human_g1k_v37.fasta.gz SRR1186252.fastq 
bwa aln ../../ref_data/whole_genome/human_g1k_v37.fasta.gz SRR1186252.fastq > SRR1186252.sai
ls
rm SRR1186252.sai
cd ..
ls -lrt
cd ..
ls -lrt
exit
ls
cd Course_Materials/
ls
cd Day2/bam/
ls
exit
ls
ls -lrt
cd Course_Materials/
ls
ls -lrt
ls -lrt ref_data/
ls -lrt ref_data/whole_genome/
ls
cd cruk-bioinf-sschool/
ls
vi  courseInstall.sh 
ls -lh Day2/
which tophat
which bowtie2
which bowtie
ls b
ls ba
lscd D
cd Day2/
ls -lnh
ls -lh bam/
ls
echo $PATH
which BOWTIE_INDEX
whereisBOWTIE_INDEX
whereis BOWTIE_INDEX
whereis BOWTIE2_INDEX
ls
cd cruk-bioinf-sschool/
cd Day2/
ls
fastqc
fastqc 16N_reads.fq 
ls
pwd
wc -l 16N_reads.fq 
less 16N_reads.fq 
ls lh
ls -lh
echo $PATH
ls /home/participant/bin/
ls /home/participant/bin/ngs_align 
ls /applications/local/binpref/
ls
ls -lh
cd ..
ls
ls ref_data/
ls ref_data/bowtie/
bowtie2 -h
bowtie2-build
cd Day2
bowtie2-build -help
bowtie2-build
ls $HOME
ls -lh "$HOME/cruk-bioinf-sschool/ref_data/bowtie/"
ls
ls ../ref_data/
tophat --help
tophat2 --help
tophat --help
tophat2 --help
ls
cp ./ref_data/chr22.fa ../ref_data/bowtie/
cp ../ref_data/chr22.fa ../ref_data/bowtie/
export BOWTIE2_INDEXES="$HOME/cruk-bioinf-sschool/ref_data/bowtie/"
echo $BOWTIE2_INDEXES
mkdir TopHat2_Alignment
tophat2 --output-dir Tophat2_Alignment chr22 16N_reads.fq
ls
du -h Tophat2_Alignment/
du -rh Tophat2_Alignment/
ls -lh Tophat2_Alignment/
pwd
R
ls
cd ..
ls
ls Day2/bam/
samtools flagstat Tophat2_Alignment/accepted_hits.bam
cd Day2
samtools flagstat Tophat2_Alignment/accepted_hits.bam
samtools flagstat bam/16N.bam
ls bam/16N_aligned.bamsamtools flagstat bam/16N_aligned.bam
samtools flagstat bam/16N_aligned.bam
ls TopHat2_Alignment/
samtools flagstat Tophat2_Alignment/unmapped.bam 
samtools flagstat bam/16N_aligned.bam
samtools flagstat Tophat2_Alignment/16N.bam
samtools flagstat Tophat2_Alignment/unmapped.bam 
ls Tophat2_Alignment/
samtools index Tophat2_Alignment/accepted_hits.bam
samtools flagstat Tophat2_Alignment/accepted_hits.bam
samtools flagstat Tophat2_Alignment/unmapped.bam
cufflinks -h
cufflinks bam/16N_aligned.bam
cufflinks Tophat2_Alignment/accepted_hits.bam
ls
less transcripts.gtf 
less genes.fpkm_tracking 
less isoforms.fpkm_tracking 
ls
R
cd cruk-bioinf-sschool/
ls
fastqc 16N_reads.fq
ls ../ref_data/bowtie
cp ../ref_data/chr22.fa ../ref_data/bowtie/
bowtie2-build #Command to build index
export BOWTIE2_INDEXES="$HOME/cruk-bioinf-sschool/ref_data/bowtie/"
echo $BOWTIE2_INDEXES
mkdir Tophat2_Alignment
tophat2 --output-dir Tophat2_Alignment chr22 16N_reads.fq
samtools index Tophat2_Alignment/accepted_hits.bam
samtools flagstat Tophat2_Alignment/accepted_hits.bam
samtools flagstat Tophat2_Alignment/unmapped.bam
ls bam/
samtools flagstat bam/16N_aligned.bam
samtools view -H bam/16N_aligned.bam
cufflinks Tophat2_Alignment/accepted_hits.bam > Tophat2_Alignment/ #This will take ~7-8 minutes
cufflinks Tophat2_Alignment/accepted_hits.bam > Tophat2_Alignment/cufflinksResult #This will take ~7-8 minutes
ls Tophat2_Alignment/
ls Tophat2_Alignment/cufflinksResult 
ls -lh  Tophat2_Alignment/
ls
ls -lh
library(Rsubread)
filesToCount <- dir("bam", pattern=".bam$", full.names=T)
R
ls
ls -lrt Course_Materials/
ls -lrt Course_Materials/ref_data/
ls -lrt
ls cruk-bioinf-sschool/
ls -lrt cruk-bioinf-sschool/
ls -lrt cruk-bioinf-sschool/Day2/
ls -lrt cruk-bioinf-sschool/
ls -lrt cruk-bioinf-sschool/Day1/
ls -lrt cruk-bioinf-sschool/Day1/alignment-demo/
ls -lrt Course_Materials/
ls -l Course_Materials/
ls -lrt cruk-bioinf-sschool/Day2/
ls -l Course_Materials/
ls -l
ls -l Course_Materials/
ls -l Course_Materials/ref_data/
ls -l Course_Materials/Day1/
cd Course_Materials/Day1/alignment-demo/
ls
rm test.reads_*.fq
ls
cd ..
ls
ls -l
ls -lrth
mv chr22.* /mounts/Open_Share/
cd Course_Materials/
ls
mkdir git
cd git/
cit clone https://github.com/bioinformatics-core-shared-training/cruk-bioinf-sschool.git
git clone https://github.com/bioinformatics-core-shared-training/cruk-bioinf-sschool.git
diff . ..
dif cruk-bioinf-sschool/ ..
diff cruk-bioinf-sschool/ ..
diff cruk-bioinf-sschool/ .. | grep -v <|>
diff cruk-bioinf-sschool/ .. | grep -v <,>
diff cruk-bioinf-sschool/ .. | grep -v <
diff cruk-bioinf-sschool/ .. | grep -v \<
diff cruk-bioinf-sschool/ .. | grep -v \<|\>
diff cruk-bioinf-sschool/ .. | grep -v \<||\>
diff cruk-bioinf-sschool/ .. | grep -v \<|grep -v \>
diff cruk-bioinf-sschool/ .. | grep -v \<|grep -v \> | less
ls
cp -r cruk-bioinf-sschool/* ..
ls ..
ls ../Day1/
diff cruk-bioinf-sschool/ .. | grep -v \<|grep -v \> | less
diff cruk-bioinf-sschool/ ..| less
cd ..
rm -r git/
ls
su brewmaster
ls
ls GenomeAnnotation-Prac/
passwd -d
passwd -d participant
su brewmaster
ls
ls -lrt
mv GenomeAnnotation-Prac/ /mounts/Open_Share/
mv cruk-bioinf-sschool/ /mounts/Open_Share/
uname -r
uname
uname -a
ls
ls -l
pwd
cd ~/.
cd /
ls
cd /home
ls
cd participant/
cd ..
clear
motd
ps auxf
clear
tmux
vim
ls
cd participant/
ls
cd bin/
ls
cd ..
l
cd R/
ls
cd x86_64-pc-linux-gnu-library/
ls
cd 3.2/
ls
mkdir ~/db
ls > 
ls > ~/db/Rpackages
ls
cd /
cd usr/
cd lib
ls
cd R/
ls
cd modules/
ls
cd ..
cd lib
LS
ls
cd ..
cd lib
cd ..
ls
cd lib
ls
cd ..
cd library/
ls
cd ..
grep DESeq
grep -R DEseq .
grep -R 'DEseq' .
cd site-library/
ls
ls > ~/db/Rsystempackages
cd ~/db
ls
mv Rpackages Rusrpackages
cat Rusrpackages 
cat Rsystempackages 
ls
cat 150727.Rhistory 
ls
ls -latr
cat ~/.ssh/
ls
ls ~/.ssh/
ls
whcih fastqc
ls
cd ..
ls
grep fastqc
grep -R 'fastqc' .
clear
fastqc
ls
cd participant/
ls
cd Course_Materials/Day1/qa/
ls
fastqc sample.fastq
ls
ls -latr
firefox
firefox sample_fastqc.html 
fastq_quality_filter -v -Q 64 -q 20 -p 75 -i sample.fastq -o sample_filtered.fastq
ls -latrt
fastx_trimmer -v -f 7 -l 36 -i sample_filtered.fastq -o sample_filtered_and_trimmed.fastq
fastqc sample2.fastqc
ls
fastqc sample2.fastq
ls
sample2_fastqc.html
grep 'TruSeq Adapter, Index 5'
grep --hrlp
grep --help
grep 'TruSeq Adapter, Index 5' sample2_fastqc.html 
ls
grep -R 'TruSeq' sample2_fastqc.html 
ls
cutadapt -m 20 -e 0.1 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACACA sample2.fastq -o sample2--cutadapt.fastq
fastqc sample2--cutadapt.fastq 
firefox sample2--cutadapt_fastqc.html 
cat ~/.bash_history 
ls
cd ..
ls
cd Course_Materials/
ls
cd Day1
ls
cd alignment-demo/
ls
ls -latr
ls
bwa --help
bwa
bwa index
bwa index -p hg19chr6bwaidx -a bwtsw hg19chr6.fa
bwa aln -t hg19chr6bwaidx SRR1186252_trimmed.fq.chr6.fq > SRR1186252_trimmed.fq.chr6.fq.bwa
ls -latr
bwa aln -t hg19chr6bwaidx.sa SRR1186252_trimmed.fq.chr6.fq > SRR1186252_trimmed.fq.chr6.fq.bwa
bwa aln -t 4 hg19chr6bwaidx.sa SRR1186252_trimmed.fq.chr6.fq > SRR1186252_trimmed.fq.chr6.fq.bwa
bwa aln -t 4 hg19chr6bwaidx SRR1186252_trimmed.fq.chr6.fq > SRR1186252_trimmed.fq.chr6.fq.bwa
ls -latr
bwa samse hg19chr6bwaidx SRR1186252_trimmed.fq.chr6.fq.bwa SRR1186252_trimmed.fq.chr6.fq > SRR1186252_trimmed.fq.chr6.fq.sam
ls -l
head SRR1186252_trimmed.fq.chr6.fq.sam
samtools view -bS SRR1186252_trimmed.fq.chr6.fq.sam > SRR1186252_trimmed.fq.chr6.fq.chr6.fq.bam
ls -latr
samtools sort -0 bam -o SRR1186252_trimmed.fq.chr6.fq.sorted.bam -T temp SRR1186252_trimmed.fq.chr6.fq.bam
samtools sort -O bam -o SRR1186252_trimmed.fq.chr6.fq.sorted.bam -T temp SRR1186252_trimmed.fq.chr6.fq.bam
samtools sort -O bam -o SRR1186252_trimmed.fq.chr6.fq.sorted.bam -T temp SRR1186252_trimmed.fq.chr6.fq.chr6.fq.bam 
samtools index SRR1186252_trimmed.fq.chr6.fq.sorted.bam 
samstat SRR1186252_trimmed.fq.chr6.fq.sorted.bam
ssh -P 58581 duncan@enttzar.co.uk
ssh duncan@enttzar.co.uk -P 58581
clear
export
print $TERM
echo $TERM
wget http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.12.tar.bz2/download
ls
ls -latr
ls
cd Downloads/
ls
cd ..
ls -latr
wget http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.12.tar.bz2/download .
wget http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.12.tar.bz2/download
ls
rm download.*
ls -latr
cd Downloads/
echo $PATH
wget http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.12.tar.bz2/download ./Download/bwa-0.7.12.tar.bz2
wget http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.12.tar.bz2/download ./Downloads/bwa-0.7.12.tar.bz2
wget http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.12.tar.bz2/download -o /Downloads/bwa-0.7.12.tar.bz2
wget http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.12.tar.bz2/download -o /Download/bwa-0.7.12.tar.bz2
wget http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.12.tar.bz2/download -o ./Download/bwa-0.7.12.tar.bz2
wget http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.12.tar.bz2/download -o Download/bwa-0.7.12.tar.bz2
wget --?
wget --help

wget http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.12.tar.bz2/download -O bwa-0.7.12.tar.bz2
ls
rm download
rm download.1 
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
ls
cp .bash_history ~/db/
cd db
ls
cp .bash_history ~/db/bash_history
ls
ssh-keygen -t rsa -b 4096 -C "e@e.com"
ssh-add ~/.ssh/id_rsa
cat ~/.ssh/id_rsa.pub 
cd ..
ls
cd Course_Materials/Day1/alignment-demo/
ls
samstat SRR1186252_trimmed.fq.chr6.fq.sorted.bam
cp .bash_history ~/db/bash_history
cd db
git clone git@github.com:haematologic/cbioinf15.git .
git clone git@github.com:haematologic/cbioinf15.git
ls
cd cbioinf15/
ls
mv ../* .
ls
git status
git add .
ls -latr
cd ..
ls
cd ..
cd Course_Materials/
ls
cd Day1/
ls
cd nki/
ls
cd ..
ls
cd qa/
ls
ls -latr
cp *.html ../../../db/cbioinf15/
ls
ls -latr
cd ..
ls
cd nki
ls
ls -latr
cd ..
ls
cd ..
ls
cd Course_Materials/Day1/alignment-demo/
ls -latr
cp *.html ../../../db/cbioinf15/
cp SRR1186252_trimmed.fq.chr6.fq.sorted.bam.tdf ../../../db/cbioinf15/
cd ../../../
ls
cd db
ls
cd cbioinf15/
ls
git status
git add .
git commit -am "day 1 working materials"
git config user.email "e@e.com"
git config user.name "My name"
git config
git info
git status
git commit -am "day 1 working materials"
git push origin master
git remove 150727.RData
git rm 150727.RData 
git status
git commit -am "day 1 working materials"
git push origin master
git rm 150727.RData 
ls
git status
sudo rm -R .git/
rm .git/
cd .git/
ls
cd ..
rm -R .git/
y
ls
ls -latr
cd ..
git clone git@github.com:haematologic/cbioinf15.git
git clone git@github.com:haematologic/cbioinf15.git ./cbioinf15/
mv cbioinf15/* .
git clone git@github.com:haematologic/cbioinf15.git
mv * cbioinf15/
cd cbioinf15/
ls -latr
git status
git add .
git commit -am "day 1 working materials"
git config user.email "e@e.com"
git config user.name "My name"
git commit -am "browser history"
git push origin master
firefox SRR1186252_trimmed.fq.chr6.fq.sorted.bam.samstat.html 
git add .
del *
rm *
ls
ls -latr
cd .git
yes 
yes rm *
yes | rm *
ls 
yes | rm -R *
ls
ls -latr
cd ..
yes | rm -R *
ls
ls -latr
rm .git/
rmdir .git/
ls
cd ..
rmdir cbioinf15/
ls
ls -latr
rm ~/.ssh/id_rsa*.*
ls ~/.ssh/
exit
ls -latr
ls -latrls
mkdir db
cd db
ssh-keygen -t rsa -b 4096 -C "e@e.com"
cat ../.ssh/id_rsa.pub
git config user.name "DB"
git clone git@github.com:haematologic/cbioinf15.git
ls
cat ../.ssh/id_rsa.pub
ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAACAQDb74BrEMiJu+upyjsA8oiVrXGrufjaPm0WWBze9DmvAgtyTitf28Nwx2WS+AXadr6tlPjHmsGy7GuJrEjJkvCj/ZEIe16nqLo+kY4942UGo6q6Z+8TSp500cgoprVExuMuk22KcJq0KW8iuybvUp7VmHU6qJNkaBIS5WxLTynPNMVrcxTC9ENKSL1BIR72hkhosY7CJITRzdALvsDkypBfIbb5KyVFSe/Ml+Lea+dtML1tz5n/kWLX+MLXM0xoK5d/gkxGtS6zwInDUlgtzCtU96WzGzVdv6BKhNMCu6s4SNkljVcQJxlKnRoem9skGw/Ivi38ELTFJxC2Eu/JV4tGcbmn5pKPBwXiUUvyPuOri8hQuLcWRlhwpDUFdDRJ8Ct6VtvDZ72MYcvEi4L1WuszySoCM+tEOKklTYpZh7jbBAlnIpMHll0CdQEOKtkx/zQ9x/b3/DtODhWbt5yddHX6F+SSWteL5leaKSxDAIPSoWQMU1ZywecHGucQUGAMaDysXP/xSxWkawZ+GXMu1fLJBfucSPQaHI7QJSdO0ZeoDKpV6aMAmC6GoC8sgLvFBecN04j0gVB7KQEi2bLXpSvmHjQwmIAKlPDWZLK/RnMF0hoR3KCi1ek9utSerhFumJacdJMXyK+zSexabNxXYuVaS50M3jkqBGO0YiT2/ftzIQ== e@e.com
git clone git@github.com:haematologic/cbioinf15.git
cd db
ls
git clone git@github.com:haematologic/cbioinf15.git
ssh-add ~/.ssh/id_rsa
git clone git@github.com:haematologic/cbioinf15.git
ls
cd cbioinf15/
git config user.name "DB"
git config user.email "haematologic+github@gmail.com"
git config user.name "DLB"
ls
git status
git add .
git commit -am "morning practical 150728 - alignment in R of genomic data"
git push origin master
cd ../..
cd Course_Materials/Day2/
ls
ls -latr
fastqc 16N_reads.fq 
firefox 16N_reads_fastqc.html 
export BOWTIE2_INDEXES="$HOME/Course_Materials/ref_data/bowtie/"
echo $BOWTIE2_INDEXES
ls
ls ../ref_data/bowtie/
cp ../ref_data/chr22.fa ../ref_data/bowtie/
ls
ls ../ref_data/
ls ../ref_data/bowtie/
export BOWTIE2_INDEXES="$HOME/Course_Materials/ref_data/bowtie/"
ls
ls -l
ls ..
ls -l ..
ls -l .
ls .
ls ..
ls /
ls ../../
ls ~/.
mkdir Tophap2_Alignment
tophat2 --output-dir Tophat2_Alignment chr22 16N_reads.fq
samtools index Tophat2_Alignment/accepted_hits.bam
samtools flagstat Tophat2_Alignment/accepted_hits.bam
samtools flagstat Tophat2_Alignment/unmapped.bam
ls
ls -latr
mkdir Tophat2_Alignment
rm Tophap2_Alignment/
rmdir Tophap2_Alignment/
cd Tophat2_Alignment/
ls
ls -latr
cd ..
ls bam/
samtools flagstat bam/16N_aligned.bam
samtools view -H bam/16N_aligned.bam
cufflinks -h
cufflinks Tophat2_Alignment/accepted_hits.bam > Tophat2_Alignment/cufflinksResult
ls Tophat2_Alignment/
vi transcripts.gtf 
vi genes.fpkm_tracking 
vi isoforms.fpkm_tracking 
library (Rsubread)
ls
ls -latr
cd ..
ls
cd ..
ls
cd db/
ls
cd cbioinf15/
ls
ls -latr
git status
git add 150728*
git status
git commit -m "end of day 2"
git push origin master
git add 150728*
git commit -m "end of day 2"
git push origin master
ssh-remove ~/.ssh/id_rsa
ssh-add --h
ssh-add -d ~/.ssh/id_rsa
cd ~/.ssh/
ls
rm id_rsa.pub 
rm id_rsa 
cd..
ls
cd ..
ls
cd db
exit
