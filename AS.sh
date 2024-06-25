##step1 Raw data is processed using fastp
------------------------------
#!/bin/bash
#SBATCH -p debug
#SBATCH -J fastp
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o fastp.o
#SBATCH -e fastp.e
cd /PUBLIC/home/xiachongzheng/AS/01_rawdata/data
for i in $(ls -d  H9*)
do
        mkdir /PUBLIC/home/xiachongzheng/AS/02_cleandata/${i}
        cd /PUBLIC/home/xiachongzheng/AS/01_rawdata/data/${i}
        fastp -w 16 -i ${i}_1.fq.gz -I ${i}_2.fq.gz -o /PUBLIC/home/xiachongzheng/AS/02_cleandata/${i}/${i}_1.fq.gz -O /PUBLIC/home/xiachongzheng/AS/02_cleandata/${i}/${i}_2.fq.gz
done

------------------------------

##step2 Use STAR to align QC data to a reference genome
------------------------------
#!/bin/bash
#SBATCH -p common
#SBATCH -J star
#SBATCH -N 2
#SBATCH -n 104
#SBATCH -o star.o
#SBATCH -e star.e
cd /PUBLIC/home/xiachongzheng/AS/02_cleandata
for i in $(ls -d H9*)
do
        cd /PUBLIC/home/xiachongzheng/AS/02_cleandata/${i}
        STAR \
        --runMode alignReads \
        --readFilesCommand zcat \
        --runThreadN 104 \
        --genomeDir /PUBLIC/home/xiachongzheng/AS/00_ref/ucsc/STAR_hg38_index \
        --readFilesIn ${i}_1.fq.gz ${i}_2.fq.gz \
        --outFileNamePrefix /PUBLIC/home/xiachongzheng/AS/03_bamSTAR/${i} \
        --outSAMtype BAM Unsorted \
        --chimSegmentMin 2 --outFilterMismatchNmax 3 --outSAMstrandField intronMotif --alignEndsType EndToEnd --alignSJDBoverhangMin 6 --alignIntronMax 299999
-----------------------

##step3 Use SamTools to sort the bam files
-----------------------
#!/bin/bash
#SBATCH -p common
#SBATCH -J samtools
#SBATCH -N 2
#SBATCH -n 104
#SBATCH -o samtools.o
#SBATCH -e samtools.e
cd /PUBLIC/home/xiachongzheng/AS/03_bamSTAR
ls -al | awk '$9~/.bam$/{print $9}'  | while read filename
do
        echo "processing ${filename}"
        samtools sort -@104 ${filename} -o /PUBLIC/home/xiachongzheng/AS/03_bamSTAR/bam/${filename%.bam}.sort.bam
        echo "done"
done 
-----------------------

##step4 Use RMATS to detect alternative splicing events
-----------------------
#!/bin/bash
#SBATCH -p common
#SBATCH -J rmats
#SBATCH -N 2
#SBATCH -n 104
#SBATCH -o rmats.o
#SBATCH -e rmats.e
cd /PUBLIC/home/xiachongzheng/AS/03_bamSTAR/bam
filenames=($(ls -al | awk '$9~/^H9/{print $9}'))
let times=${#filenames[@]}/2
for((i=0,j=0;i<times;i++,j+=2))
do
        echo "processing ${filenames[j]%Aligned.out.sort.bam}"
        ls $(pwd)/$(echo ${filenames[j]}) > b1.txt
        ls $(pwd)/$(echo ${filenames[j+1]}) > b2.txt
        mkdir /PUBLIC/home/xiachongzheng/AS/04_rmats/${filenames[j]%Aligned.out.sort.bam}
        rmats.py --b1 b1.txt --b2 b2.txt --gtf /PUBLIC/home/xiachongzheng/AS/00_ref/ucsc/hg38.refGene.gtf  \
                 --tmp /PUBLIC/home/xiachongzheng/AS/04_rmats/${filenames[j]%Aligned.out.sort.bam}/tmp \
                 -t paired \
                 --nthread 104 --od /PUBLIC/home/xiachongzheng/AS/04_rmats/${filenames[j]%Aligned.out.sort.bam} --readLength 150 --statoff
        rm b1.txt b2.txt
        echo "done"
done 
-----------------------






















