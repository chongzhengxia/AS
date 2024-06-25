# AS
1 The human reference genome was downloaded from ucsc GRCh38 #https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

2 The genome annotation file was downloaded from ucsc refGene.gtf #https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz

3 AS.sh : This script uses the SLURM job system to submit jobs. FastQ raw data is used as input, and RMATS variable spliced data is used as output

4 merge.R : The R script is used to merge the 5 different splice events output from the RMATS into a PSI matrix. If the NA value of a splicing event is less than 20% in all samples, the event will be marked as pass, otherwise it will be marked as fail
