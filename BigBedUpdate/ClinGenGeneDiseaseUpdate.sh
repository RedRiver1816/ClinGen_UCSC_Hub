#!/bin/bash
#A script for updating the BigBed Files in the ClinGen UCSC Track Hub

cd "$(dirname "$0")"

#fetch chrom sizes for each distribution
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes
chmod a+x ./fetchChromSizes
./fetchChromSizes hg19 > hg19.chrom.sizes
./fetchChromSizes hg38 > hg38.chrom.sizes

#extract ClinGen raw data and create two Bed files
printf "INFO: Downloading ClinGen Data and Converting to BED"
python ./UCSC-ClinGenBedCreator.py

#sort Bed files in prep for bedToBigBed
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedSort
chmod a+x ./bedSort
./bedSort hg19ClinGenBigBed.bed hg19ClinGenBigBed.bed
./bedSort hg38ClinGenBigBed.bed hg38ClinGenBigBed.bed

#converts Bed files into BigBed ready for upload to Github
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
chmod a+x ./bedToBigBed
./bedToBigBed hg19ClinGenBigBed.bed hg19.chrom.sizes hg19final.bb -type=bed9+9 -as=./hg19as.as -tab -extraIndex=name,geneSymbol,HGNCid,MONDOid,Classification
./bedToBigBed hg38ClinGenBigBed.bed hg38.chrom.sizes hg38final.bb -type=bed9+9 -as=./hg38as.as -tab -extraIndex=name,geneSymbol,HGNCid,MONDOid,Classification

#removes unneeded files
rm fetchChromSizes bedSort bedToBigBed hg19.chrom.sizes hg38.chrom.sizes hg19ClinGenBigBed.bed hg38ClinGenBigBed.bed