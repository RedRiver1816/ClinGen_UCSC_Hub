#!/usr/bin/env python
# coding: utf-8


#import necessary libraries

import numpy as np
import pandas as pd
import requests


#read in raw data from ClinGen

clingen = pd.read_csv('https://search.clinicalgenome.org/kb/gene-validity.csv', sep = ',', header = 4, skiprows = [5])

#assign necessary new columns and rearrange

newcols = ['DISEASE LABEL','GENE SYMBOL','GENE ID (HGNC)','DISEASE ID (MONDO)','MOI','SOP','CLASSIFICATION','ONLINE REPORT','CLASSIFICATION DATE']
clingen = clingen[newcols]
clingen.insert(1,'SCORE',1000)
clingen.insert(2,'ITEM','.')
clingen.insert(3,'RGB',0)
clingen['MOUSEOVER'] = clingen['GENE SYMBOL'] + ' - ' + clingen['DISEASE LABEL'] + ' - ' + clingen['MOI'] + ' (' + clingen['CLASSIFICATION'] + ')'
clingen['ONLINE REPORT'] = clingen['ONLINE REPORT'].str[51:]
clingen['CLASSIFICATION DATE'] = clingen['CLASSIFICATION DATE'].str[:10]
RGB = {'Definitive':'120,0,170','Strong':'0,0,102','Moderate':'0,0,255','Limited':'153,153,255','Disputed':'255,153,51','Refuted':'255,51,51','No Reported Evidence':'0,0,0'}
clingen['RGB'] = clingen.apply(lambda x: RGB[x['CLASSIFICATION']], axis=1)

#retrieve raw chromosome data for hg19 genes

hg19raw = requests.get('http://api.genome.ucsc.edu/getData/track?genome=hg19;track=ncbiRefSeqCurated;maxItemsOutput=-1')
hg19rawdata = hg19raw.json()

#extract relevant chromosomal data from raw json

transdata = pd.DataFrame(index=range(5000),columns = ['chrom','txStart','txEnd','GENE SYMBOL'])
row = 0

for chrom in hg19rawdata['ncbiRefSeqCurated']:
    if '_' not in chrom:
        for gene in hg19rawdata['ncbiRefSeqCurated'][chrom]:
            if gene['name2'] in clingen['GENE SYMBOL'].unique():
                transdata.iloc[row] = [gene['chrom'],gene['txStart'],gene['txEnd'],gene['name2']]
                row += 1
            
#filter data for longest unique transcripts for each gene            

transdata.dropna(inplace = True)
transdata['len'] = transdata['txEnd'] - transdata['txStart']
transdata.sort_values('len',ascending = False,inplace = True)
transdata.drop_duplicates('GENE SYMBOL',keep = 'first',inplace = True)

#merge chromosome info and clingen data and rearrange to fit UCSC specs

hg19final = pd.merge(clingen,transdata,on = 'GENE SYMBOL')
hg19final.drop('len',axis = 1,inplace = True)
newcols2 = ['chrom', 'txStart', 'txEnd', 'DISEASE LABEL', 'SCORE', 'ITEM', 'RGB', 'GENE SYMBOL', 'GENE ID (HGNC)', 'DISEASE ID (MONDO)', 'MOI', 'SOP', 'CLASSIFICATION', 'ONLINE REPORT', 'CLASSIFICATION DATE', 'MOUSEOVER']
hg19final = hg19final[newcols2]
hg19final.insert(6,'ExonStart',hg19final['txStart'])
hg19final.insert(7,'ExonEnd',hg19final['txEnd'])

#send final data to BED file

hg19final.to_csv('hg19ClinGenBigBed.bed',sep = '\t',header = False,index = False)

#same as above, except extracts hg38 genome data

hg38raw = requests.get('http://api.genome.ucsc.edu/getData/track?genome=hg38;track=ncbiRefSeqCurated;maxItemsOutput=-1')
hg38rawdata = hg38raw.json()

#extract relevant chromosomal data from raw json

transdata = pd.DataFrame(index=range(5000),columns = ['chrom','txStart','txEnd','GENE SYMBOL'])
row = 0

for chrom in hg38rawdata['ncbiRefSeqCurated']:
    if '_' not in chrom:
        for gene in hg38rawdata['ncbiRefSeqCurated'][chrom]:
            if gene['name2'] in clingen['GENE SYMBOL'].unique():
                transdata.iloc[row] = [gene['chrom'],gene['txStart'],gene['txEnd'],gene['name2']]
                row += 1
            
#filter data for longest unique transcripts for each gene            

transdata.dropna(inplace = True)
transdata['len'] = transdata['txEnd'] - transdata['txStart']
transdata.sort_values('len',ascending = False,inplace = True)
transdata.drop_duplicates('GENE SYMBOL',keep = 'first',inplace = True)

#merge chromosome info and clingen data and rearrange to fit UCSC specs

hg38final = pd.merge(clingen,transdata,on = 'GENE SYMBOL')
hg38final.drop('len',axis = 1,inplace = True)
newcols2 = ['chrom', 'txStart', 'txEnd', 'DISEASE LABEL', 'SCORE', 'ITEM', 'RGB', 'GENE SYMBOL', 'GENE ID (HGNC)', 'DISEASE ID (MONDO)', 'MOI', 'SOP', 'CLASSIFICATION', 'ONLINE REPORT', 'CLASSIFICATION DATE', 'MOUSEOVER']
hg38final = hg38final[newcols2]
hg38final.insert(6,'ExonStart',hg38final['txStart'])
hg38final.insert(7,'ExonEnd',hg38final['txEnd'])

#send final data to BED file

hg38final.to_csv('hg38ClinGenBigBed.bed',sep = '\t',header = False,index = False)