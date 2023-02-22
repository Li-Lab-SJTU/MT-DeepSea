#!/usr/bin/env python
from sys import argv
sc,count_file,tpm_file=argv
fin1=open(count_file)
fout=open(tpm_file,'w')
read_length=150
gene_length={}
id_read={}
gene_read={}

for line in fin1:
    if line[0]!='#' and line[0]!='G':
        temp=line.strip().split()
        name=temp[1]+'_'+temp[0].split('_')[1]
        gene_length[name]=int(temp[5])
        gene_read[name]=int(temp[6])

transcript=0.00
for i in gene_length.keys():
    transcript+=float(gene_read[i])/gene_length[i]

for i in gene_length.keys():
    tpm=(gene_read[i]/(gene_length[i]*transcript))*1000000
    fout.write(i+'\t'+str(tpm)+'\n')
