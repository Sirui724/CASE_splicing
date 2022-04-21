from __future__ import division
import argparse


f_AS=open('/picb/rnasys/zhangsirui/TCGA/test/junction/se_ann.txt')
f_junction=open('/picb/rnasys/zhangsirui/TCGA/test/junction/gdac.broadinstitute.org_LUSC.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__junction_quantification__data.Level_3.2016012800.0.0/LUSC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__junction_quantification__data.data.txt')

GFF={}
for line in f_junction.readlines():
    line_=line.strip().split('\t')
    if line_[0]=='Hybridization REF':
        name_list=line_[1:]
    elif line_[0]=='junction':
        a=0
    else:
        GFF[line_[0]]=line_[1:]

# print(GFF['chr1:178421787:+,chr1:178423582:+'])
# print(name_list)

def SE_psi(event):
    if event[-1]=='+':
        J1 = event.split(':')[0] + ':' + event.split(':')[2] +':'+event[-1]+',' +event.split(':')[0] +':'+ event.split(':')[4] + ':' + event[-1]
        J2 = event.split(':')[0] + ':' + event.split(':')[5] + ':'+event[-1]+',' +event.split(':')[0]+':' + event.split(':')[7] + ':' + event[-1]
        J3 = event.split(':')[0] + ':' + event.split(':')[2] + ':'+event[-1]+',' +event.split(':')[0] +':'+ event.split(':')[7] + ':' + event[-1]
    else:
        J1=event.split(':')[0] + ':' + event.split(':')[8] + ':'+event[-1]+',' + event.split(':')[0] +':'+event.split(':')[4] + ':' + event[-1]
        J2= event.split(':')[0] + ':' + event.split(':')[5] +':'+event[-1]+',' +event.split(':')[0] +':'+ event.split(':')[1] + ':' + event[-1]
        J3= event.split(':')[0] + ':' + event.split(':')[8] + ':'+event[-1]+',' +event.split(':')[0]+':' + event.split(':')[1] + ':' + event[-1]
    j1=j2=j3=['NA']*len(name_list)
    if GFF.has_key(J1) == True:
        j1 = GFF[J1]
    if GFF.has_key(J2) == True:
        j2 = GFF[J2]
    if GFF.has_key(J3) == True:
        j3 = GFF[J3]
    psi=[]
    # print(event)
    # print(J1)
    # print(J2)
    # print(J3)
    for i in range(len(name_list)):
        if j1[i]!='NA' and j2[i]!='NA' and j3[i]!='NA':
            if (int(j1[i]) + int(j2[i]))>=3 and int(j3[i])>=3:
                psi.append(str(((int(j1[i]) + int(j2[i])) / 2) / ((int(j1[i]) + int(j2[i])) / 2 + int(j3[i]))))
            else:
                psi.append('NA')
        else:
            psi.append('NA')
    # list = [str(psi), str(j1), str(j2), str(j3),'']
    return psi

def A5SS_psi(event):
    if event[-1]=='+':
        J1=event.split(':')[0]+':'+event.split(':')[2].split('|')[1]+ ':'+event[-1]+',' + event.split(':')[0] +':'+event.split(':')[4]+':'+event[-1]
        J2=event.split(':')[0]+':'+event.split(':')[2].split('|')[0]+ ':'+event[-1]+',' + event.split(':')[0] +':'+event.split(':')[4]+':'+event[-1]
    else:
        J1=event.split(':')[0] + ':' +event.split(':')[5]  + ':'+event[-1]+',' + event.split(':')[0] +':' + event.split(':')[2].split('|')[0] + ':' + event[-1]
        J2=event.split(':')[0] + ':' + event.split(':')[5] + ':'+event[-1]+',' + event.split(':')[0] +':' + event.split(':')[2].split('|')[1] + ':' + event[-1]
    j1=j2=['NA']*len(name_list)
    if GFF.has_key(J1) == True:
        j1 = GFF[J1]
    if GFF.has_key(J2) == True:
        j2 = GFF[J2]
    psi=[]
    # print(J1)
    # print(J2)
    for i in range(len(name_list)):
        if j1[i]!='NA' and j2[i]!='NA':
            if int(j1[i]) + int(j2[i]) >= 5:
                psi.append(str(int(j1[i])/(int(j1[i])+int(j2[i]))))
            else:
                psi.append('NA')
        else:
            psi.append('NA')

    return psi


fw=open('/picb/rnasys/zhangsirui/TCGA/test/junction/tt.txt','w')
fw.write('AS'+'\t'+'name'+'\t'+'gene_name'+'\t'+'exon_legth'+'\t'+'intron1_length'+'\t'+'intron2_length'+'\t'+"\t".join(name_list)+'\n')
for line in f_AS.readlines():
    line_=line.strip().split('\t')
    if line_[0]=='SE':
        psi_value=SE_psi(line_[1])
        # print(psi_value)
        if psi_value != ['NA']*len(name_list):
            fw.write(line.strip()+'\t'+"\t".join(psi_value)+'\n')
    else:
        psi_value = A5SS_psi(line_[1])
        fw.write(line.strip() + '\t' + "\t".join(psi_value) + '\n')
        # print(psi_value)

