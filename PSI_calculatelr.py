from __future__ import division
import argparse
import pysam
import os

def createHelp():
    epilog_string="---zsr"
    description_string = 'The program is going to calculate psi-value based on junctions'
    parser = argparse.ArgumentParser(description=description_string, epilog=epilog_string)
    parser.add_argument('-i', '--file-name', dest='fnIn', help='input files name(junction)')
    op = parser.parse_args()
    return op

op = createHelp()

f_AS=open('/picb/rnasys/zhangsirui/TCGA/test/junction/se_ann.txt')
#f_junction=open('/picb/rnasys/zhangsirui/smallexon/long-read/reference_file/minimap2/out/junc_2/HMEC.junc.txt')
a=op.fnIn
HMECTabix = pysam.Tabixfile('/picb/rnasys/zhangsirui/smallexon/long-read/reference_file/minimap2/out/junc_2/'+a+'.junc_sort.txt.gz','r')


# print(name_list)
def get_junc(junc_ann,junc_seq):
    junc1=[int(junc_ann.split(':')[1].split('-')[0]),int(junc_ann.split(':')[1].split('-')[1])]
    junc2=[int(junc_seq[1]),int(junc_seq[2])]
    if len(list(set(range(junc1[1]-6,junc1[1]+6))&set(range(junc2[1]-6,junc2[1]+6))))==0:
        return False
    else:
        return True


def SE_psi(event):
    if event[-1]=='+':
        J1 = event.split(':')[0] + ':' + event.split(':')[2] +'-'+ event.split(':')[4]
        J2 = event.split(':')[0] + ':' + event.split(':')[5] +'-' + event.split(':')[7]
        J3 = event.split(':')[0] + ':' + event.split(':')[2] + '-'+ event.split(':')[7]
    else:
        J1=event.split(':')[0] + ':' + event.split(':')[8] +'-'+event.split(':')[4]
        J2= event.split(':')[0] + ':' + event.split(':')[5] +'-'+ event.split(':')[1]
        J3= event.split(':')[0] + ':' + event.split(':')[8] +'-' + event.split(':')[1]
#chr3:152163328:152165409
    j1=j2=j3=0
    for target_site in HMECTabix.fetch(J1.split(':')[0], int(J1.split(':')[1].split('-')[0]) - 6, int(J1.split(':')[1].split('-')[0]) + 6):
        fields_af = target_site.rstrip().split('\t')
        if get_junc(J1,fields_af)==True:
            j1=j1+int(fields_af[3])

    for target_site in HMECTabix.fetch(J2.split(':')[0], int(J2.split(':')[1].split('-')[0]) - 6, int(J2.split(':')[1].split('-')[0]) + 6):
        fields_af = target_site.rstrip().split('\t')
        if get_junc(J2,fields_af)==True:
            j2=j2+int(fields_af[3])
    for target_site in HMECTabix.fetch(J3.split(':')[0], int(J3.split(':')[1].split('-')[0]) - 6, int(J3.split(':')[1].split('-')[0]) + 6):
        fields_af = target_site.rstrip().split('\t')
        if get_junc(J3,fields_af)==True:
            j3=j3+int(fields_af[3])


    # print(event)
    # print(J1)
    # print(J2)
    # print(J3)

    if (int(j1) + int(j2) +int(j3))>=3:
        psi=str(((int(j1) + int(j2)) / 2) / ((int(j1) + int(j2)) / 2 + int(j3)))
    else:
        psi='NA'

    # list = [str(psi), str(j1), str(j2), str(j3),'']
    return psi



fw=open('/picb/rnasys/zhangsirui/smallexon/long-read/reference_file/minimap2/out/psi_out/'+a+'_psi.txt','w')
fw.write('AS'+'\t'+'name'+'\t'+'gene_name'+'\t'+'exon_legth'+'\t'+'intron1_length'+'\t'+'intron2_length'+'\t'+a+'\n')
for line in f_AS.readlines()[0:1000]:
    line_=line.strip().split('\t')
    if line_[0]=='SE' and '_' not in line:
        # print(line)
        psi_value=SE_psi(line_[1])
        # print(psi_value)
        fw.write(line.strip()+'\t'+psi_value+'\n')
    # else:
    #     psi_value = A5SS_psi(line_[1])
    #     fw.write(line.strip() + '\t' + "\t".join(psi_value) + '\n')
    #     # print(psi_value)

