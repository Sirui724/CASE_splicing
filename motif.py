import math

f1=open('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\short_cancer_up300.fa')
f2=open('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\short_other_up300.fa')
fw=open('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\up300_sc_vs_so_frequency_6.txt','w')
fw2=open('D:\\wanglab\\smallexon\\new_60\\01features\\011motif\\zscore\\up300_sc_vs_so_frequency_6.fa','w')

nmer=6
N1=0
N2=0
GFF1={}
GFF2={}
for line in f1.readlines():
    if line[0]!='>':
        sequence=line.strip().upper()
        l=len(sequence)
        N1=N1+l-nmer+1
        for i in range(l-nmer+1):
            temp=sequence[i:i+nmer]
            if temp in GFF1:
                GFF1[temp]=GFF1[temp]+1
            else:
                GFF1[temp]=1

for line in f2.readlines():
    if line[0]!='>':
        sequence=line.strip().upper()
        l=len(sequence)
        N2=N2+l-nmer+1
        for i in range(l-nmer+1):
            temp=sequence[i:i+nmer]
            if temp in GFF2:
                GFF2[temp]=GFF2[temp]+1
            else:
                GFF2[temp]=1

fw.write('sequence\tf1\tf2\tp\tzscore\n')
a=1
for i in list(set(GFF1.keys()).intersection(set(GFF2.keys()))):
    f1=GFF1[i]/N1
    f2=GFF2[i]/N2
    p=(N1*f1+N2*f2)/(N1+N2)
    z=(f1-f2)/math.sqrt((1/N1+1/N2)*p*(1-p))
    # fw.write(i+'\t'+str(GFF1[i]/N1)+'\t'+str(GFF2[i]/N2)+'\n')
    fw.write(i + '\t' + str(f1) + '\t' + str(f2) + '\t' + str(p) +'\t'+str(z)+ '\n')
    if z >=4:
        fw2.write('>'+str(a)+'\n'+i+'\n')
        a=a+1


print('N1 is '+str(N1))
print('N2 is '+str(N2))
print(len(GFF1))
print(len(GFF2))
print(len(list(set(GFF1.keys()).intersection(set(GFF2.keys())))))
