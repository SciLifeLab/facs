from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement, transcribe, back_transcribe, translate
import re
import time

A_t = 0;
G_t = 0;
C_t = 0;
T_t = 0;
All = 0;
All_r = 0;

'''
x1 = open("test.fna")
for seq in SeqIO.parse(x1,'fasta'):
    All_r+=1;
    for word in seq:
        All+=1;
        if word=='A':
            A_t+=1;
        elif word=='C':
            C_t+=1;
        elif word=='G':
            G_t+=1;
        elif word=='T':
            T_t+=1;

print 'A_t->',A_t,'G_t->',G_t,'C_t->',C_t,'T_t->',T_t;
print 'ALL->',All;
print 'average length->', All//All_r;
'''
signal = 0;
New_string = '';
x1 = open("10xtest_GC44.fasta")
for seq in SeqIO.parse(x1,'fasta'):
    if signal == 0:
       New_string+='>'
       New_string+=seq.id;
       New_string+='\n'
       New_string+=seq.seq.tostring();
       New_string+='\n';
       signal = 1;
    elif signal == 1:
        signal = 0;
        
x2 = open("single.fasta",'w')
x2.write(New_string)
x2.close()
