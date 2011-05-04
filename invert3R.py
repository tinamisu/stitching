#!/bin/env python
### Tina
### 04.20.11
import sys, re, os
sys.path.append('/Genomics/grid/users/tthu/tinaLocal/lib/python')
sys.path.append('/Genomics/grid/users/tthu/tinaLocal/src/galaxy-dist/lib')
sys.path.append('/Genomics/grid/users/tthu/tinaLocal/lib/python2.6/site-packages')
from Bio import Seq
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

dmel = SeqIO.to_dict(SeqIO.parse(open('/Genomics/grid/users/afaigon/refSeqs/dmel-all-chromosome-r5.33.fasta','rb'),'fasta'))['3R'].seq
inversion_start = 3874901 - 1
inversion_end   = 17560829
dmel = dmel[0:inversion_start] + dmel[inversion_start:inversion_end].reverse_complement() + dmel[inversion_end::]
SeqIO.write(SeqRecord(dmel,id="3Rinverted"), open('3R_inverted.fa', 'w'), "fasta")
