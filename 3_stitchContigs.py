#!/bin/env python
### Tina
### 04.20.11
import sys, re, os
import itertools
sys.path.append('/Genomics/grid/users/tthu/tinaLocal/lib/python')
sys.path.append('/Genomics/grid/users/tthu/tinaLocal/src/galaxy-dist/lib')
sys.path.append('/Genomics/grid/users/tthu/tinaLocal/lib/python2.6/site-packages')
import killableprocess
from cmdline.cmdline import CommandLineApp
from Bio import Seq
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

indel_search = re.compile("-+")

class App(CommandLineApp):
	def __init__(self):
		CommandLineApp.__init__(self)

		op = self.option_parser
		op.set_usage('usage: stitchContigs.py')
		op.add_option('-c', '--chrom', dest='chrom', type='string', default='X', help='chromosome')

	def main(self):

      ####################################################################################################
      ### perform full stitching of these contigs
      # NODE_100094_length_35241_cov_17.552652  X       8389940 8721794 197     35303   +       35303   197     0.00558026230065434
		full_stitch_contigs = list()
		infile = open('unique_contigs_full.tsv','r')
		for line in infile:
		   data = line.strip().split()
		   full_stitch_contigs.append(data[0])

		infile.close()

		contigs = SeqIO.to_dict(SeqIO.parse(open('contigs.fa','rb'),'fasta'))
		dmel = SeqIO.to_dict(SeqIO.parse(open('dmel-all-chromosome-r5.33_ARMS_3Rinv.fasta','rb'),'fasta'))[self.options.chrom].seq

		contig_total = 0
		dmel_total = 0

		last_pos = 0
		new_dsim = ''
		new_dsim_status = ''

		junctions4dsim = list()
		infile = open("chr%s_mapping-1000-500.tsv" % self.options.chrom,'rb')
		#infile = open("dmel_chr%s_mapping.tsv" % self.options.chrom,'rb')
		for line in infile:
         # 2L         29127          263223        NODE_104920_length_230114_cov_19.991587    554  230106  +
         # 2L         47591           49430        NODE_11013_length_1777_cov_64.561058         0    1839  -
         # 2L         52179           52751        NODE_74998_length_4110_cov_8.511922       1149    1711  +
		   dmel_refChr, dmel_refStart, dmel_refEnd, contig, contig_start, contig_end, strand = line.strip().split()
		   dmel_refStart, dmel_refEnd, contig_start, contig_end = int(dmel_refStart), int(dmel_refEnd), int(contig_start), int(contig_end)

		   junctions4dsim.append("%s\t%d\t%d\t%d\t%d\n" % (self.options.chrom,last_pos,dmel_refStart-1,len(new_dsim)-1,len(new_dsim)-1+dmel_refStart-last_pos))

         ### just append, don't add in dmel
		   if int(dmel_refStart) < last_pos: 
		      print "Overlaps previous hit wrt dmel by %d (previous pos %d)" % (dmel_refStart-last_pos,last_pos)
		      print line
		      #sys.exit()
		   else: 
		      new_dsim += str(dmel[last_pos:dmel_refStart])
		      new_dsim_status += '!' * (dmel_refStart-last_pos)
		      dmel_total += dmel_refStart - last_pos

		   if contig in full_stitch_contigs: append_dsim = formatSeq(contigs[contig].seq,strand)
		   else: append_dsim = formatSeq(contigs[contig].seq[contig_start:contig_end],strand)
		   new_dsim += append_dsim
		   new_dsim_status += '"' * len(append_dsim)
		   contig_total += len(append_dsim)

		   last_pos = dmel_refEnd + 1

		infile.close()

      ### tack on the rest of dmel
		new_dsim += str(dmel[last_pos:])
		new_dsim_status += '!' * (len(dmel)-last_pos)
		dmel_total += len(dmel)-last_pos

		print "Final length %d" % len(new_dsim)
		print "Added dmel   %d" % dmel_total
		print "From velvet  %d" % contig_total

      ### WRITE OUT new dsim ref
      ####################################################################################################
		out = open('dsim_new.%s.fa' % self.options.chrom, 'wb')
		out.write(">%s\n%s\n" % (self.options.chrom, new_dsim))
		out.close()

		out = open('dsim_new.%s.qual' % self.options.chrom, 'wb')
		out.write(">%s\n%s\n" % (self.options.chrom, new_dsim_status))
		out.close()

      ### Report junctions that need to be filled in with dsim
      ####################################################################################################
		out = open('dsim_new.%s.junctions' % self.options.chrom, 'wb')
		out.writelines(junctions4dsim)
		out.close()



####################################################################################################
####################################################################################################
####################################################################################################

def formatSeq(seq,strand):
    if (strand == '+'): return str(seq)
    else: return str(seq.reverse_complement())

def returnBase(base,flag):
   if flag == True: return base
   else: return ''


if __name__ == '__main__':
    try:
        App().run()
    except Exception, e:
        print 'ERROR in stitchContigs.py\n'
        print '%s' % e
        sys.exit(2)
