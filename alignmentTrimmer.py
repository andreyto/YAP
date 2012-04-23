########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2012 J.Craig Venter Institute.
########################################################################################


#!/usr/bin/python
#################################################
## 	A new program
#################################################
import sys, time, random, hashlib, shlex, os, itertools, gc
from optparse import OptionParser
from Bio import SeqIO, Seq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

_author="Sebastian Szpakowski"
_date="2011/01/01"
_version="Version 1"

#################################################
##		Classes
##

#################################################
##		Functions
##
def	cleanup(recordlist, start, end):
	otpt=list()
	for rec in recordlist:
		seq = str(rec.seq)[start:end]
		tmp = SeqRecord(seq=Seq(seq), id=rec.id, name=rec.id, description=rec.id)
		otpt.append(tmp)
	return (otpt)

#################################################
##		Arguments
##
parser = OptionParser()

#parser.add_option("-f", "--file", dest="filename",
#                  help="write report to FILE", metavar="FILE")
#parser.add_option("-q", "--quiet",
#                  action="store_false", dest="verbose", default=True,
#                  help="don't print status messages to stdout")

parser.add_option("-I", "--input", dest="inputname",
                  help="load fasta FILE", metavar="FILE")
                  
parser.add_option("-s", "--start", type="int", dest="start",
                  help="trim start position [default: %default]", metavar="N", default=10)

parser.add_option("-e", "--end", type="int", dest="end",
                  help="trim end position [default: %default]", metavar="N", default=20)

(options, args) = parser.parse_args()

#################################################
##		Begin
##

inputfilename = options.inputname
extension = inputfilename.strip().split(".")[-1]

outputfilename ="%s.%s_%s.%s" % (".".join(inputfilename.strip().split(".")[:-1]), options.start, options.end, extension)
output_handle = open(outputfilename, "w")

tmp = list(SeqIO.parse(open(options.inputname), "fasta"))
SeqIO.write(cleanup(tmp, options.start, options.end), output_handle, "fasta")
output_handle.close()



#################################################
##		Finish
#################################################
