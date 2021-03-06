########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2012 J.Craig Venter Institute.
########################################################################################

#################################################
## 	A pipeline for cleaning up sff files and producing 
##	OTUs 
#################################################
import sys
from optparse import OptionParser, OptionGroup
from StepsLibrary import *
from collections import defaultdict
from Queue import Queue

_author="Sebastian Szpakowski"
_date="2012/10/01"
_version="Version 2"

#################################################
##		Classes
##
	#################################################
	### Iterator over input file.
	### every line is converted into a dictionary with variables referred to by their 
	### header name
class 	GeneralPurposeParser:
	def	__init__(self, file, skip=0, sep="\t"):
		self.filename = file
		self.fp = open(self.filename, "rU")	
		self.sep = sep
		
		self.linecounter = 0
		self.currline=""
		
	def	__iter__(self):
		return (self)
	
	def	next(self):
		otpt = dict()
		for currline in self.fp:
			currline = currline.strip().split(self.sep)
			self.currline = currline
			self.linecounter = self.linecounter + 1
			return(currline)			
		raise StopIteration
					
	def	__str__(self):
		return "%s [%s]\n\t%s" % (self.filename, self.linecounter, self.currline)

class	InfoParser:
	def	__init__(self, filename):
		self.filename = filename
		self.info = GeneralPurposeParser(filename, sep=",")
		self.store = defaultdict(list)
		for line in self.info:
			if line[0].endswith("/"):
				line[0] = line[0][:-1]
			path = "%s/%s" % (line[0].strip("\""),line[1].strip("\""))
			self.store[path].append(line[2:])
		
	def	makeOligoFile(self, x):
		tmp=defaultdict(list)
		### group barcodes and primers
		for line in self.store[x]:
			barcode = line[0].strip().strip("\"")
			forward = line[1].strip().strip("\"")
			reverse = line[2].strip().strip("\"")
			needTwoPrimers = line[3].strip().strip("\"")
			sampleID =line[4].strip().strip("\"")
						
			if needTwoPrimers in ("Yes", "yes", "Y", "YES"):
				if forward == "" or reverse =="":
					print "%s: two primers needed, only one supplied F:'%s'-R:'%s' %s\n check %s " % (x, forward, reverse, self.name)
					sys.exit(1)
				bin = "B-%s-%s" % (forward, reverse)
				tmp[bin].append((barcode, sampleID))
			else:
				if reverse =="" and forward =="":
					print "%s: please provide at least one primer for barcode:'%s' " % (x, barcode)
					sys.exit(2)
				
				if forward != "":
					bin ="F-%s" % (forward)
					tmp[bin].append((barcode, sampleID))
				
				if reverse != "":
					bin ="R-%s" % (reverse)
					tmp[bin].append((barcode, sampleID))
			
		oligofiles=list()	
		for bin in tmp:
			otpt=""
			fn_otpt=""
			if bin.startswith("B"):
				fn_otpt = "%s.B.oligos" % (x.split("/")[-1])
				T, F, R = bin.strip().split("-")
				otpt += "forward\t%s\n" % (forward)
				otpt += "reverse\t%s\n" % (forward)
				otpt += "# %s both primers\n"
			else:
				D, P = bin.strip().split("-")
				if D == "F":
					fn_otpt = "%s.F.oligos" % (x.split("/")[-1])
					otpt += "forward\t%s\n" % (P)
					otpt += "# %s this IS the forward primer\n" % (x)
				else:
					fn_otpt = "%s.R.oligos" % (x.split("/")[-1])
					otpt += "forward\t%s\n" % (P)
					otpt += "# %s this is the REVERSE primer\n" % (x)
			
			for barcode, sampleID in tmp[bin]:
				otpt+= "barcode\t%s\t%s\n" % (barcode, sampleID)
				
			otptfile = open(fn_otpt, "w")
			otptfile.write(otpt)
			otptfile.close()
			oligofiles.append(fn_otpt)
		return (oligofiles)	
				
		
	def removeNonFiles(self, x):
		return x != "path/file"
	
	def getFiles(self):
		return (filter(self.removeNonFiles,  self.store.keys()))
			
#################################################
##		Functions
##
	################################################
	### Read in a file and return a list of lines
	###
def	loadLines(x):
	try:
		fp = open(x, "rU")
		cont=fp.readlines()
		fp.close()
		#print "%s line(s) loaded."  % (len(cont))
	except:
		cont=""
		#print "%s cannot be opened, does it exist? " % ( x )	
	return cont

def	getOligoFlip(filename):
	if filename.find(".R.")>-1:
		return "T"
	else:
		return "F"
	
def	deconvolute():
	forprocessing = InfoParser(options.fn_info)
	DECONVOLUTION = list()
	
	if options.strictlevel==2:
		bdiffs="0"
		pdiffs="1"
	else:
		bdiffs="0"
		pdiffs="0"
		
	for file in forprocessing.getFiles():
		oligos = forprocessing.makeOligoFile(file)
		sffinputs = {"sff": ["%s" % (file)]}
	
		print oligos
		# #### extract fasta and qual files
		
		if options.useflows:
			A = MothurStep("sffinfo", options.nodesize, sffinputs, {"flow":"T"}, list())
			
		else:
			D = SFFInfoStep(sffinputs, dict(), list())	
	
	
		for oligofile in oligos:
			oligoinputs = {"oligos": ["../%s" % (oligofile)]}
			
			if options.useflows:
				args = { 
						"bdiffs": bdiffs, 
						"pdiffs": pdiffs, 
						"maxhomop":"7",
						"maxhomop":"7", 
						"minflows":"360", 
						"maxflows":"720"
				}		
				B = MothurStep("trim.flows", options.nodesize, oligoinputs, args, [A])
				C = MothurSHHH([B])
				D = FileMerger("fasta,name,qfile", [C])
				
			#### deconvolution
			args = {	"flip":getOligoFlip(oligofile), 
						"bdiffs": bdiffs, 
						"pdiffs": pdiffs, 
						"minlength": "%s" % (options.minlength),
						"maxlength": "2500",
						"qtrim": "F", 
						"maxambig":"0", 
						"maxhomop":"7"
			}	
			if options.useflows:
				args["force"] =  "name,fasta,oligos"
			
			E = MothurStep("trim.seqs", options.nodesize, oligoinputs, args, [D])	
			
			if  options.useflows:
				DECONVOLUTION.append(E)
			
			else:
				#### generate trimming coordinates using LUCY
				F = LUCYcheck(options.nodesize, [E])

				#### trim the reads using LUCY-generated coordinates
				G = LUCYtrim([F])
						
				DECONVOLUTION.append(G)
	
#	sys.exit(0)
	ALMOST_DONE = FileMerger("fasta,name,group,qfile", DECONVOLUTION)
	READY = MatchGroupsToFasta(dict(), [ALMOST_DONE])
	
	return (READY)

def	finalize(input):
	clean = CleanFasta(dict(), [input])
	
	####### remove sequences that are too short, and with ambiguous bases 
	args = { "minlength" : "%s" % ( options.minlength ),
			 "maxambig" : "0",
			 "force": "fasta,name,group"}
	clean2 = MothurStep("screen.seqs", options.nodesize, dict(), args, [clean])

	OutputStep("PRE454", "fasta,group,name,list,svg,pdf,tiff,taxsummary,globalsummary,localsummary", clean)

	###################### CDHIT
	#### unique and de-noise
	args= 		{ 	"c" : "1.0",
					"b" : "8",
					"aS": "1.0",
					"g" : "1",
					"M"	: "50000",
					"T" : "60"
				}
	#### de-noise/unique collapse			
	CD_1 = CDHIT_454(options.nodesize, args, [clean2])
	CD_2 = CDHIT_Mothurize(dict(), CD_1)
	
	#### add ecoli reference to the merged experiments' file
	CD_3 = FileMerger("fasta,name,group,qfile", [CD_2, REF_1, REF_2, REF_3])
	
	#### align to SILVA reference
	inputs = {"reference": ["%s/silva.bacteria.fasta" % options.dir_anno] }
	args = {"flip":"t", "ksize": "9"}	
	CD_4 = MothurStep("align.seqs", options.nodesize, inputs, args, [CD_3])
	
	#### plots (inconsequential) 
	CD_5 = AlignmentSummary(dict(),[CD_4])
	CD_6 = AlignmentPlot(dict(),[CD_5])
	supplementary.append(CD_6)
	###########################
	
	OutputStep("ALIGNED", "tre,fasta,group,name,list,svg,pdf,tiff,taxsummary,globalsummary,localsummary", CD_4)	
	
	cleanCD = cleanup(CD_4)

	OutputStep("CLEAN", "fasta,group,name,list,svg,pdf,tiff,taxsummary,globalsummary,localsummary", cleanCD)
	
	clusterCD = CDHITCluster(cleanCD)
	
	x = plotsAndStats(clusterCD)
	INS = {"annotation" : [options.fn_info]}
	ARGS = {"dist": "0.03"}
	output = R_defaultplots(INS, ARGS, x)
	output2 = AnnotateClusters(dict(), dict(), output)
		
	return (output2)

def	cleanup(input):

	### remove the "ref" group
	args = {
			"force" : "fasta,name,group",
			"groups": "ref" 
			}
			
	s15 = MothurStep("remove.groups", options.nodesize, dict(), args, [input])
	
	####### remove sequences that are too short (bad alignment?)  
	args = {
				"minlength" : "%s" % (options.minlength), 
				"maxambig" : "0",
				"force" : "fasta,name,group" ,
			}
	s16 = MothurStep("screen.seqs", options.nodesize, dict(), args, [s15])
	
	####### find chimeric sequences  
	toremove = list()
	for ch in [ "uchime" ]:
		### chimeras against reference
		args = {"force" : "fasta,reference"}
		inputs = {"reference": ["%s/silva.bacteria.fasta" % options.dir_anno] }
		
		A = MothurStep("chimera.%s" % (ch),options.nodesize, inputs, args, [s16])	
		toremove.append(A)
		
		### chimeras against self
		args ={"force": "name,group,fasta"}
		inputs = {}
		
		A = MothurStep("chimera.%s" % (ch),options.nodesize, inputs, args, [s16])	
		toremove.append(A)
		
	### merge all accnos files and remove ALL chimeras	
	allchimeras = FileMerger("accnos", toremove)
	s17 = MothurStep("remove.seqs",options.nodesize, dict(), dict(), allchimeras)

	### primer cut
	args = {
				"s" : "1044", 
				"e": "13127"
			}
	s18a = AlignmentTrim(dict(), args, s17)
	
	####### remove sequence fragments, bad alignments (?) 
	args = { "minlength" : "%s" % (options.minlength),
			 "force": "fasta,name,group"}
	s18b = MothurStep("screen.seqs", options.nodesize, dict(), args, [s18a])
	
	### build a tree
	#s18b_tree = ClearcutTree({}, s18b)
	
	####### remove empty columns
	args = {"vertical" : "T"}
	s19 = MothurStep("filter.seqs",options.nodesize, dict(), args, [s18b])
	
	####### taxonomy
	inputs = {	"reference": ["%s/trainset9_032012.pds.fasta" % options.dir_anno],
				"taxonomy": ["%s/trainset9_032012.pds.tax" % options.dir_anno]
			}
			
	args = {"iters" : "100"}
	s20 = MothurStep("classify.seqs", options.nodesize, inputs, args, [s19])
	
	### remove - and . for subsequent clustering efforts 
	s21 = CleanFasta(dict(), [s20])
	
	return (s21)

def	CDHITCluster(input):
	cdhits = list()
	for arg in ["0.99", "0.97", "0.95", "0.90"]:
		args = {"c": arg,
				"n": "8",
				"g": "1",
				"M": "10000",
				"T": "15"
				}
		
		CD_1 = CDHIT_EST(options.nodesize, args, [input])
		
		### make analogous to mothur's labels
		arg = 1.0 - float(arg)
		if arg == 0:
			arg = "unique"
		else:
			arg = "%s" % (arg)
		
		args = {"mode": arg	
				}
		CD_2 = CDHIT_Mothurize(args, CD_1)
		CD_2a = CDHIT_Perls({}, CD_2)			
		cdhits.append(CD_2)
				
	READY = FileMerger("list,rabund,sabund", cdhits)	
	SORTED = FileSort("list,rabund,sabund", READY)
	return (SORTED)

def	plotsAndStats(input):
	s23 = GroupRetriever([input])
	
	######## make a shared file 
	args = {"label" : "0.01-0.03-0.05-0.1", "find": "groups"}
	s24 = MothurStep("make.shared", options.nodesize, dict(), args, [s23])


	args = {
			"label" : "0.01-0.03-0.05-0.1",
			"basis" : "otu"
			}
			
	s25a= MothurStep("classify.otu",  options.nodesize, dict(), args, [s24])
	args = {
				"taxonomy": "otu.taxonomy",
				"taxsummary": "otu.taxsummary"				
			}
	s25aa = FileType(args, [s25a])
	
	args = {
			"label" : "0.01-0.03-0.05-0.1",
			"basis" : "sequence"
			}
			
	s25b = MothurStep("classify.otu",  options.nodesize, dict(), args, [s24])
	args = {
				"taxonomy": "seq.taxonomy",
				"taxsummary": "seq.taxsummary"				
			}
	s25bb = FileType(args, [s25b])
	
	args = {"force" : "list", "calc": "nseqs-sobs-simpson-invsimpson-chao-shannon-shannoneven-coverage"}
	s26 = MothurStep("summary.single",options.nodesize, dict(), args, [s25bb])
	
	args = {"summary": "globalsummary"}
	s26a = FileType(args, [s26])
	
	args = {"force" : "shared", "calc": "nseqs-sobs-simpson-invsimpson-chao-shannon-shannoneven-coverage"}
	s27 = MothurStep("summary.single", options.nodesize, dict(), args, [s25bb])
	
	args = {"summary": "localsummary"}
	s27a = FileType(args, [s27])
	
	args = {"force" : "shared", "calc": "thetayc-jclass-braycurtis"}
	s28 = MothurStep("tree.shared", options.nodesize, dict(), args, [s24]) 
	
	supplementary.append(s28)
	
	args = {"force" : "list", "calc": "nseqs-sobs-simpson-invsimpson-chao-shannon-shannoneven-coverage"}
	s29 = MothurStep("rarefaction.single", options.nodesize, dict(), args, [s24])
	
	args = {"force" : "shared", "calc": "nseqs-sobs-simpson-invsimpson-chao-shannon-shannoneven-coverage"}
	s30 = MothurStep("rarefaction.single",options.nodesize, dict(), args, [s24])
	
	return ([s23, s24, s25aa, s25bb, s26a, s27a, s28, s29, s30])
		
#################################################
##		Arguments
##

parser = OptionParser()

group = OptionGroup(parser, "Required", description="Will not run without these !")

group.add_option("-P", "--PROJECT", dest="project", default="",
                 help="project code", metavar="#")
group.add_option("-E", "--EMAIL", dest="email", default="",
                 help="e-mail address", metavar="@")                 
group.add_option("-i", "--info", dest="fn_info", default="",
                 help="mapping: file, barcode, primer, sample information. File should be in CSV format", metavar="allinfo.csv")

parser.add_option_group(group)

group = OptionGroup(parser, "Optional Configuration", description="parameters to alter if necessary")
group.add_option("-a", "--annotations", dest="dir_anno", default="/usr/local/devel/ANNOTATION/sszpakow/ANNOTATION/",
                 help="directory that stores auxilliary files\n[%default]", metavar="annotations")
group.add_option("-S", "--SAMPLE", dest="sampletimes", default=0, type="int",
                 help="perform sub.sampling of all reads based on the number of reads in smallest group. if 0 - all reads are used. if 1 - the sampling will be performed once, if 2 or more, then 2 or more independent samplings are going to be performed.\n[%default]", metavar="#")                 
group.add_option("-m", "--minlen", dest="minlength", default=220, type="int",
                 help="what is the minimum length of reads to process\n[%default]", metavar="#")                 
group.add_option("-x", "--strict", dest="strictlevel", default=1, type="int",
                 help="""how strict to be at demultiplexing: 
                 	1 very strict (barcode no mismatches, primer no mismatches) 
                 	2 less strict (barcode no mismatches, primer 1 mismatch allowed)
                 	[%default]""", metavar="#")                 

group.add_option("-F", "--useFlows", dest="useflows", action="store_true", default=False,
                 help="""if specified flows will be used and pyronoise will be employed. Otherwise qual values will be used with LUCY""", metavar="#") 


parser.add_option_group(group)

group = OptionGroup(parser, "Technical", description="could be useful sometimes")
group.add_option("-C", "--NODESIZE", dest="nodesize", default=30,
                 help="maximum number of grid node's CPUs to use\n[%default]", metavar="#")
parser.add_option_group(group)
	
(options, args) = parser.parse_args()

#################################################
##		Begin
##

	
if options.fn_info == "" or options.email == "" or options.project =="":
	parser.print_help()
	sys.exit(1)
		
init(options.project, options.email)

############################
######################
#### ecoli reference
inputs = {"fasta": ["%s/ecolis.fasta" % options.dir_anno] }
REF = FileImport(inputs)
REF_1 = MakeNamesFile([REF])
REF_2 = MakeGroupsFile([REF], "ref")
REF_3 = MakeQualFile  ([REF], "40" )
#############################

supplementary = list()
READY = deconvolute()
OutputStep("DECON", "fasta,group,name,list,pdf,svg,tiff,taxsummary,globalsummary,localsummary", [READY])

if options.sampletimes==0:
	tmp = finalize(READY)	
	y = R_rarefactions(dict(), dict(), tmp)
	z = R_OTUplots(dict(), dict(), tmp)
	supplementary.append(y)
	supplementary.append(z)
	OutputStep("ENTIRE", "fasta,group,name,list,pdf,svg,tiff,taxsummary,globalsummary,localsummary", [tmp])

else:
	thefinalset = list()
	for k in xrange(0, options.sampletimes):
		args = 	{
					"force" : "fasta,name,group",
					"persample": "T",
					"iter": "%s" % (k)
				}			
		sampled = MothurStep("sub.sample", options.nodesize, dict(), args, [READY])
		tmp = finalize(sampled)	
		y = R_rarefactions(dict(), dict(), tmp)
		z = R_OTUplots(dict(), dict(), tmp)
		supplementary.append(y)
		supplementary.append(z)
		OutputStep("SAMPLED_%s" % (k), "fasta,group,name,list,pdf,svg,tiff,taxsummary,globalsummary,localsummary", [tmp])
		thefinalset.append(tmp)
	
OutputStep("SUPP_PLOTS", "tre,pdf,svg,tiff,r_nseqs,rarefaction,r_simpson,r_invsimpson,r_chao,r_shannon,r_shannoneven,r_coverage", supplementary)
	
	
##########################################################################	
#  
#################################################
##		Finish
#################################################
