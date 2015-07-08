########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2013 Sebastian Szpakowski
########################################################################################

#################################################
##    A pipeline for miseq data 
##    OTUs (certain regions of 16S and ITS supported)
##    This is for demultiplexed MiSeq data
#################################################
import YAPGlobals
import sys, os, os.path
from optparse import OptionParser, OptionGroup
from StepsLibrary import *
from StepsLibrary_EXP import *
from collections import defaultdict
from Queue import Queue
##threading redefines enumerate() with no arguments. as a kludge, we drop it here
globals().pop('enumerate',None)

_author="Sebastian Szpakowski"
_date="2013/04/01"
_version="Version 5"

#################################################
##        Classes
##
 
class InfoValidator:
    def __init__(self,filename):
        self.filename = filename
        self.info = GeneralPurposeParser(filename, sep=",")
        self.URI = "http://confluence/display/~sszpakow/YAP"
        self.epilogue = "\n***\tPlease correct before continuing...\n***\t{0}\n".format(self.URI) 
        self.header = ""
        
        self.tech = ""
        
        self.files,  self.barcodes ,self.primersF, self.primersR, self.sampleIDs = self.parse()
        print ("***\tValidation complete, no obvious errors found.\n")
       
               
    def parse(self):
        counter=0;
        print ("\n***\tValidating your template\n\t{0} ...\n".format(self.filename))
        files = set()
        barcodes = set()
        primersF = set()
        primersR = set()
        sampleIDs = set()
                
        for line in pseudo_shuffle(self.info,skip=1):
            if counter == 0: 
                self.header = line
                has =  ",".join (self.header)
                needed454 = "path,file,barcode,forward,reverse,use,both,SampleID"
                neededMiSeq = "path,file1,file2,forward,reverse,SampleID"
                
                if has.lower().startswith( needed454.lower()) :
                    self.tech = "454"
                elif  has.lower().startswith( neededMiSeq.lower()) :
                    self.tech = "MiSeq"
                else:
                    self.error( "Your template's header is incorrect or missing:\nhas :\t{0}\nneed (454):\t{1}\n\t(illumina)\t{2}".format(has, needed454, neededMiSeq), 101)
                
                if not ("SampleID" in self.header):    
                    self.error( "Your template has\n\t'{0}' instead of \n\t'SampleID' in the column's header.".format(self.header[7]), 102)
                    
            else:
                file_added = False
                if line[1].strip() != "":
                    files.add("{0}/{1}".format(line[0], line[1].strip()))   
                    file_added = True
                if self.tech == "454":
                    barcodes.add(line[2])   
                    primersF.add(line[3])    
                    primersR.add(line[4])
                    sampleIDs.add(line[7])
                elif self.tech == "MiSeq":
                    if line[2].strip() != "":
                        files.add("{0}/{1}".format(line[0], line[2].strip())) 
                        file_added = True
                    primersF.add(line[3])    
                    primersR.add(line[4])
                    sampleIDs.add(line[5])
                if not file_added:
                    self.error( "This line in your template file has no files specified: {}".format(",".join(line)), 103)
                        
            counter+=1 
            
        ##### files
        for f in files:
            if not os.path.isfile(f):
                self.error("file doesn't exist\n\t{0}".format(f), 103)

        ##### F primers
        if len(primersF)>1:
            self.error("Multiple forward primers specified:\n\t{0}\n\tnot supported in the current version of YAP".format("\n\t".join(primersF)), 104)
        
        if list(primersF)[0].strip() =="" :
            self.error("Forward primer should not be empty", 104)
        
        
        ##### R primers
        if len(primersF)>1:
            self.error("Multiple reverse primers specified:\n\t{0}\n\tnot supported in the current version of YAP".format("\n\t".join(primersR)), 105)
        
        if list(primersR)[0].strip() =="" :
            self.error("Reverse primer should not be empty", 105)
        
        ##### sampleIDs
        spaces = set()
        ill = ("\\","/", "~", "-", "+", "#")
        illegalchars = set()
        digitstart = set()
        for s in sampleIDs:
            if s.count(" ")>0:
                spaces.add(s)
            for k in ill:
                if s.count(k)>0:
                    illegalchars.add(s)
            if s[0].isdigit():
                digitstart.add(s)
         
        hint = "*You could create two columns: \n\tSampleID, compliant with YAP (excel function: SUBSTITUTE()) and\n\tOriginalIDs, where any character is allowed."    
        if len(spaces)>0:
            M = "The following samplesID(s) have spaces in them:\n\t"
            for s in spaces:
                M = "{0}'{1}',".format(M, s) 
            M = "{0}\n\n\t{1}".format(M, hint)    
            self.error(M, 106)    
            
        if len(illegalchars)>0:
            M = "The following samplesID(s) have illegal chars in them {0}:\n\t".format(", ".join(ill))
            for s in illegalchars:
                M = "{0}'{1}',".format(M, s) 
            
            M = "{0}\n\n\t{1}".format(M, hint)    
            self.error(M, 107)   
            
        if len(digitstart)>0:
            M = "The following samplesID(s) start with numbers:\n\t".format(", ".join(ill))
            for s in digitstart:
                M = "{0}'{1}',".format(M, s) 
                 
            M = "{0}\n\n\t{1}".format(M, hint)    
            self.error(M, 108)  
            
            
        return (files, barcodes, primersF, primersR, sampleIDs)    
                       
                  
    def error(self, message, code):
        print "!!!\t{0}\n{1}".format(message, self.epilogue)
        sys.exit(code)
        
    def getTrimpoints(self):
        return "0", "0", "unknown"

    def getTech(self):
        return self.tech
     
class InfoParserMiSeq:
    def __init__(self, filename):
        self.filename = os.path.abspath(filename)
        self.info = GeneralPurposeParser(filename, sep=",", skip=1)
        self.store = list()
        self.IDs = defaultdict(str)
        self.primers = set()
        self.forward = ""
        self.reverse = ""
        ## scramble the order of samples so that any following partitioning would
        ## not be correlated with their original order
        for line in pseudo_shuffle(self.info,skip=1):
            path = os.path.abspath(line[0].strip())
            file1 = line[1].strip()
            file2 = line[2].strip()
            forward = line[3].strip()
            reverse = line[4].strip()
            
            if path.endswith("/"):
                path = path[:-1]
            
            path1 = "%s/%s" % (path, file1)
            path2 = "%s/%s" % (path, file2)

            ID = line[5].strip()

            paths=[]
                       
            if file1!="" and options.use_mates != "reverse_only":
                paths.append(path1)
                self.IDs[path1] = ID
            if file2!="" and options.use_mates != "forward_only":
                paths.append(path2)
                self.IDs[path2] = ID
            
            if len(paths) == 0:
                print "You have excluded both forward and reverse reads on line: {}".format(",".join(line))
                sys.exit(11)
            if options.mate_merger != "none":
                self.store.append(paths)
            else:
                for path_add in paths:
                    self.store.append([path_add])

            if file1=="" and file2=="":
                print "Both forward and reverse file fields are empty on line: {}".format(",".join(line))
                sys.exit(11) 

            if reverse =="" or forward =="":
                print "%s: please provide both primers for file(s):'%s' " % (x, ",".join(file1, file2))
                sys.exit(11) 
            else:    
                self.primers.add(">_primer_F\n%s\n" % (forward))
                self.primers.add(">_primer_F_rc\n%s\n" % (revComp(forward)))
                self.primers.add(">_primer_R\n%s\n" % (reverse))
                self.primers.add(">_primer_R_rc\n%s\n" % (revComp(reverse)))
                self.forward = forward
                self.reverse = reverse

    def getFiles(self):
        return (self.store)      

    def getSampleID(self, file):
        return self.IDs[file]
       
    def writePrimerFilename(self):
        primerfilename =  "primers.fasta"
        
        if len(self.primers)>4:
            print "The annotation file has more than 2 primers !"
            for p in self.primers:
                print "%s" % (p.strip())
            sys.exit(15)
        
        primerfile = open(primerfilename , "w")     
        
        for p in self.primers:
            primerfile.write(p) 
        primerfile.close() 

        return (primerfilename)
                        
#################################################
##        Functions
##
   
def pcr_reference(
        manifest,
        reference_file,
        out_type,
        step_import_pcr_target,
        pcr_target_idseq,
        pcr_target_padding):
    """Prepare oligos file and call Mothur pcr.seqs on the reference alignment"""
    
    inp = {out_type: [reference_file]}
    step_import_ref_ali = FileImport(inp)

    if not options.no_pcr_reference:
        
        step_pcr = PcrReferenceAlignmentStep(
                dict(pcr_target_idseq=pcr_target_idseq,
                    pcr_target_padding=pcr_target_padding,
                    primer_forward=manifest.forward,
                    primer_reverse=manifest.reverse,
                    align_type=out_type
                    ),
                [step_import_ref_ali,step_import_pcr_target]
                )
        
        step_type = FileTypeMaskInput({"fasta":out_type}, [step_pcr])
        
        step_mask = MaskExceptType(out_type,[step_type])
        
        return [step_mask]

    else:

        return [step_import_ref_ali]
    

def preprocess():
    PREPROCESS = list()

    for (i_sample,files) in enumerate(manifest.getFiles()):
        INS = {}
        
        
        if len(files) == 2:
            M1 = files[0]
            M2 = files[1]
            sampleid = manifest.getSampleID(M1)
            INS = {"mate1": ["%s~%s" % (M1, sampleid)], "mate2": ["%s~%s" % (M2, sampleid)]}
        else:
            M1 = files[0]
            sampleid = manifest.getSampleID(M1)
            INS = {"fastq": ["%s~%s" % (M1, sampleid)]}
          
        #### import files    
        if options.head == 0:
            x = FileImport(INS)
        else:
            x = FileMiniImport(INS, {"lines": options.head})    

        #### determine the encoding of fastQ
        Q = getQ(M1)
        
        if Q == "":
            print (Q)
            print "Q issues"
            print files
            sys.exit(1)
        
        if not options.no_input_qc:
            ### generate quality information:
            ARGS = {
                 "-h": options.minqual,
                 "-m": "",
                 "-v": ""
                }
            qc = SQA(ARGS, [x])
            supplementary.append(qc) 
        
        ## do not split until I implement building table file for input
        ## fastq for make.contigs (it does not understand file lists
        ## for ffastq and rfastq. Maybe this is not ever needed for a typical
        ## case of demultiplexed fastq where each file is already not too large
        if not (len(files)==2 and options.mate_merger == "make.contigs"):
            ### split into smaller files for parallelization
            ### 100,000 sequences (x4 since fastq)
            ARGS = {
                        "types": "mate1,mate2,fastq",
                        "chunk":  "400000"
            }
            P0 = FileSplit(ARGS, [x]) 
        else:
            P0 = x
            
        #### overlap mates if available
        do_final_trim = True
        if len(files)==2:
            if options.minqual_merge > 0: 
                #### trim fastQ files
                ARGS = {
                 "-h": options.minqual_merge,
                }
                P1 = SQAtrim(ARGS, [P0])
            else:
                P1 = P0
            if options.mate_merger == "flash":
                ARGS = {
                 "-M": "200",
                 "-p": Q,
                 "-r": "250"
                 #"-x":"0.15"
                }
                P2_0 = Flash({}, ARGS, [P1])
            elif options.mate_merger == "make.contigs":
                args = { "find_req":"ffastq=mate1,rfastq=mate2",
                         "force": ""}
                P2_0 = MothurStep("make.contigs", options.nodesize, dict(), args, [P1])
                ## this is to ignore scrap.contigs.fasta
                do_final_trim = False
        else:    
            P2_0 = P0
        
        if do_final_trim:
            #### final trim of fastQ files
            ARGS = {
             "-h": options.minqual,
            }
            P2 = SQAtrim(ARGS, [P2_0])
            
            #### convert fastq to fasta
            ARGS = { 
                    "-Q": Q
            }
            P3_1 = fastq2fasta(dict(), ARGS, [P2])
        else:
            P3_1 = P2_0

        P3 = MaskType("fastq",[P3_1])

        #### use fuzznuc to find cut primer sequences
        ARGS = {
                "-f": manifest.forward,
                "-r": manifest.reverse,
                "-m": "1"
        }
        P4 = PrimerClipper ( {}, ARGS, [P3])

        ### make fastA headers less problematic
        P5 = FastaHeadHash({}, { "--prefix":"{}".format(i_sample), "--id-gen":"iid" }, [P4])
        P6 = FileMerger("fasta", [P5])
        P7 = MakeGroupsFile([P6], sampleid)
        P8 = MakeNamesFile([P6])
        
        PREPROCESS.extend([P6,P7,P8])
       
    
    A1 = FileMerger("fasta,group,name", PREPROCESS)
    
    
    args = {"mingroupmembers": options.mingroupmembers, 
            "report": "failing"}
    A2 = GroupRetriever(args, [A1])
    
    args = {
            "force" : "fasta,name,group",
            "find_or_skip_step": "groups" 
            }   
    A3 = MothurStep("remove.groups", options.nodesize, dict(), args, [A2])   

    return [A3]

def finalize(input):
    
    
    clean = CleanFasta(dict(), input)
    
    ####### remove sequences that are too short, and with ambiguous bases 
    args = { "minlength" : "%s" % ( options.minlength ),
             "maxambig" : "0",
             #"maxlength" : "550",
             "maxhomop" : "8",
             "force": "fasta,name,group"}
    clean2 = MothurStep("screen.seqs", options.nodesize, dict(), args, [clean])

    args = {"mingroupmembers": 0, 
            "report": "passing"}
    clean2a = GroupRetriever(args, [clean2])
    OutputStep("2-NOISY", "groupstats,fasta,group,name,list,svg,pdf,tiff,taxsummary,globalsummary,localsummary", clean2a)

    ###################### CDHIT-454
    #### unique and de-noise
    args = {}
    
    ### strictly unique collapsing
    if options.strictlevel==1:
        args=         { 
                        "c" : "1.0",
                        "b" : "8",
                        "aS": "1.0",
                        "g" : "0",
                        "M"    : "50000",
                        "T" : "%s" % (options.nodesize)
                       
                    }
        CD_1 = CDHIT_454(options.nodesize, args, [clean2])
    
    ### aggressive de-noising:
    elif options.strictlevel==2:
        args=         { 
                        "T" : "%s" % (options.nodesize),
                        "splits" : options.preclust_splits
                    }
        CD_1 = CDHIT_Preclust(options.nodesize, args, [clean2])
        
    CD_Moth = CDHIT_Mothurize(dict(), CD_1)
    ##hide cdhit cluster output type so that if final clustering
    ##somehow does not make its own, this would not get picked up
    CD_2aa = MaskType("cdhitclstr",[CD_Moth])

    if options.min_precluster_size > 0:
        ##TODO: replace with Mothur command split.abund
        args = {"min_cluster_size": options.min_precluster_size} 
        CD_2ab =  MakeAccnosFromName(args,[CD_2aa])
        args = {}
        CD_2 = MothurStep("remove.seqs",options.nodesize, dict(), args, [CD_2ab])
    else:
        CD_2 = CD_2aa
       
    args = {"mingroupmembers": 0, 
            "report": "passing"}
    CD_2a = GroupRetriever(args, [CD_2])
    OutputStep("3-UNIQUE", "groupstats,tre,fasta,group,name,list,svg,pdf,tiff,taxsummary,globalsummary,localsummary", CD_2a)  

    
    #### add reference sequences to the merged experiments' file
    CD_3 = FileMerger("fasta,name,group,qfile", [CD_2, REF_1, REF_2, REF_3])
    
    PCR_ref = pcr_reference(
            manifest=manifest,
            reference_file=os.path.join(options.dir_anno,_alignment),
            out_type="refalign",
            step_import_pcr_target=REF,
            pcr_target_idseq=_referenceseqname,
            pcr_target_padding=options.pcr_target_padding)

    #### align to reference database
    args = {    "flip":"t", 
                "ksize": "8",
                "find_req":"reference=refalign"
            }   
     
    CD_4 = MothurStep("align.seqs", options.nodesize, {}, args, PCR_ref+[CD_3])
    
    #### AlignmentSummary determining alignment trimming options 
    #### sets trimstart and trimend variables that can be used by in subsequent steps.
    #### threshold means to keep the center part of the alignment with at least 
    #### the fraction of maximum coverage  
    args = {"ref": _referenceseqname, "thresh": options.dynthresh}
    CD_5 = AlignmentSummary(args,[CD_4])
    
    #### alignment plots  
    if _trimstart != _trimend:
        args = {"ref": _referenceseqname, 
                "trimstart" : _trimstart,  
                "trimend" : _trimend
                }
    else:  
        args = {"ref": _referenceseqname, 
                "trimstart" : "find",  
                "trimend" : "find"
                }        
    CD_6 = AlignmentPlot(args,[CD_5])
    

    #supplementary.append(CD_5)
    supplementary.append(CD_6)
    ###########################
    
    args = {"mingroupmembers": 0, 
            "report": "passing"}
    CD_4a = GroupRetriever(args, [CD_4])
    OutputStep("4-ALIGNED", "groupstats,tre,fasta,group,name,list,svg,pdf,tiff,taxsummary,globalsummary,localsummary", CD_4a)    
       
    cleanCD = cleanup(CD_5)
    args = {"mingroupmembers": 0, 
            "report": "passing"}
    cleanCDa = GroupRetriever(args, cleanCD)
    OutputStep("5-CLEAN", "groupstats,fasta,group,name,list,svg,pdf,tiff,taxsummary,globalsummary,localsummary", cleanCDa)
    
    clusterCD = CDHITCluster(cleanCD)
    
    x = plotsAndStats(clusterCD)
    if not options.no_statistics:
        INS = {"annotation" : [options.fn_info]}
        ARGS = {"dist": "0.03"}
        output1 = R_defaultplots(INS, ARGS, x)
        output2 = AnnotateClusters(dict(), dict(), output1)
    else:
        output2 = x
        
    return (output2)
    

def remove_chimera(input):
    ####### find chimeric sequences  
    toremove = list()
    ## we are already looking at input sequences that were aligned to reference
    ## alignment and cut to the aligned region, so there is no point for us to
    ## consider chimeras that involve other regions, therefore, at first glance
    ## it would make snse to pcr-cut the reference alignment used for chimera
    ## detection. However, what if chimera involves a small chunk of sequence
    ## outside of the targeted region, small enough that it still aligns to the
    ## target region? We keep the full length reference alignment for chimera
    ## detection for now
    for ch in [ "uchime" ]:
        ### chimeras against reference
        args = {"force" : "fasta,reference", "dereplicate" : "t"}
        inputs = {"reference": ["%s/%s" % (options.dir_anno, _alignment_chimera)] }
        
        A = MothurStep("chimera.%s" % (ch),options.nodesize, inputs, args, input)    
        
        toremove.append(A)
        
        if not options.quickmode:
            ### chimeras against self
            ### Seems like a bug in Mothur - it prints errors that non-centroid
            ### (non-unique) sequences from name file are not found in fasta file (and they
            ### should not be). To work around it, we create a count file just for this
            ### chimera.uchime
            A1 = MothurStep("count.seqs",1, {}, {"force":"name,group"}, input)    

            args ={"force_exclude": "name,group", "force": "count,fasta", "dereplicate" : "t"}
            inputs = {}
            
            A2 = MothurStep("chimera.%s" % (ch),options.nodesize, inputs, args, A1)    
            ### hide count file from later steps
            A = MaskType("count",A2)

            toremove.append(A)
        
    ### merge all accnos files and remove ALL chimeras    
    allchimeras = FileMerger("accnos", toremove)


    args ={"force_exclude": "fastq"}
    B = MothurStep("remove.seqs",options.nodesize, dict(), args, allchimeras)
    
    ####### remove empty columns after chimeras were removed
    args = {"vertical" : "T"} #somehow trump=. here removes everything
    out = MothurStep("filter.seqs",options.nodesize, dict(), args, [B]) 
    return [out]

def cleanup(input):

    ### remove the "ref" group
    args = {
            "force" : "fasta,name,group",
            "groups": "ref" 
            }
            
    s15 = MothurStep("remove.groups", options.nodesize, dict(), args, input)
    
    ####### remove sequences that are too short (bad alignment?) 
    ##A.T. Note that screen.seqs() expunges gaps before calculating sequence length,
    ##so if you trimmed the alignment region, then things that would be aligned to
    ##other regions are going to be dropped here because they are too short within the 
    ##alignment window.
    args = {
                "minlength" : "%s" % (options.minlength), 
                "maxambig" : "0",
                "force" : "fasta,name,group" ,
            }
    s16 = MothurStep("screen.seqs", options.nodesize, dict(), args, [s15])
    
    
    #### if primer trimming points are not unknown
    if _trimstart!=_trimend:
        ### primer cut
        args = {
                    "s" : _trimstart, 
                    "e": _trimend,
                }
    else:
        args = {
                "s" : "find:trimstart",  
                "e" : "find:trimend"
        }
     
    ## A.T. This step has no counterpart in the SOP. Here it cuts the alignment
    ## at both ends at 75% of the max observed coverage, and screen.seqs() after 
    ## that drops everything that was
    ## aligned beyond the trimmed window. In the original YAP, the combination
    ## of aligning to the full length database, finding where coverage drops
    ## below 75% of the maximum, cutting and dropping by length must be the way
    ## of getting rid of sequences that align to a wrong region.
    ## In the SOP, they do pcr.seqs() (optional), followed by alligning and dropping 
    ## by mode start and end points, and filter.seqs() that removes all-gap columns and 
    ## and columns that have at least one dot. Thus, there is a difference in effects with YAP: if there is a batch
    ## of 100bp sequences that are aligned in the middle of the alignment, SOP
    ## will drop them too. YAP will keep them. In fact, they increase the max coverage
    ## that YAP computes, and the correspondng 75% coverage cutoff point that it uses
    ## to trim the alignment, thus cutting legitimate alignment regions of longer
    ## sequences. It is also not clear why YAP wants to trim so aggressively the longer alignments that overlap
    ## the target region. Cutting at 75% of the max coverage means, for the start point, that only 25%
    ## of all sequences start to the right of that point.
    if not options.no_trim_alignment:
        s18a = AlignmentTrim(dict(), args, [s16])
    else:
        s18a = s16
            
    ####### remove sequence fragments, bad alignments (?) 
    ##dynamic means alignments were trimmed by now, and bad became short?
    args = {"minlength" : "{}".format(50 if options.dynamic else options.minlength) ,
            "force": "fasta,name,group",
            "force_exclude":"summary"}
    if options.alignment_screen_start_end_criteria < 100:
        args.update({
            "optimize":"start-end",
            "criteria":"{}".format(options.alignment_screen_start_end_criteria)}
            )
    s18b = MothurStep("screen.seqs", options.nodesize, dict(), args, [s18a])

    
    ### build a tree
    #s18b_tree = ClearcutTree({}, s18b)
    
    ####### remove empty columns
    args = {"vertical" : "T"} #sometimes trump=. here removes everything, which means
    ## we still have a dot at every position from at least one sequence
    s19 = MothurStep("filter.seqs",options.nodesize, dict(), args, [s18b]) 
    
    s19a = remove_chimera(s19)
    
    ####### taxonomy
    inputs = {    "reference": ["%s/%s" % (options.dir_anno,_trainset)],
                "taxonomy": ["%s/%s" % (options.dir_anno, _taxonomy )]
            }
            
    args = {    "iters" : "100",
            "cutoff":  "{}".format(options.classify_seqs_cutoff)
            }
    s20 = MothurStep("classify.seqs", options.nodesize, inputs, args, s19a)
    
    ### remove - and . for subsequent clustering efforts 
    s21 = CleanFasta(dict(), [s20])
    
    return [s21]

def CDHITCluster(input):
    cdhits = list()
    for arg in ["0.99", "0.97", "0.95", "0.90"]:
        args = {"c": arg,
                "d" : "0",
                "n": "8",
                "g": "1",
                "M": "0",
                "T": "%s" % (options.nodesize) 
                }
        
        CD_1 = CDHIT_EST(options.nodesize, args, input)
        
        ### make sth. analogous to mothur's labels
        arg = 1.0 - float(arg)
        if arg == 0:
            arg = "unique"
        else:
            arg = "%s" % (arg)
        
        args = {"mode": arg    
                }
        CD_2 = CDHIT_Mothurize(args, CD_1)
        CD_2aa = MothurStep("get.sabund", 1, dict(), {}, [CD_2])
        CD_2ab = MothurStep("get.rabund", 1, dict(), {}, [CD_2])
        CD_2a = CDHIT_Perls({}, [CD_2aa,CD_2ab])            
        cdhits.append(CD_2)
                
    READY = FileMerger("list,rabund,sabund", cdhits)    
    SORTED = FileSort("list,rabund,sabund", READY)
    return [SORTED]

def plotsAndStats(input):
    import re
    
    ### all groups!
    args = {"mingroupmembers": 0, 
            "report": "passing"}
    s23 = GroupRetriever(args, input)
    
    ######## make a shared file 
    labels = ["0.01","0.03","0.05","0.1"]
    args = {"label" : "-".join(labels), "find": "groups"}
    s24_1 = MothurStep("make.shared", options.nodesize, dict(), args, [s23])
    s24 = FileMerger("shared", [s24_1],
            cut_header_lines_others=1,
            order=[re.escape(lab) for lab in labels]) 

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
    
    if options.no_rarefaction:
        return ([s23, s24, s24_1, s25aa, s25bb, s26a, s27a, s28])
    else:
        args = {"force" : "list", "calc": "nseqs-sobs-simpson-invsimpson-chao-shannon-shannoneven-coverage", "freq": "0.01"}
        s29 = MothurStep("rarefaction.single", options.nodesize, dict(), args, [s24])
        args = {"force" : "shared", "calc": "nseqs-sobs-simpson-invsimpson-chao-shannon-shannoneven-coverage", "freq": "0.05"}
        s30 = MothurStep("rarefaction.single",options.nodesize, dict(), args, [s24]) 
        return ([s23, s24, s24_1, s25aa, s25bb, s26a, s27a, s28, s29, s30])
    
    
#################################################
##        Arguments
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

group.add_option("-Y", "--Yap", dest="mode", default="16S",
                 help="""Which Pipeline: 16S ITS [%default]""", metavar="#") 

group.add_option("-D", "--dynamic", dest="dynamic", action = "store_true", default=True,
                 help="""If specified, alignment will be scanned for primer locations and trimmed accordingly. Otherwise a database of known primers and trimming points will be used. [%default]""", metavar="#") 

group.add_option("-d", "--trim-alignment-thresh", dest="dynthresh", default=0.75, type="float",
                 help="""in conjunction with -D, otherwise this is ignored. This allows to specify how much of the alignment to keep using the per-base coverage. The [%default] value indicates that ends of the alignment are trimmed until a base has a coverage of [%default] * peak coverage.""", metavar="#") 

group.add_option("-a", "--annotations", dest="dir_anno", default=os.environ["YAP_DATA"]+"/",
                 help="directory that stores auxilliary files\n[%default]", metavar="annotations")
group.add_option("-S", "--SAMPLE", dest="sampletimes", default=0, type="int",
                 help="perform sub.sampling of all reads based on the number of reads in smallest group. if 0 - all reads are used. if 1 - the sampling will be performed once, if 2 or more, then 2 or more independent samplings are going to be performed.\n[%default]", metavar="#")                 
group.add_option("-m", "--minlen", dest="minlength", default=200, type="int",
                 help="what is the minimum length of reads to process\n[%default]", metavar="#")     

group.add_option("-g", "--mingroupsize", dest="mingroupmembers", default=100, type="int",
                 help="after demultiplexing, discard groups with fewer reads than #\n[%default]", metavar="#")

group.add_option("--min-precluster-size", dest="min_precluster_size", default=2, type="int",
                 help="after pre-clustering, discard clusters with fewer sequences than #\n[%default]. Set to 2 to discard singletons.", metavar="#")

group.add_option("-Z", "--minqual-before-pair-merge", dest="minqual_merge", default=3, type="int",
                 help="Keep stretches of reads this good or better before merging paired reads (zero means no trimming)#\n[%default]", metavar="#")

group.add_option("-M", "--mate-merger", dest="mate_merger", default="make.contigs", type="choice",
                 choices=("make.contigs","flash","none"),
                 help="Method for merging paired-end reads into contigs\n[%default]", metavar="mate_merger")

group.add_option("--use-mates", dest="use_mates", default="both", type="choice",
                 choices=("both","forward_only","reverse_only"),
                 help="Choice to use 'both','forward_only' or 'reverse_only' paired end reads. Combined with --mate-merger option\n[%default]", metavar="use_mates")

group.add_option("-Q", "--minqual", dest="minqual", default=30, type="int",
                 help="Keep stretches of reads this good or better (if merging paired reads, this is done after merging and only if merging method produces FASTQ - see also --minqual-before-pair-merge) #\n[%default]", metavar="#")

group.add_option("--classify-seqs-cutoff", dest="classify_seqs_cutoff", default=60, type="int",
                 help="Mothur classify.seqs(cutoff) parameter#\n[%default]", metavar="#")

group.add_option("--pcr-target-padding", dest="pcr_target_padding", default=10, type="int",
                 help="When trimming reference alignment, pad target reference sequence (ecoli) by that many bases beyond primer locations#\n[%default]", metavar="#")

group.add_option("-q", "--quick", dest="quickmode", action = "store_true", default=False,
                 help="""If specified, only single, reference DB based chimera checking will be used. [%default]""", metavar="#") 

group.add_option("-s", "--no-statistics", dest="no_statistics", action = "store_true", default=False,
                 help="""If set, do not do statistical analysis (future default). [%default]""", metavar="#") 
 
group.add_option("--no-rarefaction", dest="no_rarefaction", action = "store_true", default=False,
                 help="""If set, do not do run rarefaction analysis (that can take a very long time for large sample sizes). [%default]""", metavar="#") 
 
group.add_option("--no-input-qc", dest="no_input_qc", action = "store_true", default=False,
                 help="""If set, do not run QC report on the input sequence files. [%default]""", metavar="#") 
 
group.add_option("--pcr-reference", dest="no_pcr_reference", action = "store_false", default=True,
                 help="""Cut reference alignment to the region defined by the primers in the manifest. [False]""", metavar="#") 
 
group.add_option("--no-trim-alignment", dest="no_trim_alignment", action = "store_true", default=False,
                 help="""Do not trim alignment of sequences to the reference alignment based on coverage. [%default]""", metavar="#") 

group.add_option("--alignment-screen-start-end-criteria", dest="alignment_screen_start_end_criteria", default=100, type="int",
                 help="screen.seqs(fasta=alignment,optimize=start-end,criteria=?). Set to 100 to keep all sequences.#\n[%default].", metavar="#")
 
parser.add_option("-H", "--head", dest="head", default=0, type="int",
                 help="For dry runs, import only # of lines from the input files")

              
group.add_option("-x", "--strict", dest="strictlevel", default=2, type="int",
                 help="""how strict to be at pre-clustering: 
1 very strict, conservative denoising (precluster identical sequences) 
2 less strict, aggresive denoising (precluster using 98% similarity)
[%default]""", metavar="#")                 

parser.add_option_group(group)

group = OptionGroup(parser, "Technical", description="could be useful sometimes")
group.add_option("-C", "--NODESIZE", dest="nodesize", default=16,
                 help="maximum number of grid node's CPUs to use\n[%default]", metavar="#")
group.add_option("-G", "--debug-grid-tasks", dest="debug_grid_tasks", action = "store_true", default=False,
                 help="Debug GridTasks by default\n[%default]", metavar="#")
group.add_option("-T", "--step-dummy-thread", dest="step_dummy_thread", action = "store_true", default=False,
                 help="Use dummy threads inside the main thread for StepXXX classes (for interactive debugging)\n[%default]", metavar="#")
group.add_option("--dummy-grid-tasks", dest="dummy_grid_tasks", action = "store_true", default=False,
                 help="Use dummy grid tasks that run tasks inside the current process in a blocking subprocess (for debugging). Probably use with dummy threads or you can flood the current node from multiple threads\n[%default]", metavar="#")
group.add_option("--large-run", dest="large_run", action = "store_true", default=False,
                 help="This will be a large scale run, modify behaviour in some places for scalability\n[%default]", metavar="#")
group.add_option("--preclust-splits", dest="preclust_splits", default=10,
                 help="Number of data splits in pre-clustering step. Might be useful if in a large scale run the CDHIT 454 preclust is killed\n[%default]", metavar="#")
parser.add_option_group(group)

(options, args) = parser.parse_args()

YAPGlobals.debug_grid_tasks = options.debug_grid_tasks
YAPGlobals.step_dummy_thread = options.step_dummy_thread
YAPGlobals.dummy_grid_tasks = options.dummy_grid_tasks
YAPGlobals.large_run = options.large_run

#################################################
##        Begin
##

    
if options.fn_info == "" or options.email == "" or options.project =="":
    parser.print_help()
    sys.exit(1)
     
if not options.mode in ("16S", "ITS"):
    parser.print_help()
    sys.exit(2)    

### parameters specific to YAP incarnations

### 16S V1-V3    
if options.mode=="16S":
    ### file in the annotations directory that has reference sequences
    _referenceseq = "ecolis.fasta"
    ### which fasta ID use as the reference (if file has more than one)
    _referenceseqname = "e_coli2_genbank"
    ### mothur's compendium of ALIGNED 16S sequences
    _alignment = "silva.seed_v119.align"
    #DEBUG:
    #_alignment = "silva.bacteria.fasta"
    ### mothur's compendium of ALIGNED 16S sequences for chimera detection
    _alignment_chimera = "silva.gold.align"
    ### mothur's curated version of RDP's curated train set and corresponding taxonomy
    _trainset = "trainset10_082014.pds.fasta"
    _taxonomy = "trainset10_082014.pds.tax"
    
### ITS NSI1 - NLB4 (barcoded)   
elif options.mode=="ITS":
    _referenceseq = "yeastITS.fasta"
    _referenceseqname = "AF293_reference"
    _alignment = "FungalITSseed.092012.1.aln.fasta"
    _alignment_chimera = _alignment
    _trainset = "FungalITSdb.092012.1.fasta"
    _taxonomy = "FungalITSdb.092012.1.tax"

else:
    parser.print_help()
    sys.exit(2)
                   
validator = InfoValidator(options.fn_info)  
_trimstart , _trimend, _region = validator.getTrimpoints()    
_tech = validator.getTech()

O = list()
init_res = init(options.project, options.email)
BOH = init_res["BOH"]
QS = init_res["QS"]
MOTHUR = init_res["MOTHUR"]

try:
    try:
        BOH.toPrint("-----", "GLOBAL",  "We are in %s mode" % (options.mode)) 
        BOH.toPrint("-----", "GLOBAL",  "We will be processing %s data" % (_tech)) 

        if options.dynamic or _region == "unknown":
            BOH.toPrint("-----", "GLOBAL",  "Dynamic alignment trimming enabled")
            BOH.toPrint("-----", "GLOBAL",  "Alignment will be trimmed using %s * peak coverage threshold" % (options.dynthresh))
            _trimstart = "0"
            _trimend = "0"
        else:
            BOH.toPrint("-----", "GLOBAL",  "Alignment trimming predefined: %s - %s" % (_trimstart, _trimend))

        manifest = InfoParserMiSeq(options.fn_info)

#############################
#######################
##### reference: 
        inputs = {"fasta": ["%s/%s" % (options.dir_anno, _referenceseq)] }
        REF = FileImport(inputs)
        REF_1 = MakeNamesFile([REF])
        REF_2 = MakeGroupsFile([REF], "ref")
        REF_3 = MakeQualFile  ([REF], "40" )
##############################

        supplementary = list()
        READY = preprocess()
        O.append(OutputStep("1-PREPROCESS", "groupstats,fasta,group,name,list,pdf,svg,tiff,taxsummary,globalsummary,localsummary", READY))

        if options.sampletimes==0:
            fin = finalize(READY)   
            if not options.no_rarefaction:
                y = R_rarefactions(dict(), dict(), fin)
                z = R_OTUplots(dict(), dict(), fin)
                supplementary.append(y)
                supplementary.append(z)
            O.append(OutputStep("6-ENTIRE", "taxonomy,shared,groupstats,fasta,group,name,list,pdf,svg,tiff,taxsummary,globalsummary,localsummary,phylotax", fin))
            O.append(OutputStep("8-TBC", "phylotax,group,list,fasta", fin))
            
#else:
#    thefinalset = list()
#    for k in xrange(0, options.sampletimes):
#        args =     {
#                    "force" : "fasta,name,group",
#                    "persample": "T",
#                    "iter": "%s" % (k)
#                }            
#        sampled = MothurStep("sub.sample", options.nodesize, dict(), args, [READY])
#        tmp = finalize(sampled)    
#        y = R_rarefactions(dict(), dict(), tmp)
#        z = R_OTUplots(dict(), dict(), tmp)
#        supplementary.append(y)
#        supplementary.append(z)
#        OutputStep("SAMPLED_%s" % (k), "groupstats,fasta,group,name,list,pdf,svg,tiff,taxsummary,globalsummary,localsummary", [tmp])
#        thefinalset.append(tmp)
#    
        O.append(OutputStep("7-SUPP_PLOTS", "tre,pdf,png,svg,tiff,r_nseqs,rarefaction,r_simpson,r_invsimpson,r_chao,r_shannon,r_shannoneven,r_coverage", supplementary))
    except:
        BOH.stop()
        QS.stop()
        raise

finally:
    for o in O:
        o.join()

    #if YAPGlobals.step_dummy_thread:
    BOH.stop()
    QS.stop()

    BOH.join()
    QS.join()
    
###########################################################################    
##  
##################################################
###        Finish
##################################################
