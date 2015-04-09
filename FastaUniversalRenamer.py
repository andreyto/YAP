########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2013 Sebastian Szpakowski
########################################################################################


#################################################
##     Rename fastaIDs with "x" +  "hash (id)"
##     new
#################################################
import sys, argparse, hashlib

_author="Sebastian Szpakowski"
_date="Feb 7, 2013"
_version="V1.00"

#################################################
##        Classes
##
    #################################################
    ### Iterator over input fasta file.
    ### Only reading when requested
    ### Useful for very large FASTA files
    ### with many sequences
    
class FastaParser:
    def __init__ (self, x, quals=False):
        self.filename = x
        self.fp = open(x, "r")    
        self.currline = "" 
        self.currentFastaName = ""
        self.currentFastaSequence = ""
        self.lastitem=False    
        if quals:
            self.linesep=" "    
        else:
            self.linesep=""    
    def __iter__(self):
        return(self)    
                
        ##### 
    def next(self):
        for self.currline in self.fp:
            if self.currline.startswith(">"):
                self.currline = self.currline[1:]
                if self.currentFastaName == "":
                    self.currentFastaName = self.currline
                else:
                    otpt = (self.currentFastaName.strip(), self.currentFastaSequence.strip())
                    self.currentFastaName = self.currline
                    self.currentFastaSequence = ""    
                    self.previoustell = self.fp.tell()
                    return (otpt)
                
            else:
                self.addSequence(self.currline)    
        
        if not self.lastitem:
            self.lastitem=True            
            return (self.currentFastaName.strip(), self.currentFastaSequence.strip())
        else:
            raise StopIteration    
                                   
    def addSequence(self, x):
               self.currentFastaSequence = "%s%s%s" % (self.currentFastaSequence,self.linesep, x.strip())            
    def __str__():
        return ("reading file: %s" %self.filename) 
        

class iid_gen:
    
    def __init__(self,start=0):
        self.curr = start

    def __call__(self,head):
        self.curr += 1
        return str(self.curr)

#################################################
##        Functions
##

#################################################
##        Arguments
##

_arguments = argparse.ArgumentParser(description="""replace the entire line of fasta header with a UUID \
        or a combination of user provided prefix and incremented integer index""")
_arguments.add_argument('--fasta', metavar='N', type=str, dest = "fn_fasta")
_arguments.add_argument('--prefix', metavar='N', type=str, dest = "prefix", default="")
_arguments.add_argument('--id-gen', metavar='N', type=str, dest = "id_gen", default="uuid32", choices = ("iid","uuid32","sha256","md5"))
args = _arguments.parse_args()


#################################################
##        Begin
##

newfileFASTA   = "{0}.idhash.fasta".format( ".".join(args.fn_fasta.strip().split(".")[:-1]))
newfileMAPPING = "{0}.idhash.mapping".format( ".".join(args.fn_fasta.strip().split(".")[:-1]))

otptseq = open(newfileFASTA, "w")
otptmap = open(newfileMAPPING, "w")

if args.id_gen == "sha256":
    id_gen = lambda h: hashlib.sha256(head).hexdigest()
elif args.id_gen == "md5":
    id_gen = lambda h: hashlib.md5(head).hexdigest()
elif args.id_gen == "uuid32":
    import uuid
    id_gen = lambda h: uuid.uuid4().hex
elif args.id_gen == "iid":
    id_gen = iid_gen()
else:
    raise ValueError(args.id_gen)

prefix = args.prefix
if prefix:
    prefix = prefix + "_"

for head, seq in FastaParser(args.fn_fasta):
    ## sha256 hex is 64 byte; md5 is 32 byte
    newhead =  id_gen(head)
    otptseq.write(">{2}{0}\n{1}\n".format(newhead,seq, prefix) )
    otptmap.write("{2}{0}\t{1}\n".format(newhead, head, prefix) )          
    
otptseq.close()
otptmap.close()

#################################################
##        Finish
#################################################

