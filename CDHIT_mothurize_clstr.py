########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2013 Sebastian Szpakowski
########################################################################################

#################################################
##  create or merge with existing names file an output from
##  CD-HIT
#################################################
import sys
from optparse import OptionParser
import itertools

_author="Andrey Tovchigrechko"
_date="2015/04/08"
_version="Version 2"

#################################################
##      Methods
##

def clstr_reader(file_name):
    """Iterate over sequence IDs in CDHIT cluster file.
    Yield tuples (representative ID, list of other IDs 
    including the representative."""

    with open(file_name,"r") as inp:
        line = inp.next()
        assert line[0] == ">"
        repr = None
        clust = []
        for line in inp:
            if line[0] == ">":
                assert repr and clust
                yield (repr,clust)
                repr = None
                clust = []
            else:
                id_seq = line.split(">",1)[1].split("...",1)[0]
                clust.append(id_seq)
                if line.rstrip()[-1] == "*":
                    repr = id_seq
        if clust:
            yield (repr,clust)

def iter_str_split(s,sep):
    """Constant memory (iterator) string splitter.
    Note that this will not work correctly on empty fields (e.g. ',,')
    """

    return (x.group(0) for x in re.finditer(r"[^{}]+".format(sep), s))

def iter_elems_with_names(elems,names=None,no_split=True, repr=None):
    """Join elems iterable with optional names dict, where elems are
    keys in the dict, and iterate over all values.
    @param no_split If True, do not split names values by commas"""
    if names:
        if no_split:
            if repr:
                return itertools.chain.from_iterable(
                        (iter_cut_word(names[elem],repr,",") if elem == repr else (names[elem],) for elem in elems)
                        )
            else:
                return (( names[elem] for elem in elems ))
        else:
            if repr:
                return itertools.chain.from_iterable(( iter_str_split(names[elem],",") for elem in elems if elem != repr ))
            else:
                return itertools.chain.from_iterable(( iter_str_split(names[elem],",") for elem in elems ))
    else:
        return elems

def iter_cut_word(s,w,sep,required=True):
    """'aaa,bbb,xxx,zzz' -> ('aaa','xxx,zzz')"""
    if s == w:
        pass
    else:
        ind = s.find(sep+w+sep)
        if ind > 0:
            yield s[0:ind]
            yield s[ind+len(sep+w+sep):]
        elif ind < 0:
            if s.endswith(sep+w):
                yield s[:-len(sep+w)]
            elif s.startswith(w+sep):
                yield s[len(w+sep):]
            elif required:
                raise ValueError("Match for {} not found".format(w))
        else:
            raise ValueError("String should not start with a separator")


#################################################
##      Arguments
##
parser = OptionParser()

parser.add_option("-c", "--cluster", dest="fn_clstr",
                 help="clustering from CDHIT", metavar="FILE")

parser.add_option("-n", "--names", dest="fn_name", default ="",
                 help="names from CDHIT", metavar="FILE")
                 
parser.add_option("-o", "--output_mode", dest="out_mode", default = "name",
                 help="output mode, i.e. either output list or names files", metavar="FILE")
                 
#parser.add_option("-q", "--quiet",
#                  action="store_false", dest="verbose", default=True,
#                  help="don't print status messages to stdout")



(options, args) = parser.parse_args()

#################################################
##      Begin
##


if options.fn_name != "" :
    with open(options.fn_name,"r") as inp_names:
        names = dict(( line.strip().split("\t",1) for line in inp_names ))

if options.out_mode =="name":
    
    prefix = options.fn_clstr.strip().split("/")[-1]
    
    with open("%s.name" % (prefix), "w") as out:
        
        for (repr,elems) in clstr_reader(options.fn_clstr):
            out.write("{}\t{}\n".format(repr,",".join(
                iter_elems_with_names(elems,names)
                )))
    
else:
    ## we only write list file - sabund and rabund should be created from list (and optional count)
    ## files by Mothur commands get.sabund and get.rabund

    ## we make two passes to avoid loading entire cluster file into RAM
    n_otus = 0
    for rec in clstr_reader(options.fn_clstr):
        n_otus += 1
    prefix = options.fn_clstr.strip().split("/")[-1]
    with open("%s.list" % (prefix), "w" ) as out:
        out.write("{}\t{}".format(options.out_mode, n_otus))
        for (repr,elems) in clstr_reader(options.fn_clstr):
            ## repr should be listed first, so we write it out and exclude from the elems
            out.write("\t{}".format(repr))
            ## if something else is in either new or old cluster aside from repr
            if len(elems) > 1 or names[repr] != repr:
                out.write(",")
                out.write(",".join(
                    iter_elems_with_names(elems,names,repr=repr)
                    ))
        out.write("\n")
    
#################################################
##      Finish
#################################################

