########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2013 Sebastian Szpakowski
########################################################################################

#################################################
##  traverser - create a DOT (grahviz) 
##  representation of the performed pipeline steps.
#################################################
import sys, glob
from optparse import OptionParser
from collections import defaultdict
import re

_author="Sebastian Szpakowski"
_date="2011/01/01"
_version="Version 1"

#################################################
##      Classes
##

class   Node:
    def __init__(self, file):
        self.path = file
        self.workpathid = file[1:].strip().split("_")[0]
        self.label =  "_".join(file[1:].strip().split("_")[1:])
        self.has_manifest = False       
        if self.workpathid == "tep":
            self.workpathid = file[1:].strip().split("_")[-1]
            self.label =  "_".join(file[1:].strip().split("_")[1:-1])
        
        #### mapping type - path for files
        self.inputs = defaultdict(set)
        
        #### mapping type - name for files
        self.outputs = defaultdict(set)
        self.outputids = defaultdict(str)
        #### mapping arg  val for program's arguments
        self.arguments= dict()
    
        self.parseManifest()
        
    
    def parseManifest(self):
        try:
            fp = open("%s/%s.manifest" % (self.path, self.workpathid), "r")
            lines=fp.readlines()
            fp.close()
            self.has_manifest=True
        except:
            lines = list()
        counter=0
        for line in lines:
            line = line.strip("\n").split("\t")
            
            if line[0] == "output":
                type = line[1]
                if type not in ("e", "o", "pe", "po", "r"):
                    files = line[2].split(",")
                    
                    if len(files)<5:
                        for file in files:
                            if not file.startswith("[var]"):
                                file = file.split("/")[-1] 
                            
                            else:
                                file = file[5:]
                                counts = len(file.split("-"))
                                if counts>1:
                                    file = "%s\ \(%s\)\ %s" % (type, counts, file)
                                else:
                                    file = "%s\ %s" % (type, file)  
                            
                            
                                
                            if len(file)>120:
                                file = "...%s" % file[-96:]
                            
                            
                                
                            self.outputs[type].add(file)    
                            self.outputids[file] = "<f%s>" % (counter)
                            counter+=1
                    else:
                        id = files[0]
                        id = id.split("/")[-1] 
                        id = id.split(".")[0] 
                        file = "%s [%s files] %s" % ( id, len(files), type )
                        self.outputs[type].add(file)    
                        self.outputids[file] = "<f%s>" % (counter)
                        counter+=1
                    
            elif line[0] == "input":
                type = line[1] 
                files = line[2].split(",")
                
                if len(files)<5:
                    for file in files:
                        if file.startswith("[var]"):
                            file = file[5:]
                            counts = len(file.split("-"))
                            if counts>1:
                                file = "%s\ \(%s\)\ %s" % (type, counts, file)
                            else:
                                file = "%s\ %s" % (type, file)
                            
                        elif len(file.split("~"))>1:
                            file = file.split("~")[0]   
                        else:   
                            file = file.split("/")[-1] 
                            
                            
                        if len(file)>120:
                                file = "...%s" % file[-96:]
                            
                        self.inputs[type].add(file) 
                else:
                
                    ### group by id
                    ids = defaultdict(list)
                    for file in files:
                        id = file
                        id = id.split("/")[-1] 
                        id = id.split(".")[0] 
                        ids[id].append(file.split("/")[-1])
                    
                    
                    for id in ids.keys():
                        files = ids[id]
                        if len(files)<5:
                            for file in files:
                                self.inputs[type].add(file)
                        else:
                            file = "%s [%s files] %s" % ( id, len(files), type )    
                            self.inputs[type].add(file) 
                        #print (file)
                        
            elif line[0] == "argument":
                
                
                if len(line)==2:
                    self.arguments[line[1]] = " "
                    
                elif line[1] in ["postprocess", "awk"]:
                    val = re.escape(line[2].replace("-", "_"))                  
                    #print "huh", val, line[2]
                    self.arguments[line[1]] = val
                    
                else:
                    self.arguments[line[1]] = line[2]   
                    
                
    def getIns(self):
        otpt = list()
        for type, values in self.inputs.items():
            for value in values :
                otpt.append((value, type))
        return (otpt)       
    
    def getOuts(self):
        otpt = list()
        types = self.outputs.keys()
        types.sort()
        
        #for type, values in self.outputs.items():
        for type in types:
            values = self.outputs[type]
            for value in values :
                id = self.outputids[value]
                otpt.append((value, id,  type))
        return (otpt)
    
    def getLabel(self):
        otpt = "%s [%s]" % (self.label, self.workpathid)
        for arg, val in self.arguments.items():
                        
            tmp = list()
            for v in val.split("-"):
                for vv in v.split(","):
                    vv = vv.strip()
                    if len(vv)>0:
                        tmp.append(vv)
                
            otpt = "%s\\n%s -\\> %s " % (otpt, arg, "\\n-\\> ".join(tmp))
        return (otpt)
    
    def getNodeID(self):
        tmp = "%s_%s" % (self.workpathid, self.label)
        #print tmp
        tmp = tmp.replace(".","_").replace("-", "_") 
        return (tmp)
    
    def __str__(self):
        otpt = "[%s]\n%s\n" % (self.workpathid, self.getLabel())
        
        return (otpt)
        
#################################################
##      Functions
##

#################################################
##      Arguments
##
parser = OptionParser()

#parser.add_option("-f", "--file", dest="filename",
#                  help="write report to FILE", metavar="FILE")
#parser.add_option("-q", "--quiet",
#                  action="store_false", dest="verbose", default=True,
#                  help="don't print status messages to stdout")


(options, args) = parser.parse_args()

#################################################
##      Begin
##

nodes = dict()
outs = dict()

for file in glob.glob("./S*_*"):
    file = file.strip("./")
    tmp =  Node(file)
    if tmp.has_manifest:
        nodes [file] = tmp
    
    

# print """digraph Workflow {
#               
# """
#   
# for node in nodes.values():
#   x =  node.getOuts()
#   if len(x)==0:
#       print """node [ width=3 height=3 label="%s" shape=circle style=filled color="green" ] N_%s ;""" % (node.getLabel(), node.getNodeID())
#   else:
#       print """node [ width=3 height=3 label="%s" shape=box style=filled color="gray" ] N_%s ;""" % (node.getLabel(), node.getNodeID())
#   
#   
#   for file, type in node.getOuts():
#       if outs.has_key(file):
#           #print node, file, type
#           pass
#       else:
#           outs[file]=node
#   
# for node in nodes.values():
#   for file, type in node.getIns():
#       if outs.has_key(file):
#           print """edge [ label="%s [%s]" arrowhead=normal penwidth=5 color="black"] N_%s -> N_%s ;""" % (file.replace(".", "\\."), type.replace(".", "\\."), outs[file].getNodeID() , node.getNodeID())          
#       #else:
#       #   print file
#   
#   
# print "}"

print """digraph Workflow {
            graph 
            [
                rankdir = "LR"
                fontsize = "15"
            ];  
"""
    
for node in nodes.values():
    x =  node.getOuts()
    
    label = "<origin> %s" % (node.getLabel())

    
    for file, id, type in node.getOuts():
        if outs.has_key(file):
            #print node, file, type
            pass
        else:
            outs[file]=node
            label = """%s | %s \\"%s\\" """ % (label, id, file )  
    
    if len(x)==0:
        print """
                    "N_%s" 
                    [ 
                            label="%s" 
                            shape="Mrecord"
                            fillcolor="black"
                            penwidth="5"
                            fontcolor="white"
                            fontsize=20
                            color="lightgray"
                            style="filled"
                    ] ;""" % (node.getNodeID(), label )
                    
    else:
        print """
                    "N_%s" 
                    [ 
                            label="%s" 
                            shape="Mrecord"
                            fillcolor="lightgray"
                            penwidth="2"
                            style="filled"
                    ]; """ % (node.getNodeID(), label)
    
    
    
for node in nodes.values():
    for file, type in node.getIns():
        #print file, type
        if outs.has_key(file):
            print """edge [ label="[%s]" arrowhead=normal penwidth=3 color="lightgray"] N_%s : %s -> N_%s : <origin> ;""" % (type.replace(".", "\\."), outs[file].getNodeID(), outs[file].outputids[file], node.getNodeID() )           
        else:
            dummynodeid =file.replace(".", "").replace("/", "").replace("-", "")
            label = file.replace(".", "\.")
            print   """ "N_%s" 
                    [ 
                        label="%s" 
                        shape="folder"
                        fillcolor="red"
                        penwidth="5"
                        style="filled"
                    ] ;""" % ( dummynodeid, label  )
            print """edge [ label="[%s]" arrowhead=normal penwidth=5 color="black"] N_%s -> N_%s : <origin> ;""" % (type.replace(".", "\\."), dummynodeid , node.getNodeID() )          

    
    
print "}"


#################################################
##      Finish
#################################################
