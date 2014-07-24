<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](http://doctoc.herokuapp.com/)*

- [Introduction to YAP](#introduction-to-yap)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

# Introduction to YAP

YAP is a workflow engine that links several independent modules and scripts into a contiguous analytic pipeline.

YAP instances are composed of Steps. YAP orchestrates a specific order of Steps using a multithreaded python interface. Each Step has its own folder named:

```
Step__NAME___hash 
e.g. (Step_classify.otu_1a9385fddedfe918900ef1ba5f84ff8e) 
```

The NAME string comes from the task performed (e.g. input, sed, awk, CD-HIT, Mothur's classify.otu etc.). The 32 character hash is generated based on the sequence of steps that were required to get to the current step. Since the hash is a digest of input files and arguments from the previous steps, YAP knows which step needs to be updated, e.g. if a parameter changed. Upon completion of any given step, a manifest file is generated inside of each Step directory. The manifest contains all the input files and output files.

Upon completion of the workflow, YAP generates a workflow.dot file that lists all the steps that the workflow executed. It highlights the input files and their transition through the workflow to the output steps. The workflow also lists all the parameters that were used in each of the steps.

YAP is currently tailored around the JCVI computational resources. In most cases, it can schedule a grid task in a queue with most available CPU slots, or with most memory, if required.

Yap runs in two modes:

1. The 454 mode requires original SFF files, prior to demultiplexing.
2. The MiSeq mode requires the demultiplexed illumina FASTQ files.

