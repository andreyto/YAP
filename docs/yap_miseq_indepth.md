<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](http://doctoc.herokuapp.com/)*

- [YAP: The MiSeq mode](#yap-the-miseq-mode)
- [Introduction](#introduction)
      - [Using and citing YAP](#using-and-citing-yap)
  - [Overview](#overview)
  - [Location](#location)
  - [Example Usage](#example-usage)
  - [YAP Flags](#yap-flags)
  - [Input File Format](#input-file-format)
    - [The required columns and brief description of what they do.](#the-required-columns-and-brief-description-of-what-they-do)
  - [What to expect in output](#what-to-expect-in-output)
    - [What to keep and why](#what-to-keep-and-why)
    - [Re-running the pipeline](#re-running-the-pipeline)
- [YAP Walkthrough](#yap-walkthrough)
    - [Before you begin...](#before-you-begin)
    - [Running YAP](#running-yap)
  - [YAP Flags](#yap-flags-1)
- [YAP: Methods and References](#yap-methods-and-references)
    - [Blurb](#blurb)
    - [References](#references)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

# YAP: The MiSeq mode


In the YAP MiSeq mode, YAP processes demultiplexed Illumina FASTQ files amplified from bacterial 16S and fungal ITS regions.


# Introduction
There is no single analytic workflow that will produce a paper and publication-ready graphs. Please bear in mind that that the YAP workflow is intended as a starting point in the exploration of your data.

#### Using and citing YAP ####
If you use YAP, don't forget to cite this source as well the authors of all the independent modules within the workflow. Please see YAP: Methods for a complete list.

## Overview

The YAP_MiSeq.sh script wraps around the following utilities and resources:

*   [Mothur v.1.29.2](http://www.mothur.org/ "Mothur Website")
*   [CD-HIT v.4.6.1](http://weizhong-lab.ucsd.edu/cd-hit/ "CD-HIT")
*   [LUCY](http://lucy.sourceforge.net/ "LUCY software")
*   Custom [R](http://cran.r-project.org/ "R Statistical Language Website") and [python](http://python.org/ "python programming language website") scripts by Drs. Szpakowski and Tovchigrechko. to facilitate a smooth integration and visualization of the output.

*   The following are installed as a part of the YAP package
*   
        *   The R Statistical Language  
        *   ade4 (for PCA analyses)
        *   RColorBrewer (for color palettes)
        *   fdrtool 
        *   Biostrings for fastA file manipulation
        *   ggplot2 for graphing
        *   **Python 2.7+** is required.
    
## Location

`/usr/local/devel/ANNOTATION/sszpakow/YAP/{scripts,bin,example}`

In general, to run YAP the user needs to edit the **~/.profile** and **~/.bashrc** files so that they contain:

```
export LD_LIBRARY_PATH=/usr/local/packages/gcc/lib64:$LD_LIBRARY_PATH</pre>
```

## Example Usage

```
/usr/local/devel/ANNOTATION/sszpakow/YAP/bin/YAP_MiSeq.sh -Y 16S -P ##### -E youremail@email.com -i /path/to/allinfo.csv
```

where:

* `YAP_MiSeq.sh` is a convenience wrapper and the orchestrator of steps,
* `-Y 16S` selects the mode - 16S. 16S is the default so this really does not have to be specified, unless some other mode is requested (e.g. ITS)
* `-P #####` is the project ID, required for accounting on the grid
* `-E [XXXXX@jcvi.org](mailto:youremail@email.com)` is your e-mail address, to notify you if a problem occurred,
* `-i (...)allinfo.csv` is the path to the input configuration and annotation file

*Please note*: The entire path (or the absolute path) needs to be spelled out for YAP: beginning with the root or "/", traversing through all the directories, e.g. usr/local/projects/aProject/ and ending with the file name and the ".csv" suffix. The absolute path is required so that the nodes of the GRID will know exactly where to look for the file.


## YAP Flags

Running this command:

```
/usr/local/devel/ANNOTATION/sszpakow/YAP/bin/YAP_454.py
```
will simply print out available parameters that one can modify if necessary.


## Input File Format

The input file is a [comma-separated-value (CSV) file](http://en.wikipedia.org/wiki/Comma-separated_values "information about CSV format"). It's openable and editable by spreadsheet programs, but it has to be saved as CSV. Binary spreadsheet format such as ".xlsx" will NOT work in the current implementation of the program.

Several columns are required. The header names are fixed in both wording and casing based on the expectations of programs down the line. Any deviation from the defaults in the required column names and letter cases will cause errors.

### The required columns and brief description of what they do.

path |  file1  | file2| forward | reverse |	**SampleID** |	GroupingA |	GroupingB 
-----|---------|------|---------|-------- | -------------|----------- |-----------
REQUIRED |	REQUIRED |	REQUIRED |	REQUIRED |	REQUIRED |	REQUIRED |	optional | optional 


Where:

* path: Absolute path to where the fastq file is located, the path should end with "/"
* file1: the name of the mate1 fastq file
* the name of the mate2 fastq file. Leave blank if only interested in single mate analysis.
* forward primer sequence relevant to the the particular sample 
* reverse primer sequence relevant to the particular sample 
* sample name. This column requires strict formatting!! 
    + It is crucial that:
    + no sample name start with a digit. 
    + no sample name contain spaces. 
    + If two samples have the same sampleID, that means that the information associated with these two barcodes will be averaged and treated as a single sample.


An example input *.csv file should look like this:

path |  file1  | file2| forward | reverse |SampleID     |Experimentname|Samplegroup
---  |---------|------|---------|-------- | ------------|--------------|-----------
... |A_R1.fastq|A_R2.fastq	|AGAGTTTGATYMTGGCTCAG|ATTACCGCGGCTGCTGG|Mouse766|A	| Control		 
... |B_R1.fastq|B_R2.fastq	|AGAGTTTGATYMTGGCTCAG|ATTACCGCGGCTGCTGG|Mouse767|B	| Control	 
... |C_R1.fastq|C_R2.fastq	|AGAGTTTGATYMTGGCTCAG|ATTACCGCGGCTGCTGG|Mouse768|A	| Treatment		
... |D_R1.fastq|D_R2.fastq	|AGAGTTTGATYMTGGCTCAG|ATTACCGCGGCTGCTGG|Mouse769|B	| Treatment		

_Note_: ... indicates the absolute path.

This example file represents a sequencing effort of four mouse samples (four different entries in the SampleID column). Specifying forward primer and reverse primer sequences will essentially search and trim the data twice, once requiring the forward primer, second requiring the reverse primer. In the end the clean results will be merged.

In this example, the four samples are divided among two different experiments A and B (two different unique entries in the Experiment column). Each experiment has one control sample treated with a placebo and one treatment sample.

The YAP workflow will cleanup all the data, average the data and create some preliminary plots labeled with the SampleID, the optional Experiment and Samplegroup columns. Please note that the optional columns are completely up to the user both in name and in content. Columns names *should not* start with numbers.

## What to expect in output

* To see the order of steps an overview of the workflow performed, please review workflow.dot.
* The file called `logfile.txt` contains a log of the entire procedure.

Briefly, there are a number of directories that store data for the intermediate steps of the analysis. The convention for the directory naming is as follows:

```
Step_{Step_Name}_{step_sequence_id}
```

The most relevant contain the word: "OUTPUT" in them.

As the workflow progresses, the following OUTPUT directories are created.

*   **STEP_OUTPUT_1-PREPROCESSING** - stores all relevant files after demultiplexing step and quality trimming (fasta, group, name, groupstats - counts)
*   **STEP_OUTPUT_2-NOISY** - captures the first fragment check so that sequences with ambiguous (N) bases are filtered out, only sequences 220-2500 bases in length remain. (fasta, group, name, groupstats)
*   **STEP_OUTPUT_3-UNIQUE** - files post pre-clustering step featuring a unique set sequences (cdhit/fasta, name, group and groupstats) remain
*   **STEP_OUTPUT_4-ALIGNED** - files storing sequences after unique sequences are aligned against a database of rDNA sequences (fasta/align, group, names)
*   **STEP_OUTPUT_5-CLEAN **- after chimera Slayer removes potentially chimeric sequences, bad alignments are removed, sequences outside of the primer ranges are clipped, and erroneous, misaligned sequences are filtered out. The sequences are then processed so that any &quot;.&quot; or &quot;-&quot; characters are removed and RDP classification occurs. The files in this directory are the cleanest yet (fasta, group, name, taxonomy)
*   **STEP_OUTPUT_6-ENTIRE** - contains pdf plots of rarefaction curves, and the results of sequence/similarity clustering.
*   **STEP_OUTPUT_7-SUPP_PLOTS** - stores pdf files with read alignment statistics, OTU size histograms, as well as rarefaction results local (per sample group) and global (cumulative), in text and visual forms. Furthermore files
*   **STEP_OUTPUT_8-TBC **- files for post processing (to be continued) &nbsp;

Based on the example file earlier:
 **OUTPUT_6-ENTIRE** will contain a whole set of pdf files containing graphics that were generated based on the columns and markup of the experiments within the input file.

### What to keep and why

1.  The days of the files on scratch are numbered, up to seven days since their creation they WILL BE DELETED without a warning and without a backup.
2.  `OUTPUT_ENTIRE` and `OUTPUT_SUPP_PLOTS` are probably the two most important directories that you will want to copy from scratch and store in a safe, archivable location
3.  `logfile.txt` and `workflow.dot` store the information about the pipeline's run (parameters and commands), so it's a good idea to keep them with the OUTPUTs as well.
    * logfile stores all steps and commands run in the process of the analysis. The very command to initiate the pipeline is also stored in the logfile.
    * workflow.dot has information of how steps relate to each other and which command uses what parameters.

### Re-running the pipeline

1.  Unless parameters change, rerunning the pipeline will not overwrite any files. The software checks which steps have completed and will not repeat them.
2.  If you ever modify the csv input file (e.g. to add more columns to stratify the samples for more informative plots), please delete the `R_plots` folder and `OUTPUT_ENTIRE *OUTPUT_SUPP_PLOTS` and then rerun the pipeline with the exact parameters as before (check your logfile). If everything works properly, the missing steps will be detected and recreated.

# YAP Walkthrough

### Before you begin...

1.  YAP spawns other programs to the grid. In order to achieve that, however, you need to run YAP from a computer (or host) that is "submit-enabled". `lserver` qlogin machines should work fine.
2.  If a command `qstat` or `qsub` returns a "command not found" error your computer needs a small configuration change.

        1.  Please see your grid page for adding your host as a submit host.
        2.  If you dont want to run the `source` command outlined above, you can edit ~/.bashrc script and add the following JCVI grid specific line:
        ```
        source /usr/local/sge_current/jcvi/common/settings.sh
        ```
        This way every time you log in to your submit host, configuration will be applied automatically.
        3.  To run YAP the user needs to edit the **~/.profile** and **~/.bashrc** files so that they contain:

```
export LD_LIBRARY_PATH=/usr/local/packages/gcc/lib64:$LD_LIBRARY_PATH</pre>
```

The following are only necessary if you are copying YAP somewhere else...

1.  This program requires a few R packages. The following procedure can be used to semi-automatically install them:

    1.  Start R (type **`R`** in the terminal)
    2.  **`source("/usr/local/devel/ANNOTATION/sszpakow/YAP/scripts/prereqs.R")`**
    3.  if asked to install libraries locally, type **`y`** for yes
    4.  when all is done , type **`q()`** to quit R
  
  
2. This program requires Python 2.7.3, Numpy and Biopython.

### Running YAP

1.  Create a new directory in **`/usr/local/scratch/MYUSERNAME/MYDIRECTORY`** that will house your analysis.
2.  Within the newly created directory make two more: **`data`** and **`analysis`**
3.  Copy FASTQ files from original location to the data directory
4.  If you are making symbolic links, please make sure that the linked location is reachable from the grid. (level 3 storage is not, for example)
4.  Create a spreadsheet that matches the specifications described in the previous section.
    1.  Update the path fields to contain **`/usr/local/scratch/MYUSERNAME/MYDIRECTORY/data/`** (or wherever you saved your FASTQ files)
    2.  Update the filename fields so that they contains the filenames of SFF file at the specified path
    3.  Update primer sequences if necessary (most likely just copy and paste )
    4.  Update the barcodes to match contents of the sff file. barcode ACTG in sff file 1 may contain different samples than the same barcode in sff file 2.
    5.  Create additional columns that will be used to generate plots of interest. There is no defined limit how many columns to include or how many plots will be generated. See more details below about what you can do.
6.  Save the spreadsheet in the data directory as a CSV file. In excel use the `>save as` option and make sure that the format is called `Comma Separated Value` or `csv`
7.  Enter the analysis directory
8.  Run pipeline:

```
`/absolute/path/to/YAP/YAP.py -P ###### -E [MYUSERNAME@jcvi.org](mailto:myemail@email.com) -i /usr/local/scratch/MYUSERNAME/MYDIRECTORY/data/allinfo.csv
```
where ###### is the project id used for grid accounting.

Depending on the number of samples, analysis could take anywhere from a few hours (<50 samples) to nearly 30 hours (~400 samples) or more.

## YAP Flags

YAP_MiSeq.py flags:

* Options:
    * -h, --help            show this help message and exit
    * -H HEAD, --head=HEAD  For dry runs, import only # of lines from the input
                        files
* Required: Will not run without these !
    * -P #, --PROJECT=#   project code
    * -E @, --EMAIL=@     e-mail address
    * -i allinfo.csv, --info=allinfo.csv
                        mapping: file, barcode, primer, sample information.
                        File should be in CSV format
* Optional Configuration: parameters to alter if necessary
    * -Y #, --Yap=#       Which Pipeline: 16S ITS [16S]
    * -D, --dynamic       If specified, alignment will be scanned for primer
                        locations and trimmed accordingly. Otherwise a
                        database of known primers and trimming points will be
                        used. [False]
    * -d #, --thresh=#    in conjunction with -D, otherwise this is ignored.
                        This allows to specify how much of the alignment to
                        keep using the per-base coverage. The [0.75] value
                        indicates that ends of the alignment are trimmed until
                        a base has a coverage of [0.75] * peak coverage.
    * -a annotations, --annotations=annotations
                        directory that stores auxilliary files
                        [/usr/local/devel/ANNOTATION/sszpakow/ANNOTATION/]
    * -S #, --SAMPLE=#    perform sub.sampling of all reads based on the number
                        of reads in smallest group. if 0 - all reads are used.
                        if 1 - the sampling will be performed once, if 2 or
                        more, then 2 or more independent samplings are going
                        to be performed. [0]
    * -m #, --minlen=#    what is the minimum length of reads to process [200]
    * -g #, --mingroupsize=#
                        after demultiplexing, discard groups with fewer reads
                        than # [100]
    * -Q #, --minqual=#   Keep stretches of reads this good or better # [25]
    * -q, --quick         If specified, only single, reference DB based chimera
                        checking will be used. [False]
    * -x #, --strict=#    how strict to be at pre-clustering:  1 very strict,
                        conservative denoising (precluster identical
                        sequences)  2 less strict, aggresive denoising
                        (precluster using 98% similarity) [2]
                        
* Technical: could be useful sometimes
    * -C #, --NODESIZE=#  maximum number of grid node's CPUs to use.

# YAP: Methods and References

### Blurb

The 16S sequence-processing pipeline used for this study consists of a myriad of state-of-the-art tools carefully selected based on their literature-documented accuracy robustness and speed attributes. 

- SolexaQA trimming
- FLASH overlap
- NUCMER for primer removal

Subsequent screen.seqs function of MOTHUR was used to remove sequences shorter than **220 **bases (Engelbrektson et al., 2010). Furthermore CD-HIT (**version 4.6.1**) (Huang, Niu, Gao, Fu, & Li, 2010; Niu, Fu, Sun, & Li, 2010) was used to collapse duplicate reads and potential 454 hybrids, while retaining their count for subsequent enrichment statistics, (analogously to the functionality of unique.seqs of MOTHUR but orders of magnitude faster and less demanding on the computer hardware). At this point the sequences were aligned against **SILVA database of 16S** sequences (Pruesse et al., 2007; Schloss et al., 2009) to verify 1) the **16S** origin of noise-filtered sequences; and 2) the correct positioning of the reads with respect to the expectation of which variable regions should have been amplified and sequenced. Thereafter, the remaining sequences were subjected to MOTHUR's implementation of chimera detections (Haas et al., 2011; Schloss et al., 2009) to filter out hybrid reads. The chimeras were identified in two ways: 1) against the reads and 2) against the **SILVA database**. The alignment was trimmed to remove ends with less than **90%** of max coverage and to remove any remaining sequences that could not be aligned reliably. Taxonomical classification of the final set of reads was performed using MOTHUR&rsquo;s version of Bayesian classifier (Wang, Garrity, Tiedje, & Cole, 2007) using normalized RDP training dataset (Cole et al., 2009). The final step of the pipeline clustered the sequences based on their similarity to produce operational taxonomic units (OTU). Customarily, a similarity threshold of 97% has been used to identify OTUs at approximately the species level (see Hamady &amp; Knight, 2009). A module of CD-HIT suite (Huang et al., 2010) called CD-HIT-EST was employed to perform species-level read-clustering for subsequent analyses.

The orchestration and automation of steps has been achieved using a custom set of in-house utilities written in python and R languages known as YAP ([github.com/shpakoo/yap](http://github.com/shpakoo/yap)). JCVI grid infrastructure based on Sun Grid Engine (SGE) was used for all steps described. Enrichment and diversity statistics were calculated within MOTHUR. Further statistical analyses were accomplished using R programming language and specific R packages (ade4 for PCA).

### References

*   Chou, H. H., &amp; Holmes, M. H. (2001). DNA sequence quality trimming and vector removal. Bioinformatics (Oxford, England), 17(12), 1093;1104.
*   Cole, J. R., Wang, Q., Cardenas, E., Fish, J., Chai, B., Farris, R. J., et al. (2009). The Ribosomal Database Project: improved alignments and new tools for rRNA analysis. Nucleic acids research, 37(Database issue), D141. <a>doi:10.1093/nar/gkn879</a>
*   Engelbrektson, A., Kunin, V., Wrighton, K. C., Zvenigorodsky, N., Chen, F., Ochman, H., &amp; Hugenholtz, P. (2010). Experimental factors affecting PCR-based estimates of microbial species richness and evenness. The ISME Journal, 4(5), 642&ndash;647. <a>doi:10.1038/ismej.2009.153</a>
*   Haas, B. J., Gevers, D., Earl, A. M., Feldgarden, M., Ward, D. V., Giannoukos, G., et al. (2011). Chimeric 16S rRNA sequence formation and detection in Sanger and 454-pyrosequenced PCR amplicons. Genome research, 21(3), 494&ndash;504. <a>doi:10.1101/gr.112730.110</a>
*   Hamady, M., &amp; Knight, R. (2009). Microbial community profiling for human microbiome projects: Tools, techniques, and challenges. Genome research, 19(7), 1141&ndash;1152. <a>doi:10.1101/gr.085464.108</a>
*   Huang, Y., Niu, B., Gao, Y., Fu, L., &amp; Li, W. (2010). CD-HIT Suite: a web server for clustering and comparing biological sequences. Bioinformatics (Oxford, England), 26(5), 680&ndash;682. <a>doi:10.1093/bioinformatics/btq003</a>
*   Kunin, V., Engelbrektson, A., Ochman, H., &amp; Hugenholtz, P. (2010). Wrinkles in the rare biosphere: pyrosequencing errors can lead to artificial inflation of diversity estimates. Environmental Microbiology, 12(1), 118&ndash;123. <a>doi:10.1111/j.1462-2920.2009.02051.x</a>
*   Niu, B., Fu, L., Sun, S., &amp; Li, W. (2010). Artificial and natural duplicates in pyrosequencing reads of metagenomic data. BMC Bioinformatics, 11, 187. <a>doi:10.1186/1471-2105-11-187</a>
*   Pruesse, E., Quast, C., Knittel, K., Fuchs, B. M., Ludwig, W., Peplies, J., &amp; Gl&ouml;ckner, F. O. (2007). SILVA: a comprehensive online resource for quality checked and aligned ribosomal RNA sequence data compatible with ARB. Nucleic acids research, 35(21), 7188&ndash;7196. <a>doi:10.1093/nar/gkm864</a>
*   Quince, C., Lanz&eacute;n, A., Davenport, R. J., &amp; Turnbaugh, P. J. (2011). Removing noise from pyrosequenced amplicons. BMC Bioinformatics, 12, 38. <a>doi:10.1186/1471-2105-12-38</a>
*   Schloss, P. D., Gevers, D., &amp; Westcott, S. L. (2011). Reducing the effects of PCR amplification and sequencing artifacts on 16S rRNA-based studies. PLoS ONE, 6(12), e27310. <a>doi:10.1371/journal.pone.0027310</a>
*   Schloss, P. D., Westcott, S. L., Ryabin, T., Hall, J. R., Hartmann, M., Hollister, E. B., et al. (2009). Introducing mothur: open-source, platform-independent, community-supported software for describing and comparing microbial communities. Applied and environmental microbiology, 75(23), 7537&ndash;7541. <a>doi:10.1128/AEM.01541-09</a>
*   Wang, Q., Garrity, G. M., Tiedje, J. M., &amp; Cole, J. R. (2007). Naive Bayesian classifier for rapid assignment of rRNA sequences into the new bacterial taxonomy. Applied and environmental microbiology, 73(16), 5261&ndash;5267. <a>doi:10.1128/AEM.00062-07</a>
