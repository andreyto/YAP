#!/bin/bash

this_dir=$(cd $(dirname $0) && pwd)

## Environment variables that YAP code needs to be defined on submit host
export YAP_DEPS="/usr/local/projects/GATES/jshankar/YAPCOPY/sszpakow/YAP/bin"
## This will be used if --annotations argument is not passed
export YAP_DATA="/usr/local/devel/ANNOTATION/sszpakow/ANNOTATION/"
## You probably only ever change this if you are a developer
export YAP_SCRIPTS="$this_dir"

export PATH=$YAP_DEPS:/usr/local/packages/graphviz/bin:$PATH

# make sure all the libraries are linked
export LD_LIBRARY_PATH=/usr/local/packages/mysql/lib:/usr/local/packages/gcc/lib64:/usr/local/packages/curl/lib:$LD_LIBRARY_PATH

ulimit -n hard

#screen $YAP_DEPS/bin/python $(which rpdb2) -s $this_dir/YAP_MiSeq.py "$@"
#exec screen rpdb2 -s $this_dir/YAP_MiSeq.py "$@"
$YAP_DEPS/python $this_dir/YAP_MiSeq.py "$@"

