#!/bin/bash

YAP_DEPS=/usr/local/devel/ANNOTATION/sszpakow/YAP

this_dir=$(cd $(dirname $0) && pwd)

export PATH=$YAP_DEPS/bin:$PATH

# make sure all the libraries are linked
export LD_LIBRARY_PATH=/usr/local/packages/mysql/lib:/usr/local/packages/gcc/lib64:/usr/local/packages/curl/lib:$LD_LIBRARY_PATH

ulimit -n hard

#screen $YAP_DEPS/bin/python $(which rpdb2) -s $this_dir/YAP_454.py "$@"
#exec screen rpdb2 -s $this_dir/YAP_454.py "$@"
$YAP_DEPS/bin/python $this_dir/YAP_454.py "$@"

