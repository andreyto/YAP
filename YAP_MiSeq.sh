#!/bin/bash

this_dir=$(dirname $0)

. "$this_dir/YAP.rc"

#screen $YAP_DEPS/bin/python $(which rpdb2) -s $this_dir/YAP_MiSeq.py "$@"
#exec screen rpdb2 -s $this_dir/YAP_MiSeq.py "$@"
exec $YAP_DEPS/python $YAP_SCRIPTS/YAP_MiSeq.py "$@"

