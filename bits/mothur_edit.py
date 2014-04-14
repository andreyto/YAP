#!/usr/bin/env python
import logging_helper
from subprocess import check_call
import sh
import sys, os, glob, shutil, tempfile
from contextlib import nested
from os.path import join as pjoin

from argh import ArghParser,arg

log = logging.getLogger(os.path.basename(sys.argv[0]))


def remove_missing_seqs(list_file,group_file):
    """Remove from group file sequence names not found in list file"""
    pass

def main():
    logging_helper.logging_config(detail="high")
    parser = ArghParser()
    parser.add_commands([
        remove_missing_seqs
        ]
        )
    parser.dispatch()

if __name__ == "__main__":
    main()

