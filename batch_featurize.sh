#!/bin/bash

CONFS_DIR=/path/to/dir
ETSC_DIR=/other/path/to/dir

for file in $CONFS_DIR/*.py; do
  python2 $ETSC_DIR/featurize.py $file
done
