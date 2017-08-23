#!/bin/sh

COUNTER=0
SCRIPT_NAME=$1

while [ $COUNTER -lt 14 ]; do
    python2 $SCRIPT_NAME $COUNTER
    echo $COUNTER
    let COUNTER=COUNTER+1
done
