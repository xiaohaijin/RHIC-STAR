#! /bin/bash

mkdir blank
rsync -a --delete blank/ ./err
rsync -a --delete blank/ ./out
rsync -a --delete blank/ ./list
rsync -a --delete blank/ ./report
rsync -a --delete blank/ ./csh
rm -rf blank
