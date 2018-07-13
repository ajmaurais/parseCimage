
#/bin/bash

rsync -vr -e 'ssh -p 22022' --exclude='.*' --exclude='*.pyc' src mauraisa@pleiades.bc.edu:~/scripts/parseCimage
