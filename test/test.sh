#!/bin/bash

set -o errexit
set -o nounset
set -o pipefail

# Arguments are:
# -c Specify Chothia datafile (Default: chothia.dat)
# -v Verbose; give explanations when no canonical found
# -n The sequence file has Chothia (rather than Kabat) numbering
    
rm -f ./test?.out

../chothia -c ./chothia.dat.ex1 -v ./numbered.kabat.dat > test1.out 2>&1 
../chothia -c ./chothia.dat.ex2 -v ./numbered.kabat.dat > test2.out 2>&1 
../chothia -c ./chothia.dat.ex3 -v ./numbered.kabat.dat > test3.out 2>&1 

echo "chothia tests passed"

