#!/bin/bash
module load amber/openmpi/intel/20.11
mkdir ${1%.pdb}
cd ${1%.pdb}
antechamber -i ../$1 -fi pdb -o ${1/pdb/ac} -fo ac -nc ${2-0}  -c bcc

