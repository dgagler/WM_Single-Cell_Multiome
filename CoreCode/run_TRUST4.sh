#!/bin/bash

# Script for running TRUST4 on each patient individually using their bams. Note that you must supply a new patient bam for each run or integrate this into a smarter workflow (sorry lol)
# To be run from within the TRUST4 directory which contains executables
# For more info on how to run, see the TRUST4 github at https://github.com/liulab-dfci/TRUST4

./run-trust4 \
-f /Users/gagled01/morganLab/Waldenstroms/TRUST4_working/TRUST4/hg38_bcrtcr.fa \
--ref /Users/gagled01/morganLab/Waldenstroms/TRUST4_working/TRUST4/human_IMGT+C.fa \
-b /Users/gagled01/morganLab/Waldenstroms/TRUST4_working/bams/ExamplePatient_bam.bam \
--barcode CB \
--barcodeRange 0 15 + \
--read1Range 16 -1