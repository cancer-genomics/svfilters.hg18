#!/bin/bash
#module load samtools
#module load bedtools
#module load perl
R_LIBS_USER=~/Library/R/3.3-bioc-release R-bioc-stable --vanilla CMD BATCH --no-save binAssemblyGaps.R
