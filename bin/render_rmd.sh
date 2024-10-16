#!/usr/bin/env bash

module purge all
module load fhR/4.1.2-foss-2021b
export R_LIBS=/fh/fast/_SR/Genomics/user/fwu/R/x86_64-pc-linux-gnu-library/4.1 # 
export PATH=/app/software/RStudio-Server/1.4.1717-foss-2021b-Java-11-R-4.1.2/bin/pandoc/:$PATH

R --vanilla "--args $@" <<code
args <- as.vector(commandArgs(T)); 
if(length(args) < 2){
    rmarkdown::render(args[1])
}else{
    rmarkdown::render(args[1], output_file= args[2])
}

code

