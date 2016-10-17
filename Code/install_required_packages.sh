#!/bin/bash

echo "Downloading & installing package biostrings"
Rscript biostrings.R

echo "Downloading tRap package"
wget -nd -nH -P ${HOME}/Downloads http://trap.molgen.mpg.de/download/TRAP_R_package/tRap_0.5.tar.gz

echo "Downloading evd package"
wget -nd -nH -P ${HOME}/Downloads https://cran.rstudio.com/src/contrib/evd_2.3-2.tar.gz

EVDPCK=${HOME}/Downloads/evd_2.3-2.tar.gz
TRAPPCK=${HOME}/Downloads/tRap_0.5.tar.gz

echo "Installing evd package"
R CMD INSTALL ${EVDPCK}

echo "Installing tRap package"
R CMD INSTALL ${TRAPPCK}

rm ${HOME}/Downloads/evd_2.3-2.tar.gz
rm ${HOME}/Downloads/tRap_0.5.tar.gz
