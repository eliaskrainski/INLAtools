rm -rf ~/temp/INLAtools*
mkdir ~/temp/INLAtools/
cp -r DESCRIPTION NAMESPACE R/ man/ src/ demo/ ~/temp/INLAtools/
# cd vignettes/
# mkdir ~/temp/INLAtools/vignettes/
# cp preamble.tex references.bib *.Rmd ~/temp/INLAtools/vignettes/
cd ~/temp/
R CMD build INLAtools
R CMD check INLAtools_*.tar.gz --as-cran


