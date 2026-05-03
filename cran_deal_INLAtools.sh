rm -rf ~/temp/INLAtools*
mkdir ~/temp/INLAtools/
cp -r DESCRIPTION NAMESPACE R/ man/ demo/ ~/temp/INLAtools/
mkdir ~/temp/INLAtools/src/
cd src
cp INLAtools.h INLAtools_init.c cgeneric.h cgeneric_generic0.c cgeneric_kronecker.c cgeneric_get.c Makevars ~/temp/INLAtools/src/
#cd vignettes/
#mkdir ~/temp/INLAtools/vignettes/
#cp cgenericKronecker.pdf cgenericKronecker.pdf.asis ~/temp/INLAtools/vignettes/
#cd ~/temp/INLAtools/vignettes/
#R -e "tools::compactPDF('cgenericKronecker.pdf', gs_quality = 'ebook')"
#cd ../../
cd ~/temp/
R CMD build INLAtools
R CMD check INLAtools_*.tar.gz --as-cran


