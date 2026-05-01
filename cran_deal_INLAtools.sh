rm -rf ~/temp/INLAtools*
mkdir ~/temp/INLAtools/
cp -r DESCRIPTION NAMESPACE R/ man/ demo/ ~/temp/INLAtools/
mkdir ~/temp/INLAtools/src/
cp src/INLAtools.h src/INLAtools_init.c src/cgeneric.h src/cgeneric_get.c src/Makevars ~/temp/INLAtools/src/
#cd vignettes/
#mkdir ~/temp/INLAtools/vignettes/
#cp cgenericKronecker.pdf cgenericKronecker.pdf.asis ~/temp/INLAtools/vignettes/
#cd ~/temp/INLAtools/vignettes/
#R -e "tools::compactPDF('cgenericKronecker.pdf', gs_quality = 'ebook')"
#cd ../../
cd ~/temp/
R CMD build INLAtools
R CMD check INLAtools_*.tar.gz --as-cran


