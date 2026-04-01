rm -rf ~/temp/INLAtools*
mkdir ~/temp/INLAtools/
cp -r DESCRIPTION NAMESPACE R/ man/ src/ demo/ ~/temp/INLAtools/
cd vignettes/
mkdir ~/temp/INLAtools/vignettes/
cp cgenericKronecker.pdf cgenericKronecker.pdf.asis ~/temp/INLAtools/vignettes/
cd ~/temp/INLAtools/vignettes/
R -e "tools::compactPDF('cgenericKronecker.pdf', gs_quality = 'ebook')"
cd ../../
R CMD build INLAtools
R CMD check INLAtools_*.tar.gz --as-cran


