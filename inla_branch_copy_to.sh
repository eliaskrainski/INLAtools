cd src/
cp cgeneric_generic0.c cgeneric_kronecker.c INLAtools.h ../../inla_branch_INLAtools/src/
cd ../inla_branch_INLAtools/
git add src/* 
git commit -m 'update inla branch'
git push
cd ../INLAtools/
