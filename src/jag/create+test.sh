#!/bin/sh
cd ..
rm  src/test/checksum.md5
rm *
cp src/jag.py jag && chmod +x jag
echo "step 1"
./jag --group src/test/total.txt --bfile=src/test/hapmap_CEU_r23_60 --pheno=src/test/pheno60.txt 

echo "step 2"
./jag --group src/test/total.txt --bfile=src/test/hapmap_CEU_r23_60 --pheno=src/test/pheno60.txt --perm=5 --seed=11
./jag --group src/test/total.txt --bfile=src/test/hapmap_CEU_r23_60 --pheno=src/test/pheno60.txt --perm=2 --seed=12
echo "step 3"
./jag --merge
echo "step 4"
./jag --draw 100 --group src/test/total.txt --neff FG2 --data=src/test/all_snps_in_gene_NCBI.txt --bfile=src/test/hapmap_CEU_r23_60 --emp empp.pheno1.out --seed=11 
./jag --draw 100 --group src/test/total.txt --ngenes FG2 --data=src/test/all_snps_in_gene_NCBI.txt --bfile=src/test/hapmap_CEU_r23_60 --emp empp.pheno1.out --seed=11

md5sum  * |grep -v jag |grep -v perm.merged.P2 > src/test/checksum.md5 
rm *
cp src/jag.py jag && chmod +x jag
cd ..

