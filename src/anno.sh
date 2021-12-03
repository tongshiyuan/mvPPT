#!/bin/bash
file=$1
tag=`basename $file .vcf`
perl ../annovar/convert2annovar.pl -format vcf4 $file > ${file}.avinput
perl ../annovar/table_annovar.pl \
	${file}.avinput ../annodb --buildver hg19 \
	-remove -out $tag \
	-protocol ensGene,CCRs,HIP,gnomad211exoms_allpop,dbnsfp35a,p6b \
	-operation g,r,r,f,f,f -nastring . --thread 12

awk -F "\t" '{if ($1=="Chr" || $9== "nonsynonymous SNV") {print}}' ${tag}.hg19_multianno.txt > ${tag}.mssnv.txt
awk -F "\t" '{if($7!~/^E.*;.*/){print}}' ${tag}.mssnv.txt > ${tag}GeneUniqAnno.txt

rm ${file}.avinput ${tag}.hg19_multianno.txt ${tag}.mssnv.txt

