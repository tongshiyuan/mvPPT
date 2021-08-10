#!/bin/bash
file=$1
tag=`basename $file .txt`
perl /home/biowork/biodatabase/annotation/annovar/table_annovar.pl \
	$file /bio-analysis3/mvPPT/version_2/annodb --buildver hg19 \
	-remove -out $tag \
	-protocol ensGene,CCRs,HIP,gnomad211exoms_allpop,dbnsfp35a,p6b,CADD16,PrimateAI,capice,ClinPred,fathmmxf,mcap14,mistic,ReVe,VEST4,MVP \
	-operation g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . --thread 12
cut -f 6- $file > ${tag}.tags
sed -i '2d' ${tag}.hg19_multianno.txt
paste ${tag}.hg19_multianno.txt ${tag}.tags > ${tag}.anno.txt

rm ${tag}.hg19_multianno.txt ${tag}.tags