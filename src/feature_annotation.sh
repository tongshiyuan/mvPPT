#!/bin/bash
# note: CCRs,GeVIR_per,LOEUF_per,VIRLoF_per,HIP_score,gnomad211exoms_allpop,ReVe,cadd14_snp,clinpred,mcap14,PrimateAI was built by ourselves
# annotate training set demo
perl table_annovar.pl ../example/var1_demo.txt /path/of/annovar/database/ --buildver hg19 -remove -out var1_demo \
-protocol refGene,CCRs,GeVIR_per,LOEUF_per,VIRLoF_per,HIP_score,gnomad211exoms_allpop,dbnsfp35a,ReVe,cadd14_snp,clinpred,mcap14,PrimateAI \
-operation g,r,r,r,r,r,f,f,f,f,f,f,f -nastring .
perl table_annovar.pl ../example/var2_demo.txt /path/of/annovar/database/ --buildver hg19 -remove -out var2_demo \
-protocol refGene,CCRs,GeVIR_per,LOEUF_per,VIRLoF_per,HIP_score,gnomad211exoms_allpop,dbnsfp35a,ReVe,cadd14_snp,clinpred,mcap14,PrimateAI \
-operation g,r,r,r,r,r,f,f,f,f,f,f,f -nastring .
perl table_annovar.pl ../example/var3_demo.txt /path/of/annovar/database/ --buildver hg19 -remove -out var3_demo \
-protocol refGene,CCRs,GeVIR_per,LOEUF_per,VIRLoF_per,HIP_score,gnomad211exoms_allpop,dbnsfp35a,ReVe,cadd14_snp,clinpred,mcap14,PrimateAI \
-operation g,r,r,r,r,r,f,f,f,f,f,f,f -nastring .
cut -f 6 ../example/var3_demo.txt > ../example/var3_demo.tags
sed -i '2d' ../example/var3_demo.hg19_multianno.txt
paste ../example/var3_demo.hg19_multianno.txt ../example/var3_demo.tags > ../example/var3_demo.anno.txt
# annotation test set demo
perl table_annovar.pl ../example/ClinHGMD_demo.txt /path/of/annovar/database/ --buildver hg19 -remove -out ClinHGMD_demo \
-protocol refGene,CCRs,GeVIR_per,LOEUF_per,VIRLoF_per,HIP_score,gnomad211exoms_allpop,dbnsfp35a,ReVe,cadd14_snp,clinpred,mcap14,PrimateAI \
-operation g,r,r,r,r,r,f,f,f,f,f,f,f -nastring .
cut -f 6 ../example/ClinHGMD_demo.txt > ../example/ClinHGMD_demo.tags
sed -i '2d' ../example/ClinHGMD_demo.hg19_multianno.txt
paste ../example/ClinHGMD_demo.hg19_multianno.txt ../example/ClinHGMD_demo.tags > ../example/ClinHGMD_demo.anno.txt