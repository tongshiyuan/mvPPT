# mvPPT

mvPPT website is at: [http://mvppt.fudan.edu.cn/](http://mvppt.fudan.edu.cn/)  
A comprehensive prediction tool, mvPPT (Pathogenicity Prediction Tool for missense variants).

## Training data
Three training sets based on different combinations of variants from [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/), 
[HGMD](http://www.hgmd.cf.ac.uk/), [1000 Genomes Project](https://www.internationalgenome.org), 
and [Genome Aggregation Database (gnomAD)](https://gnomad.broadinstitute.org).  

## Annotation
Variants were annotated by the [ANNOVAR](https://doc-openbio.readthedocs.io/projects/annovar/en/latest/).  
- A) AFs and GFs of variants estimated from 125,748 exomes in gnomAD (version 2.1.1); 
- B) genomic context of the variant, i.e., region/gene-based information from [GeVIR](https://www.gevirank.org), 
[VIRLoF](https://gnomad.broadinstitute.org), [LOEUF](https://gnomad.broadinstitute.org), 
[HIP](https://decipher.sanger.ac.uk/about/downloads/data), [CCRs](https://s3.us-east-2.amazonaws.com/ccrs/ccr.html), and 
Interpro domain; 
- C) pathogenicity likelihood scores assessed by different component tools, including MutationAssessor, fitCons, phyloP, 
GERP, phastCons, PROVEAN, SiPhy, and GenoCanyon.

We annotated datasets with ANNOVAR using dbNSFP v.3.5a to generate some of the required prediction scores from different 
component tools, including Interpro domain, MutationAssessor, fitCons, phyloP, GERP, phastCons, PROVEAN, SiPhy, and 
GenoCanyon. Mutations located in the interpro domains were recorded as 1 and the rest were recorded as 0. AFs and GFs of 
each variant in different populations were obtained from the gnomAD database. AFs and GFs were computed from different 
populations: all, African/African American (AFR), Latino/Admixed American (AMR), Ashkenazi Jewish (ASJ), East Asian 
(EAS), Finnish (FIN), Non-Finnish European (NFE), South Asian (SAS), and others (OTH). AFs, HomFs, and HetFs were 
assigned 0 and WtFs were assigned 1 if the variant was not present in the database. The GeVIR, VIRLoF, LOEUF, HIP, 
and CCR scores were downloaded from their respective websites. All the features were selected to provide complementary 
information, and they either did not require training or their training data are publicly available to allow exclusion 
from our data.

The REVEL, MetaSVM/MetaLR, SIFT, PolyPhen2, and VEST3 scores were obtained from dbNSFP v3.5a. The 
[M-CAP](http://bejerano.stanford.edu/MCAP/) (version 1.4), [ClinPred](https://sites.google.com/site/clinpred/), 
[ReVe](http://varcards.biols.ac.cn), [PrimateAI](https://basespace.illumina.com/s/cPgCSmecvhb4), and 
[CADD](https://cadd.gs.washington.edu/) (version 1.4) scores were downloaded from their respective websites.

## Training
mvPPT was trained using the python package [LightGBM](https://github.com/microsoft/LightGBM) (version 2.3.1), and 
parameters were tuned by [Bayesian optimization](https://github.com/fmfn/BayesianOptimization). The random status was 
set as `2020` throughout the model training process.

