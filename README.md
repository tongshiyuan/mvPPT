# mvPPT

mvPPT website is at: [http://mvppt.fudan.edu.cn/](http://mvppt.fudan.edu.cn/)  
or you can get scores at [google drive](https://drive.google.com/file/d/1zDT1e4B_-hQs4i-BLzXlcOWjCBdl7kjs/view?usp=sharing). 
A comprehensive prediction tool, mvPPT (Pathogenicity Prediction Tool for missense variants).

## Training data
Three training sets based on different combinations of variants from [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/), 
[HGMD](http://www.hgmd.cf.ac.uk/), [Uniprot](https://www.uniprot.org), 
and [Genome Aggregation Database (gnomAD)](https://gnomad.broadinstitute.org).  

## Annotation
Variants were annotated by the [ANNOVAR](https://doc-openbio.readthedocs.io/projects/annovar/en/latest/).  
- A) AFs, AAFs, and GFs of variants estimated from 125,748 exomes in gnomAD (version 2.1.1); 
- B) genomic context of the variant, i.e., region/gene-based information from [GeVIR](https://www.gevirank.org), 
[VIRLoF](https://gnomad.broadinstitute.org), [oe mis upper](https://gnomad.broadinstitute.org), 
[HIP](https://decipher.sanger.ac.uk/about/downloads/data), [CCRs](https://s3.us-east-2.amazonaws.com/ccrs/ccr.html), 
Interpro domain, and amino acid sequence before and after mutation; 
- C) pathogenicity likelihood scores assessed by different component tools, including MutationAssessor, SIFT, PROVEAN, GERP++ RS, phyloP, phastCons, and SiPhy.

We annotated datasets with ANNOVAR using dbNSFP (v.4.1a, see URLs) to generate some of the required prediction scores from different component tools, including Interpro domain, MutationAssessor, phyloP, GERP, phastCons, PROVEAN, and SiPhy. Mutations located in the interpro domains were recorded as 1 and the rest were recorded as 0. AFs, GFs, and AAFs of each variant in different populations were obtained from the gnomAD exomes database. AFs, AAFs, HomFs, and HetFs were assigned 0 and WtFs were assigned 1 if the variant was not present in the database. The GeVIR, VIRLoF, oe mis upper, HIP, and CCR scores were downloaded from their respective websites (see URLs). One-hot encoding has been applied to amino acid sequence, representing each amino acid with a binary vector of length 20 with a single non-zero value. All the features were selected to provide complementary information, and they either did not require training or their training data are publicly available to allow exclusion from our data.

The MVP, REVEL, PrimateAI, FATHMM-XF, ClinPred, MetaSVM/MetaLR, PolyPhen2, and VEST4 scores were obtained from dbNSFP v4.1a. The 
[M-CAP](http://bejerano.stanford.edu/MCAP/) (version 1.4), [MISTIC](http://lbgi.fr/mistic), [CAPICE](https://zenodo.org/record/3928295#.YFRaGi21FpQ)
[ReVe](http://varcards.biols.ac.cn),  and 
[CADD](https://cadd.gs.washington.edu/) (version 1.6) scores were downloaded from their respective websites.

## Training
mvPPT was trained using the python package [LightGBM](https://github.com/microsoft/LightGBM) (version 2.3.1), and 
parameters were tuned by [Bayesian optimization](https://github.com/fmfn/BayesianOptimization)(version 1.2.0). The random status was 
set as `1` throughout the model training process.

## Environments
The environments of mvPPT built in our study: 
- python 3.7.4
- sklearn 0.22.1
- numpy 1.17.3
- scipy 1.4.1
- pandas 0.25.3
- matplotlib 3.1.2
- lightGBM 2.3.1
- bayesian-optimization 1.1.0

