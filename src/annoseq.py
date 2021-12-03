import os
import sys
import pandas as pd

infile = sys.argv[1]
outfile = sys.argv[2]

dfGene = pd.read_table('../annodb/GeVIR_Gnomad.txt')
dfClin = pd.read_table(infile, low_memory=False)
dfClin.rename(columns={'Gene.ensGene': 'gene_id'}, inplace=True)
dfAnno = pd.merge(dfClin, dfGene, how='left', on='gene_id')
dfAnno.to_csv('anno.tmp1', index=False, sep='\t')

canonicalDict = {}
with open('../annodb/ens_long_transcript.txt') as f:
    for i in f:
        canonicalDict[i.split('\t')[0]] = i.split('\t')[1].strip()

of = open('anno.tmp2', 'w')
with open('anno.tmp1') as f:
    head = f.readline()
    of.write(head.rstrip() + '\ttranscript_select\tcanonical_AAChange.ensGene\n')
    # canindex =head.split('\t').index('canonical_transcript')
    excindex = head.split('\t').index('AAChange.ensGene')
    for i in f:
        line = i.split('\t')
        trs = line[excindex].split(',')
        for j in trs:
            if j.split(':')[1] == canonicalDict[j.split(':')[0]]:  # line[canindex]:
                of.write(i.rstrip() + '\t' + j.split(':')[1] + '\t' + j + '\n')
                break
        else:
            of.write(i.rstrip() + '\t' + trs[0].split(':')[1] + '\t' + trs[0] + '\n')
of.close()

dfEns = pd.read_table('../annodb/cds.seq', low_memory=False, header=None)
dfEns.rename(columns={0: 'transcript_select', 1: 'CDSSEQ'}, inplace=True)
dfClin = pd.read_table('anno.tmp2', low_memory=False)
dfMerge = pd.merge(dfClin, dfEns, how='left', on='transcript_select')
dfPep = pd.read_table('../annodb/ensGeneAA.seq', low_memory=False, header=None)
dfPep.rename(columns={0: 'transcript_select', 1: 'AASEQ'}, inplace=True)
dfMerge2 = pd.merge(dfMerge, dfPep, how='left', on='transcript_select')

dfMerge2.to_csv(outfile, sep='\t', index=False)
os.remove('anno.tmp1')
os.remove('anno.tmp2')