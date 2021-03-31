import pandas as pd
import matplotlib.pyplot as plt


def getCount(file):
    df = pd.read_table(file, low_memory=False)
    df = df[df['ExonicFunc.ensGene'] == 'nonsynonymous SNV']
    counts = df['canonical_transcript'].value_counts()
    df = pd.DataFrame(counts)
    df = df.reset_index()
    df.rename(columns={'canonical_transcript': file.split('.')[0]}, inplace=True)


dfGnomadExomes = getCount('gnomadExomes.txt')  # annotated by selectCanonicalTranscript.py, same as below
dfGnomaGenomes = getCount('gnomadGenomes.txt')
dfKGP = getCount('kgp.txt')
dfClinP = getCount('ClinVarPathogenic.txt')
dfClinB = getCount('ClinVarBenign.txt')
dfHGMD = getCount('HGMD.txt')
dfUniprotP = pd.read_table('uniprotPathogenic.txt', low_memory=False)
dfUniprotB = pd.read_table('uniprotBenign.txt', low_memory=False)

dfCounts = pd.merge(dfGnomadExomes, dfGnomaGenomes, how='outer', on='index')
dfCounts = pd.merge(dfCounts, dfKGP, how='outer', on='index')
dfCounts = pd.merge(dfCounts, dfClinP, how='outer', on='index')
dfCounts = pd.merge(dfCounts, dfClinB, how='outer', on='index')
dfCounts = pd.merge(dfCounts, dfHGMD, how='outer', on='index')
dfCounts = pd.merge(dfCounts, dfUniprotP, how='outer', on='index')
dfCounts = pd.merge(dfCounts, dfUniprotB, how='outer', on='index')

dfEns = pd.read_table('ensCanonicalTranscript.txt', low_memory=False, header=None)
dfEns.rename(columns={0: 'gene', 1: 'index', 2: 'length'}, inplace=True)
dfMerge = pd.merge(dfCounts, dfEns, how='left', on='index')
dfMerge.dropna(subset=['length'], inplace=True)
dfMerge.fillna(0, inplace=True)

dfCounts['gnomadExomes'] = dfCounts['gnomadExomes'] / dfCounts['length']
dfCounts['gnomadGenomes'] = dfCounts['gnomadGenomes'] / dfCounts['length']
dfCounts['kpg'] = dfCounts['kpg'] / dfCounts['length']
dfCounts['ClinVarPathogenic'] = dfCounts['ClinVarPathogenic'] / dfCounts['length']
dfCounts['ClinVarBenign'] = dfCounts['ClinVarBenign'] / dfCounts['length']
dfCounts['HGMD'] = dfCounts['HGMD'] / dfCounts['length']
dfCounts['UniprotPathogenic'] = dfCounts['UniprotPathogenic'] / dfCounts['length']
dfCounts['UniprotBenign'] = dfCounts['UniprotBenign'] / dfCounts['length']

df = dfCounts.sort_values(['length'], ascending=False)
df['length'] = df['length'] / max(df['length'])
df.index = pd.RangeIndex(len(df.index))

plt.figure(figsize=(24, 16))
plt.plot(df.index, df['gnomadExomes'], alpha=0.8, label='gnomadExomes')
plt.plot(df.index, df['gnomadGenomes'], alpha=0.8, label='gnomadGenomes')
plt.plot(df.index, df['kpg'], alpha=0.8, label='kpg')
plt.plot(df.index, df['length'], alpha=1, label='length')
plt.legend(loc="upper right")

plt.figure(figsize=(24, 16))
plt.plot(df.index, df['HGMD'], alpha=0.8, label='HGMD')
plt.plot(df.index, df['UniprotPathogenic'], alpha=0.8, label='UniprotPathogenic')
plt.plot(df.index, df['ClinVarPathogenic'], alpha=0.8, label='ClinVarPathogenic')
plt.plot(df.index, df['length'], alpha=1, label='length')
plt.legend(loc="upper right")

plt.figure(figsize=(24, 16))
plt.plot(df.index, df['ClinVarBenign'], alpha=0.8, label='ClinVarBenign')
plt.plot(df.index, df['UniprotBenign'], alpha=0.8, label='UniprotBenign')
plt.plot(df.index, df['length'], alpha=1, label='length')
plt.legend(loc="upper right")
