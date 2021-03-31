import pandas as pd


def counts(df, tag):
    counts = df['gene_id'].value_counts()
    dfc = pd.DataFrame(counts)
    dfc = dfc.reset_index()
    dfc.rename(columns={'gene_id': tag}, inplace=True)
    return dfc


def cutoff(data, df, th, tag):
    geneList = df[df[tag] > th]['index'].tolist()
    geneIndex = df[df[tag] > th].index.tolist()
    df2 = pd.DataFrame(columns=data.columns.tolist())
    for i, j in zip(geneList, geneIndex):
        dftmp = data[data.gene_id == i]
        dfds = dftmp.sample(n=int(th * df['length'][j]))  # 抽取80%样本的比例
        df2 = df2.append(dfds)
    df2_other = data[~data.gene_id.isin(geneList)]
    df2 = df2.append(df2_other)
    return df2


def downsampling(df, rate=0.2, th=[0, 0], file='ensCanonicalTranscript.txt', nos=False,
                 keep=False):
    dfens = pd.read_table(file, low_memory=False, header=None)
    dfens.rename(columns={0: 'index', 2: 'length'}, inplace=True)

    dfVarAAF_0 = df[df.tags == 0]
    dfVarAAF_1 = df[df.tags == 1]
    dfVarAAF_1c = counts(dfVarAAF_1, 'Pathogenic')
    dfVarAAF_0c = counts(dfVarAAF_0, 'Benign')
    dfc = pd.merge(dfVarAAF_1c, dfVarAAF_0c, how='outer', on='index')
    dfMerge = pd.merge(dfc, dfens, how='left', on='index')
    dfMerge.fillna(0, inplace=True)
    if nos and (not keep):
        keepList = dfMerge[dfMerge['Pathogenic'] > 1]['index'].tolist()
        dfVarAAF_1 = dfVarAAF_1[dfVarAAF_1['gene_id'].isin(keepList)]
    elif nos and keep:
        keepList = dfMerge[dfMerge['Pathogenic'] > 1]['index'].tolist()
        dfVarAAF_1 = dfVarAAF_1[dfVarAAF_1['gene_id'].isin(keepList)]
        df2 = pd.concat([dfVarAAF_0, dfVarAAF_1])
        return df2

    dfMerge['BR'] = dfMerge['Benign'] / dfMerge['length']
    dfMerge['PR'] = dfMerge['Pathogenic'] / dfMerge['length']
    dfB = dfMerge[dfMerge['BR'] > 0]
    dfP = dfMerge[dfMerge['PR'] > 0]
    dfB = dfB.sort_values(['BR'], ascending=False)
    dfB.index = pd.RangeIndex(len(dfB.index))
    dfP = dfP.sort_values(['PR'], ascending=False)
    dfP.index = pd.RangeIndex(len(dfP.index))
    if rate:
        thB = dfB['BR'][int(dfB.shape[0] * rate)]
        thP = dfP['PR'][int(dfP.shape[0] * rate)]
    else:
        thB = th[0]
        thP = th[1]
    dsB = cutoff(dfVarAAF_0, dfB, thB, 'BR')
    dsP = cutoff(dfVarAAF_1, dfP, thP, 'PR')

    df2 = pd.concat([dsB, dsP])

    return df2
