import sys
import pickle
import pandas as pd

model = pickle.load(open('../modelFile/singleTon.m', 'rb'))

trainFeatures = [
    'Gnomad_exomes_AF', 'Gnomad_exomes_wtf', 'Gnomad_exomes_hetf', 'Gnomad_exomes_homf', 'Gnomad_AAF',
    'MutationAssessor_score', 'SIFT_score', 'PROVEAN_score', 'GERP++_RS',
    'phyloP100way_vertebrate.1', 'phyloP30way_mammalian', 'phyloP17way_primate',
    'phastCons100way_vertebrate.1', 'phastCons30way_mammalian', 'phastCons17way_primate', 'SiPhy_29way_logOdds',
    'CCRs', 'Interpro_domain', 'HIP', 'oe_mis_upper', 'gevir_percentile', 'virlof_percentile',
    'refAA_A', 'refAA_C', 'refAA_D', 'refAA_E', 'refAA_F', 'refAA_G', 'refAA_H', 'refAA_I', 'refAA_K', 'refAA_L',
    'refAA_M', 'refAA_N', 'refAA_P', 'refAA_Q', 'refAA_R', 'refAA_S', 'refAA_T', 'refAA_V', 'refAA_W', 'refAA_Y',
    'altAA_A', 'altAA_C', 'altAA_D', 'altAA_E', 'altAA_F', 'altAA_G', 'altAA_H', 'altAA_I', 'altAA_K', 'altAA_L',
    'altAA_M', 'altAA_N', 'altAA_P', 'altAA_Q', 'altAA_R', 'altAA_S', 'altAA_T', 'altAA_V', 'altAA_W', 'altAA_Y']


def getAA(df):
    # 氨基酸信息
    aachange_col = list(df['canonical_AAChange.ensGene'])
    ref = []
    alt = []
    pos = []
    aa = []
    for i, v in enumerate(aachange_col):
        ref.append(v.split(':')[4][2])
        alt.append(v.split(':')[4][-1])
        pos.append(v.split(':')[4][3:-1])
        aa.append(v.split(':')[4])
    df['refAA'] = ref
    df['altAA'] = alt
    df['posAA'] = pos
    df['canonical_AAChange'] = aa
    return df


def feature_edit(data):
    for k, v in data.items():
        if k in ['CCRs', 'HIP']:
            my_column = list(data[k])
            for i, val in enumerate(my_column):
                if isinstance(val, str):
                    vals = val.lstrip('Name=').split(',')
                    removed = ['NA', '.']
                    for item in removed:
                        if item in vals:
                            vals.remove(item)
                    if vals:
                        val = max(vals)
                    else:
                        val = ''
                    my_column[i] = val
            data[k] = my_column
            data[k].replace('', pd.NA, inplace=True)
    data['CCRs'].replace(pd.NA, 0, inplace=True)
    # gnomad
    for i in ['Gnomad_exomes_AF', 'Gnomad_exomes_hetf', 'Gnomad_exomes_homf']:
        data[i].replace('.', 0, inplace=True)
    for i in ['Gnomad_exomes_wtf']:
        data[i].replace('.', 1, inplace=True)
    # domain
    my_col = list(data['Interpro_domain'])
    for index, i in enumerate(my_col):
        i = i.strip('.;')
        if i:
            my_col[index] = 1
        else:
            my_col[index] = 0
    data['Interpro_domain'] = my_col
    for c in ['SIFT_score', 'GERP++_RS', 'phyloP100way_vertebrate', 'phyloP20way_mammalian',
              'phastCons100way_vertebrate', 'phastCons20way_mammalian', 'SiPhy_29way_logOdds',
              'phyloP100way_vertebrate.1', 'phyloP30way_mammalian',
              'phyloP17way_primate', 'phastCons100way_vertebrate.1', 'phastCons30way_mammalian',
              'phastCons17way_primate'
              ]:
        data[c].replace('.', pd.NA, inplace=True)
    for c in ['oe_mis', 'mu_mis', 'oe_mis_lower', 'oe_mis_upper', 'mis_z', 'gevir_percentile',
              'virlof_percentile', 'gevir_ad_enrichment', 'virlof_ad_enrichment', 'gevir_ar_enrichment',
              'virlof_ar_enrichment']:
        data[c].replace('', pd.NA, inplace=True)
        data[c].replace('.', pd.NA, inplace=True)

    ref_list = []
    refString = []
    my_col = list(data['AASEQ'])
    aaPos = list(data['posAA'])
    for index, i in enumerate(my_col):
        if int(aaPos[index]) < 5:
            refn = i[:11]
        elif len(i) - int(aaPos[index]) < 5:
            refn = i[-11:]
        else:
            refn = i[int(aaPos[index]) - 5:int(aaPos[index]) + 6]
        refn2 = refn
        for index2, j in enumerate(refn2):
            if j not in ['G', 'A', 'V', 'L', 'I', 'F', 'W', 'Y', 'D', 'N', 'E', 'K', 'Q', 'M', 'S', 'T', 'C', 'P', 'H',
                         'R']:
                # refn=refn.replace(j,' ')
                refn = refn[:index2] + ' ' + refn[index2 + 1:]
        ref_list.append(list(refn))
        refString.append(refn)

    df = pd.DataFrame(ref_list)
    df.replace(' ', pd.NA, inplace=True)
    dfSEQ = pd.get_dummies(df)

    dfseq = pd.DataFrame(refString, columns=['AAString'])

    dummies = pd.get_dummies(data[['refAA', 'altAA']])
    dummies.index = pd.RangeIndex(len(dummies.index))
    data.index = pd.RangeIndex(len(data.index))
    dfSEQ.index = pd.RangeIndex(len(dfSEQ.index))
    # print(data.shape, dummies.shape, dfSEQ.shape)
    dfMerge = pd.concat([data, dummies, dfSEQ, dfseq], axis=1)
    # print(dfMerge.shape)
    # print(dfMerge.isnull().sum() / dfMerge.shape[0])
    return dfMerge


df = pd.read_table(sys.argv[1], low_memory=False)
dfAAF = pd.read_table('../annodb/gnomad_AAF.txt', low_memory=False, header=None)
dfAAF.rename(columns={0: 'transcript_select', 1: 'canonical_AAChange', 2: 'Gnomad_AAF'}, inplace=True)
df = df[~pd.isna(df.AASEQ)]
df = getAA(df)
df = pd.merge(df, dfAAF, how='left', on=['transcript_select', 'canonical_AAChange'])
df['Gnomad_AAF'].fillna(0, inplace=True)
df = df[df['altAA'] != '*']
df = feature_edit(df)
df = df.dropna(subset=trainFeatures)
df.index = pd.RangeIndex(len(df.index))
df = df.apply(pd.to_numeric, errors='ignore')
df['mvPPT'] = model.predict_proba(df[trainFeatures[:-1]])[:, 1]
df[['Chr', 'Start', 'End', 'Ref', 'Alt', 'mvPPT']].to_csv('result.txt', sep='\t', index=False)
