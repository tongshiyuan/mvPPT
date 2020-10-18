import numpy as np
import pandas as pd

categoryABC = ['GeVIR_per', 'LOEUF_per', 'VIRLoF_per', 'HIP_score', 'CCRs', 'Interpro_domain', 'Gnomad_exomes_wtf',
               'Gnomad_exomes_EAS_wtf', 'Gnomad_exomes_NFE_wtf',
               'Gnomad_exomes_AMR_wtf', 'Gnomad_exomes_ASJ_wtf', 'Gnomad_exomes_FIN_wtf',
               'Gnomad_exomes_AFR_wtf', 'Gnomad_exomes_OTH_wtf', 'Gnomad_exomes_SAS_wtf',
               'Gnomad_exomes_AF', 'Gnomad_exomes_EAS_AF', 'Gnomad_exomes_NFE_AF',
               'Gnomad_exomes_AMR_AF', 'Gnomad_exomes_ASJ_AF', 'Gnomad_exomes_FIN_AF',
               'Gnomad_exomes_AFR_AF', 'Gnomad_exomes_OTH_AF', 'Gnomad_exomes_SAS_AF',
               'Gnomad_exomes_hetf', 'Gnomad_exomes_homf', 'Gnomad_exomes_EAS_hetf',
               'Gnomad_exomes_EAS_homf', 'Gnomad_exomes_NFE_hetf', 'Gnomad_exomes_NFE_homf',
               'Gnomad_exomes_AMR_hetf', 'Gnomad_exomes_AMR_homf', 'Gnomad_exomes_ASJ_hetf',
               'Gnomad_exomes_ASJ_homf', 'Gnomad_exomes_FIN_hetf', 'Gnomad_exomes_FIN_homf',
               'Gnomad_exomes_AFR_hetf', 'Gnomad_exomes_AFR_homf', 'Gnomad_exomes_OTH_hetf',
               'Gnomad_exomes_OTH_homf', 'Gnomad_exomes_SAS_hetf', 'Gnomad_exomes_SAS_homf', 'MutationAssessor_score',
               'PROVEAN_score', 'GERP++_RS', 'integrated_fitCons_score',
               'phyloP100way_vertebrate', 'phyloP20way_mammalian', 'phastCons100way_vertebrate',
               'phastCons20way_mammalian', 'SiPhy_29way_logOdds', 'GenoCanyon_score', 'tags']
gnomad_list1 = ['Gnomad_exomes_AF', 'Gnomad_exomes_hetf', 'Gnomad_exomes_homf',
                'Gnomad_exomes_EAS_AF', 'Gnomad_exomes_EAS_hetf', 'Gnomad_exomes_EAS_homf',
                'Gnomad_exomes_NFE_AF', 'Gnomad_exomes_NFE_hetf', 'Gnomad_exomes_NFE_homf',
                'Gnomad_exomes_AMR_AF', 'Gnomad_exomes_AMR_hetf', 'Gnomad_exomes_AMR_homf',
                'Gnomad_exomes_ASJ_AF', 'Gnomad_exomes_ASJ_hetf', 'Gnomad_exomes_ASJ_homf',
                'Gnomad_exomes_FIN_AF', 'Gnomad_exomes_FIN_hetf', 'Gnomad_exomes_FIN_homf',
                'Gnomad_exomes_AFR_AF', 'Gnomad_exomes_AFR_hetf', 'Gnomad_exomes_AFR_homf',
                'Gnomad_exomes_OTH_AF', 'Gnomad_exomes_OTH_hetf', 'Gnomad_exomes_OTH_homf',
                'Gnomad_exomes_SAS_AF', 'Gnomad_exomes_SAS_hetf', 'Gnomad_exomes_SAS_homf']
gnomad_list2 = ['Gnomad_exomes_wtf', 'Gnomad_exomes_EAS_wtf', 'Gnomad_exomes_NFE_wtf',
                'Gnomad_exomes_AMR_wtf', 'Gnomad_exomes_ASJ_wtf', 'Gnomad_exomes_FIN_wtf',
                'Gnomad_exomes_AFR_wtf', 'Gnomad_exomes_OTH_wtf', 'Gnomad_exomes_SAS_wtf']
test_list = ['SIFT_score', 'cadd_snp_phred', 'PrimateAI', 'ReVe', 'mcap_v14', 'REVEL_score', 'ClinPred_Score',
             'MetaSVM_score', 'MetaLR_score', 'Polyphen2_HDIV_score', 'Polyphen2_HVAR_score', 'VEST3_score']


def feature_edit(data, gnomad_list1, gnomad_list2, features):
    # deal with region_based features
    for k, v in data.items():
        if k in ['CCRs', 'GeVIR_per', 'LOEUF_per', 'HIP_score', 'VIRLoF_per']:
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
            data[k].replace('', np.nan, inplace=True)
    # gnomad frequency
    for i in gnomad_list1:
        data[i].replace('.', 0, inplace=True)
    for i in gnomad_list2:
        data[i].replace('.', 1, inplace=True)
    # domain
    my_col = list(data['Interpro_domain'])
    for index, i in enumerate(my_col):
        if i == '.':
            my_col[index] = 0
        else:
            my_col[index] = 1
    data['Interpro_domain'] = my_col
    # total
    for i in features:
        data[i].replace('.', np.nan, inplace=True)
    data = data.apply(pd.to_numeric, errors='ignore')
    # print information
    print(data[features].shape)
    print(data[features].isnull().sum() / data[features].shape[0])
    return data[features]


def get_training1(var3, var1, var2, gnomad_list1, gnomad_list2, features):
    df_vs3 = pd.read_table(var3, low_memory=False)
    df_vs1 = pd.read_table(var1, low_memory=False)
    df_vs2 = pd.read_table(var2, low_memory=False)
    df_vs1['tags'] = 0
    df_vs2['tags'] = 1
    # vs3 tags
    my_column = list(df_vs3['ClinicalSignificance'])
    for index, i in enumerate(my_column):
        # Benign, Likely benign, and Benign/Likely benign
        if i.find('enign') != -1:
            my_column[index] = 'Benign'
        # Pathogenic, Likely pathogenic, and Pathogenic/Likely pathogenic
        # but not Conflicting interpretations of pathogenicity
        elif i.find('athogenic') != -1 and i.find('Conflicting interpretations of pathogenicity') == -1:
            my_column[index] = 'Pathogenic'
    df_vs3['tags'] = my_column
    df_vs3 = df_vs3[df_vs3['tags'].isin(['Benign', 'Pathogenic'])]
    df_vs3['tags'].replace({'Benign': 0, 'Pathogenic': 1}, inplace=True)
    data = pd.concat([df_vs2, df_vs1, df_vs3], sort=False)
    data.index = pd.RangeIndex(len(data.index))

    return feature_edit(data, gnomad_list1, gnomad_list2, features)


def get_training2(var2, var1, gnomad_list1, gnomad_list2, features):
    # data
    df_vs2 = pd.read_table(var2, low_memory=False)
    df_vs1 = pd.read_table(var1, low_memory=False)
    # tags
    df_vs1['tags'] = 0
    df_vs2['tags'] = 1
    # concat
    data = pd.concat([df_vs1, df_vs2], sort=False)
    data.index = pd.RangeIndex(len(data.index))

    return feature_edit(data, gnomad_list1, gnomad_list2, features)


def get_training3(var2, var3, gnomad_list1, gnomad_list2, features):
    df_vs2 = pd.read_table(var2, low_memory=False)
    df_vs3 = pd.read_table(var3, low_memory=False)
    # tags
    df_vs2['tags'] = 1
    my_column = list(df_vs3['ClinicalSignificance'])
    for index, i in enumerate(my_column):
        if i.find('enign') != -1:
            my_column[index] = 'Benign'
        elif i.find('athogenic') != -1 and i.find('Conflicting interpretations of pathogenicity') == -1:
            my_column[index] = 'Pathogenic'
    df_vs3['tags'] = my_column
    df_vs3 = df_vs3[df_vs3['tags'].isin(['Benign', 'Pathogenic'])]
    df_vs3['tags'].replace({'Benign': 0, 'Pathogenic': 1}, inplace=True)

    data = pd.concat([df_vs2, df_vs3], sort=False)
    data.index = pd.RangeIndex(len(data.index))

    return feature_edit(data, gnomad_list1, gnomad_list2, features)


def get_testset(data, gnomad_list1, gnomad_list2, features):
    df_test = pd.read_table(data)
    return feature_edit(df_test, gnomad_list1, gnomad_list2, features)


train1 = get_training1('../example/var3_demo.anno.txt', '../example/var1_demo.hg19_multianno.txt',
                       '../example/var2_demo.hg19_multianno.txt', gnomad_list1, gnomad_list2, categoryABC)
train2 = get_training2('../example/var2_demo.hg19_multianno.txt', '../example/var1_demo.hg19_multianno.txt',
                       gnomad_list1, gnomad_list2, categoryABC)
train3 = get_training3('../example/var2_demo.hg19_multianno.txt', '../example/var3_demo.anno.txt', gnomad_list1,
                       gnomad_list2, categoryABC)
test_set = get_testset('../example/ClinHGMD_demo.anno.txt', gnomad_list1, gnomad_list2, categoryABC + test_list)

# All the features were selected to provide complementary information
train1s = train1.dropna()
train2s = train2.dropna()
train3s = train3.dropna()
test_sets = test_set.dropna()
train1s.index = pd.RangeIndex(len(train1s.index))
train2s.index = pd.RangeIndex(len(train2s.index))
train3s.index = pd.RangeIndex(len(train3s.index))
test_sets.index = pd.RangeIndex(len(test_sets.index))
# print information
print('training set 1')
print(train1s['tags'].value_counts())
print('training set 2')
print(train2s['tags'].value_counts())
print('training set 3')
print(train3s['tags'].value_counts())
print('ClinHGMD')
print(test_sets['tags'].value_counts())

train1s.to_csv('../example/train1s_demo.txt', sep='\t', index=False)
train2s.to_csv('../example/train2s_demo.txt', sep='\t', index=False)
train3s.to_csv('../example/train3s_demo.txt', sep='\t', index=False)
test_sets.to_csv('../example/test_demo.txt', sep='\t', index=False)
