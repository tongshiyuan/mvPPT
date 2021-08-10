import pandas as pd

dfClinVar = pd.read_table('variant_summary_2020-07.txt.gz', compression='gzip', low_memory=False)
print('ClinVar 总共记录', dfClinVar.shape[0], '行数据（GRCh38 + GRCh37）')
dfClinGRCh37Snv = dfClinVar[
    (dfClinVar['Type'].isin(['single nucleotide variant'])) & (dfClinVar['Assembly'] == 'GRCh37')]
print('其中 GRCh37 的 SNV 有', dfClinGRCh37Snv.shape[0], '条')
print('-' * 60)
print('各标签数据统计：')
print(dfClinGRCh37Snv.loc[:, 'ClinicalSignificance'].value_counts())
dfClinGRCh37Snv[['Chromosome', 'Start', 'Stop', 'ReferenceAllele', 'AlternateAllele', 'GeneSymbol',
                 'ReviewStatus', 'ClinicalSignificance', 'LastEvaluated']].to_csv(
    '1_getSNVs/clinvarGRCh37snv.txt', sep='\t', index=False)

dfHgmd = pd.read_table('hgmd_hg19_vcf.txt', low_memory=False)
print('Hgmd vcf 总共记录', dfHgmd.shape[0], '行数据（hg19）')
dfAllMut = pd.read_table('hgmd_allmut.txt', low_memory=False)
dfHgmdDate = dfAllMut[['acc_num', 'new_date']].copy()
dfHgmdDate.rename(columns={'acc_num': 'id'}, inplace=True)
dfHgmdMerge = pd.merge(dfHgmd, dfHgmdDate, how='left', on='id')
dfHgmdSnv = dfHgmdMerge[(dfHgmdMerge.ref.str.len() == 1) & (dfHgmdMerge.alt.str.len() == 1)]
print('其中 SNV 有', dfHgmdSnv.shape[0], '条')
dfHgmdSnv[['chrom', 'pos', 'pos', 'ref', 'alt', 'info', 'new_date']].to_csv('1_getSNVs/hgmdHg19snv.txt', sep='\t',
                                                                            index=False)
