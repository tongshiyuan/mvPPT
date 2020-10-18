import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score

# model was built in our study
model1 = pickle.load(open("../model/model1.m", 'rb'))
model2 = pickle.load(open("../model/model2.m", 'rb'))
model3 = pickle.load(open("../model/modelABC.m", 'rb'))
modelA = pickle.load(open("../model/modelA.m", 'rb'))
modelB = pickle.load(open("../model/modelB.m", 'rb'))
modelC = pickle.load(open("../model/modelC.m", 'rb'))
modelAB = pickle.load(open("../model/modelAB.m", 'rb'))
modelAC = pickle.load(open("../model/modelAC.m", 'rb'))
modelBC = pickle.load(open("../model/modelBC.m", 'rb'))
test_set = pd.read_table('../example/test_demo.txt')

GeneRegionScore = ['GeVIR_per', 'LOEUF_per', 'VIRLoF_per', 'HIP_score', 'CCRs', 'Interpro_domain']
GFsAFs = ['Gnomad_exomes_wtf', 'Gnomad_exomes_EAS_wtf', 'Gnomad_exomes_NFE_wtf',
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
          'Gnomad_exomes_OTH_homf', 'Gnomad_exomes_SAS_hetf', 'Gnomad_exomes_SAS_homf']
VariantPredictionScore = ['MutationAssessor_score', 'PROVEAN_score', 'GERP++_RS',
                          'integrated_fitCons_score', 'phyloP100way_vertebrate',
                          'phyloP20way_mammalian', 'phastCons100way_vertebrate',
                          'phastCons20way_mammalian', 'SiPhy_29way_logOdds', 'GenoCanyon_score']
categoryA = GFsAFs + ['tags']
categoryB = GeneRegionScore + ['tags']
categoryC = VariantPredictionScore + ['tags']
categoryAB = GeneRegionScore + GFsAFs + ['tags']
categoryAC = GFsAFs + VariantPredictionScore + ['tags']
categoryBC = GeneRegionScore + VariantPredictionScore + ['tags']
categoryABC = GeneRegionScore + GFsAFs + VariantPredictionScore + ['tags']


def plot_roc(dataset):
    np.random.seed(2020)
    X_test = dataset.iloc[:, :-13]
    y_test = dataset['tags']

    y_pred1s_proba = model1.predict_proba(X_test)
    y_pred2s_proba = model2.predict_proba(X_test)
    y_pred3s_proba = model3.predict_proba(X_test)

    fpr1s, tpr1s, threshold1s = roc_curve(y_test,
                                          y_pred1s_proba[:, 1], pos_label=1, drop_intermediate=False)
    roc_auc1s = auc(fpr1s, tpr1s)
    fpr2s, tpr2s, threshold2s = roc_curve(y_test,
                                          y_pred2s_proba[:, 1], pos_label=1, drop_intermediate=False)
    roc_auc2s = auc(fpr2s, tpr2s)
    fpr3s, tpr3s, threshold3s = roc_curve(y_test,
                                          y_pred3s_proba[:, 1], pos_label=1, drop_intermediate=False)
    roc_auc3s = auc(fpr3s, tpr3s)

    # revel
    fpr_r, tpr_r, threshold_r = roc_curve(y_test, dataset['REVEL_score'], pos_label=1)
    roc_auc_r = auc(fpr_r, tpr_r)
    # mcap
    fpr_m, tpr_m, threshold_m = roc_curve(y_test, dataset['mcap_v14'], pos_label=1)
    roc_auc_m = auc(fpr_m, tpr_m)
    # clinpred
    fpr_c, tpr_c, threshold_c = roc_curve(y_test, dataset['ClinPred_Score'], pos_label=1)
    roc_auc_c = auc(fpr_c, tpr_c)
    # ReVe
    fpr_reve, tpr_reve, threshold_reve = roc_curve(y_test, dataset['ReVe'], pos_label=1)
    roc_auc_reve = auc(fpr_reve, tpr_reve)
    # PrimateAI
    fpr_pri, tpr_pri, threshold_pri = roc_curve(y_test, dataset['PrimateAI'], pos_label=1)
    roc_auc_pri = auc(fpr_pri, tpr_pri)
    # MetaSVM
    fpr_MetaSVM, tpr_MetaSVM, threshold_MetaSVM = roc_curve(y_test, dataset['MetaSVM_score'], pos_label=1)
    roc_auc_MetaSVM = auc(fpr_MetaSVM, tpr_MetaSVM)
    # MetaLR
    fpr_MetaLR, tpr_MetaLR, threshold_MetaLR = roc_curve(y_test, dataset['MetaLR_score'], pos_label=1)
    roc_auc_MetaLR = auc(fpr_MetaLR, tpr_MetaLR)
    # SIFT
    fpr_SIFT, tpr_SIFT, threshold_SIFT = roc_curve(y_test, 1 - dataset['SIFT_score'], pos_label=1)
    roc_auc_SIFT = auc(fpr_SIFT, tpr_SIFT)
    # CADD
    fpr_cadd, tpr_cadd, threshold_cadd = roc_curve(y_test, dataset['cadd_snp_phred'], pos_label=1)
    roc_auc_cadd = auc(fpr_cadd, tpr_cadd)
    # HDIV
    fpr_hdiv, tpr_hdiv, threshold_hdiv = roc_curve(y_test, dataset['Polyphen2_HDIV_score'], pos_label=1)
    roc_auc_hdiv = auc(fpr_hdiv, tpr_hdiv)
    # HVAR
    fpr_hvar, tpr_hvar, threshold_hvar = roc_curve(y_test, dataset['Polyphen2_HVAR_score'], pos_label=1)
    roc_auc_hvar = auc(fpr_hvar, tpr_hvar)
    # VEST3
    fpr_vest, tpr_vest, threshold_vest = roc_curve(y_test, dataset['VEST3_score'], pos_label=1)
    roc_auc_vest = auc(fpr_vest, tpr_vest)

    lw = 2
    plt.figure(figsize=(10, 10.5))

    plt.plot(fpr1s, tpr1s, color='#459B77',
             lw=lw, label='Model_1(%0.3f)' % roc_auc1s)
    plt.plot(fpr2s, tpr2s, color='#3B76B0',
             lw=lw, label='Model_2(%0.3f)' % roc_auc2s)
    plt.plot(fpr3s, tpr3s, color='r',
             lw=lw, label='Model_3(%0.3f)' % roc_auc3s)

    plt.plot(fpr_cadd, tpr_cadd, color='#CC3399',
             lw=lw, label='CADD(%0.3f)' % roc_auc_cadd)
    plt.plot(fpr_c, tpr_c, color='#CC99CC',
             lw=lw, label='ClinPred(%0.3f)' % roc_auc_c)
    plt.plot(fpr_m, tpr_m, color='#9933CC',
             lw=lw, label='M-CAP(%0.3f)' % roc_auc_m)
    plt.plot(fpr_hdiv, tpr_hdiv, color='#006600',
             lw=lw, label='Polyphen-2 HDIV(%0.3f)' % roc_auc_hdiv)
    plt.plot(fpr_hvar, tpr_hvar, color='#009933',
             lw=lw, label='Polyphen-2 HVAR(%0.3f)' % roc_auc_hvar)
    plt.plot(fpr_reve, tpr_reve, color='#66CC66',
             lw=lw, label='ReVe(%0.3f)' % roc_auc_reve)
    plt.plot(fpr_SIFT, tpr_SIFT, color='#0099CC',
             lw=lw, label='SIFT(%0.3f)' % roc_auc_SIFT)
    plt.plot(fpr_vest, tpr_vest, color='#0099FF',
             lw=lw, label='VEST3(%0.3f)' % roc_auc_vest)
    plt.plot(fpr_pri, tpr_pri, color='#66CCFF',
             lw=lw, label='PrimateAI(%0.3f)' % roc_auc_pri)
    plt.plot(fpr_r, tpr_r, color='#CC6600',
             lw=lw, label='REVEL(%0.3f)' % roc_auc_r)
    plt.plot(fpr_MetaLR, tpr_MetaLR, color='#FF9933',
             lw=lw, label='MetaLR(%0.3f)' % roc_auc_MetaLR)
    plt.plot(fpr_MetaSVM, tpr_MetaSVM, color='#FFCC99',
             lw=lw, label='MetaSVM(%0.3f)' % roc_auc_MetaSVM)

    plt.plot([0, 1], [0, 1], color='silver', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.00])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False positive rate', size=16)
    plt.ylabel('True positive rate', size=16)
    plt.legend(loc="lower right")
    plt.xticks(fontproperties='Times New Roman', size=16)
    plt.yticks(fontproperties='Times New Roman', size=16)
    plt.locator_params()
    plt.show()


def plot_prc(dataset):
    np.random.seed(2020)
    X_test = dataset.iloc[:, :-13]
    y_test = dataset['tags']

    y_pred1s_proba = model1.predict_proba(X_test)
    y_pred2s_proba = model2.predict_proba(X_test)
    y_pred3s_proba = model3.predict_proba(X_test)

    precision1s, recall1s, _1s = precision_recall_curve(y_test, y_pred1s_proba[:, 1])
    average_precision1s = average_precision_score(y_test, y_pred1s_proba[:, 1])
    precision2s, recall2s, _2s = precision_recall_curve(y_test, y_pred2s_proba[:, 1])
    average_precision2s = average_precision_score(y_test, y_pred2s_proba[:, 1])
    precision3s, recall3s, _3s = precision_recall_curve(y_test, y_pred3s_proba[:, 1])
    average_precision3s = average_precision_score(y_test, y_pred3s_proba[:, 1])

    # revel
    precision_r, recall_r, _r = precision_recall_curve(y_test, dataset['REVEL_score'])
    average_precision_r = average_precision_score(y_test, dataset['REVEL_score'])
    # mcap
    precision_m, recall_m, _m = precision_recall_curve(y_test, dataset['mcap_v14'])
    average_precision_m = average_precision_score(y_test, dataset['mcap_v14'])
    # clinpred
    precision_c, recall_c, _c = precision_recall_curve(y_test, dataset['ClinPred_Score'])
    average_precision_c = average_precision_score(y_test, dataset['ClinPred_Score'])
    # ReVe
    precision_reve, recall_reve, _reve = precision_recall_curve(y_test, dataset['ReVe'])
    average_precision_reve = average_precision_score(y_test, dataset['ReVe'])
    # PrimateAI
    precision_pri, recall_pri, _pri = precision_recall_curve(y_test, dataset['PrimateAI'])
    average_precision_pri = average_precision_score(y_test, dataset['PrimateAI'])
    # MetaSVM
    precision_MetaSVM, recall_MetaSVM, _MetaSVM = precision_recall_curve(y_test, dataset['MetaSVM_score'])
    average_precision_MetaSVM = average_precision_score(y_test, dataset['MetaSVM_score'])
    # MetaLR
    precision_MetaLR, recall_MetaLR, _MetaLR = precision_recall_curve(y_test, dataset['MetaLR_score'])
    average_precision_MetaLR = average_precision_score(y_test, dataset['MetaLR_score'])
    # SIFT
    precision_SIFT, recall_SIFT, _SIFT = precision_recall_curve(y_test, 1 - dataset['SIFT_score'])
    average_precision_SIFT = average_precision_score(y_test, 1 - dataset['SIFT_score'])
    # CADD
    precision_cadd, recall_cadd, _cadd = precision_recall_curve(y_test, dataset['cadd_snp_phred'])
    average_precision_cadd = average_precision_score(y_test, dataset['cadd_snp_phred'])
    # HDIV
    precision_hdiv, recall_hdiv, _hdiv = precision_recall_curve(y_test, dataset['Polyphen2_HDIV_score'])
    average_precision_hdiv = average_precision_score(y_test, dataset['Polyphen2_HDIV_score'])
    # HVAR
    precision_hvar, recall_hvar, _hvar = precision_recall_curve(y_test, dataset['Polyphen2_HVAR_score'])
    average_precision_hvar = average_precision_score(y_test, dataset['Polyphen2_HVAR_score'])
    # VEST3
    precision_vest, recall_vest, _vest = precision_recall_curve(y_test, dataset['VEST3_score'])
    average_precision_vest = average_precision_score(y_test, dataset['VEST3_score'])

    plt.figure(figsize=(10, 10.5))
    lw = 2

    plt.step(recall1s, precision1s, lw=lw, color='#459B77',
             label='Model_1({0:0.3f})'.format(average_precision1s), where='post')
    plt.step(recall2s, precision2s, lw=lw, color='#3B76B0',
             label='Model_2({0:0.3f})'.format(average_precision2s), where='post')
    plt.step(recall3s, precision3s, lw=lw, color='r',
             label='Model_3({0:0.3f})'.format(average_precision3s), where='post')

    plt.step(recall_cadd, precision_cadd, lw=lw, color='#CC3399',
             label='CADD({0:0.3f})'.format(average_precision_cadd), where='post')
    plt.step(recall_c, precision_c, lw=lw, color='#CC99CC',
             label='ClinPred({0:0.3f})'.format(average_precision_c), where='post')
    plt.step(recall_m, precision_m, lw=lw, color='#9933CC',
             label='M-CAP({0:0.3f})'.format(average_precision_m), where='post')
    plt.step(recall_hdiv, precision_hdiv, lw=lw, color='#006600',
             label='Polyphen-2 HDIV({0:0.3f})'.format(average_precision_hdiv), where='post')
    plt.step(recall_hvar, precision_hvar, lw=lw, color='#009933',
             label='Polyphen-2 HVAR({0:0.3f})'.format(average_precision_hvar), where='post')
    plt.step(recall_reve, precision_reve, lw=lw, color='#66CC66',
             label='ReVe({0:0.3f})'.format(average_precision_reve), where='post')
    plt.step(recall_SIFT, precision_SIFT, lw=lw, color='#0099CC',
             label='SIFT({0:0.3f})'.format(average_precision_SIFT), where='post')
    plt.step(recall_vest, precision_vest, lw=lw, color='#0099FF',
             label='VEST3({0:0.3f})'.format(average_precision_vest), where='post')
    plt.step(recall_pri, precision_pri, lw=lw, color='#66CCFF',
             label='PrimateAI({0:0.3f})'.format(average_precision_pri), where='post')
    plt.step(recall_r, precision_r, lw=lw, color='#CC6600',
             label='REVEL({0:0.3f})'.format(average_precision_r), where='post')
    plt.step(recall_MetaLR, precision_MetaLR, lw=lw, color='#FF9933',
             label='MetaLR({0:0.3f})'.format(average_precision_MetaLR), where='post')
    plt.step(recall_MetaSVM, precision_MetaSVM, lw=lw, color='#FFCC99',
             label='MetaSVM({0:0.3f})'.format(average_precision_MetaSVM), where='post')

    plt.xlabel('Recall', size=16)
    plt.ylabel('Precision', size=16)
    plt.ylim([0.0, 1.05])
    plt.xlim([0.0, 1.0])
    plt.legend(loc="lower left")
    plt.xticks(fontproperties='Times New Roman', size=16)
    plt.yticks(fontproperties='Times New Roman', size=16)
    plt.locator_params()
    plt.show()


def category_roc(data):
    lw = 2
    plt.figure(figsize=(10, 10.5))
    plt.plot([0, 1], [0, 1], color='silver', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False positive rate', size=16)
    plt.ylabel('True positive rate', size=16)
    plt.xticks(fontproperties='Times New Roman', size=16)
    plt.yticks(fontproperties='Times New Roman', size=16)
    plt.locator_params()
    np.random.seed(2020)

    # A
    test_set = data[categoryA]
    X_test = test_set
    y_test = test_set['tags']
    y_pred_proba = modelA.predict_proba(X_test.iloc[:, :-1])
    fpr, tpr, thrauprcold = roc_curve(y_test, y_pred_proba[:, 1], pos_label=1, drop_intermediate=False)
    roc_auc = auc(fpr, tpr)
    plt.plot(fpr, tpr, color='#CC3399',
             lw=lw, label='A (%0.3f)' % roc_auc)

    # B
    test_set = data[categoryB]
    X_test = test_set
    y_pred_proba = modelB.predict_proba(X_test.iloc[:, :-1])
    fpr, tpr, threshold = roc_curve(y_test, y_pred_proba[:, 1], pos_label=1, drop_intermediate=False)
    roc_auc = auc(fpr, tpr)

    plt.plot(fpr, tpr, color='#CC99CC',
             lw=lw, label='B(%0.3f)' % roc_auc)
    # C
    test_set = data[categoryC]
    X_test = test_set
    y_pred_proba = modelC.predict_proba(X_test.iloc[:, :-1])
    fpr, tpr, threshold = roc_curve(y_test, y_pred_proba[:, 1], pos_label=1, drop_intermediate=False)
    roc_auc = auc(fpr, tpr)

    plt.plot(fpr, tpr, color='#006600',
             lw=lw, label='C(%0.3f)' % roc_auc)

    # AB
    test_set = data[categoryAB]
    X_test = test_set
    y_pred_proba = modelAB.predict_proba(X_test.iloc[:, :-1])
    fpr, tpr, threshold = roc_curve(y_test, y_pred_proba[:, 1], pos_label=1, drop_intermediate=False)
    roc_auc = auc(fpr, tpr)
    plt.plot(fpr, tpr, color='#66CC66',
             lw=lw, label='AB(%0.3f)' % roc_auc)

    # AC
    test_set = data[categoryAC]
    X_test = test_set
    y_pred_proba = modelAC.predict_proba(X_test.iloc[:, :-1])
    fpr, tpr, threshold = roc_curve(y_test, y_pred_proba[:, 1], pos_label=1, drop_intermediate=False)
    roc_auc = auc(fpr, tpr)

    plt.plot(fpr, tpr, color='#0099CC',
             lw=lw, label='AC(%0.3f)' % roc_auc)
    # BC
    test_set = data[categoryBC]

    X_test = test_set
    y_pred_proba = modelBC.predict_proba(X_test.iloc[:, :-1])

    fpr, tpr, threshold = roc_curve(y_test, y_pred_proba[:, 1], pos_label=1, drop_intermediate=False)
    roc_auc = auc(fpr, tpr)

    plt.plot(fpr, tpr, color='#66CCFF',
             lw=lw, label='BC(%0.3f)' % roc_auc)

    # ABC
    test_set = data[categoryABC]
    X_test = test_set
    y_pred_proba = model3.predict_proba(X_test.iloc[:, :-1])

    fpr, tpr, threshold = roc_curve(y_test, y_pred_proba[:, 1], pos_label=1, drop_intermediate=False)
    roc_auc = auc(fpr, tpr)
    plt.plot(fpr, tpr, color='r',
             lw=lw, label='ABC(%0.3f)' % roc_auc)

    plt.legend(loc="lower right")
    plt.show()


def category_pr(data):
    np.random.seed(2020)
    lw = 2
    plt.figure(figsize=(10, 10.5))
    # plt.plot([0, 1], [0, 1], color='silver', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False positive rate', size=16)
    plt.ylabel('True positive rate', size=16)
    plt.xticks(fontproperties='Times New Roman', size=16)
    plt.yticks(fontproperties='Times New Roman', size=16)
    plt.locator_params()
    # A
    test_set = data[categoryA]
    X_test = test_set
    y_test = test_set['tags']
    y_pred_proba = modelA.predict_proba(X_test.iloc[:, :-1])

    precision, recall, _ = precision_recall_curve(y_test, y_pred_proba[:, 1])
    average_precision = average_precision_score(y_test, y_pred_proba[:, 1])
    plt.step(recall, precision, lw=lw, color='#CC3399',
             label='A({0:0.3f})'.format(average_precision), where='post')

    # B
    test_set = data[categoryB]
    X_test = test_set
    y_pred_proba = modelB.predict_proba(X_test.iloc[:, :-1])

    precision, recall, _ = precision_recall_curve(y_test, y_pred_proba[:, 1])
    average_precision = average_precision_score(y_test, y_pred_proba[:, 1])
    plt.step(recall, precision, lw=lw, color='#CC99CC',
             label='B({0:0.3f})'.format(average_precision), where='post')

    # C
    test_set = data[categoryC]

    X_test = test_set
    y_pred_proba = modelC.predict_proba(X_test.iloc[:, :-1])

    precision, recall, _ = precision_recall_curve(y_test, y_pred_proba[:, 1])
    average_precision = average_precision_score(y_test, y_pred_proba[:, 1])
    plt.step(recall, precision, lw=lw, color='#006600',
             label='C({0:0.3f})'.format(average_precision), where='post')

    # AB
    test_set = data[categoryAB]

    X_test = test_set
    y_pred_proba = modelAB.predict_proba(X_test.iloc[:, :-1])

    precision, recall, _ = precision_recall_curve(y_test, y_pred_proba[:, 1])
    average_precision = average_precision_score(y_test, y_pred_proba[:, 1])
    plt.step(recall, precision, lw=lw, color='#66CC66',
             label='AB({0:0.3f})'.format(average_precision), where='post')

    # AC
    test_set = data[categoryAC]

    X_test = test_set
    y_pred_proba = modelAC.predict_proba(X_test.iloc[:, :-1])

    precision, recall, _ = precision_recall_curve(y_test, y_pred_proba[:, 1])
    average_precision = average_precision_score(y_test, y_pred_proba[:, 1])
    plt.step(recall, precision, lw=lw, color='#0099CC',
             label='AC({0:0.3f})'.format(average_precision), where='post')

    # BC
    test_set = data[categoryBC]

    X_test = test_set
    y_pred_proba = modelBC.predict_proba(X_test.iloc[:, :-1])

    precision, recall, _ = precision_recall_curve(y_test, y_pred_proba[:, 1])
    average_precision = average_precision_score(y_test, y_pred_proba[:, 1])
    plt.step(recall, precision, lw=lw, color='#66CCFF',
             label='BC({0:0.3f})'.format(average_precision), where='post')

    # ABC
    test_set = data[categoryABC]

    X_test = test_set
    y_pred_proba = model3.predict_proba(X_test.iloc[:, :-1])

    precision, recall, _ = precision_recall_curve(y_test, y_pred_proba[:, 1])
    average_precision = average_precision_score(y_test, y_pred_proba[:, 1])
    plt.step(recall, precision, lw=lw, color='r',
             label='ABC({0:0.3f})'.format(average_precision), where='post')

    plt.legend(loc="lower left")
    plt.show()


def relative_importance(features):
    importance_df = pd.DataFrame({
        'Features': features,
        'Model': model3.feature_importances_,

    })
    importance_df['Model'] = importance_df['Model'] / np.sum(importance_df['Model'])
    importance_df['total'] = importance_df.apply(lambda x: x['Model'], axis=1)
    df_f = importance_df.sort_values(['total'], ascending=False)
    df_f.drop(['total'], axis=1, inplace=True)
    plt.figure(figsize=(12, 7))
    plt.bar(x=df_f['Features'], height=df_f['Model'], color='grey')
    plt.tick_params(axis='x', labelsize=10)
    tmp = plt.xticks(rotation=90)
    plt.locator_params()
    plt.ylabel('Relative importance')
    plt.xlabel('Features')
    plt.show()


def single_plot(test_set, features):
    y_test = test_set['tags']
    auc_list = []
    prc_list = []
    for feature in features[:-1]:
        model = pickle.load(open('../model/model_%s.m' % feature, 'rb'))
        y_pred_proba = model.predict_proba(test_set[feature].values.reshape(-1, 1))
        fpr, tpr, threshold = roc_curve(y_test, y_pred_proba[:, 1], pos_label=1, drop_intermediate=False)
        roc_auc = auc(fpr, tpr)
        precision, recall, _ = precision_recall_curve(y_test, y_pred_proba[:, 1])
        average_precision = average_precision_score(y_test, y_pred_proba[:, 1])
        auc_list.append(roc_auc)
        prc_list.append(average_precision)
    plt.figure(figsize=(6, 12))
    plt.xlim([0.45, 1.0])
    plt.xticks(fontproperties='Times New Roman', size=16)
    plt.yticks(fontproperties='Times New Roman')
    plt.locator_params()
    plt.plot(auc_list, features[:-1], 'o', color='k')

    plt.figure(figsize=(6, 12))
    plt.xlim([0.45, 1.0])
    plt.xticks(fontproperties='Times New Roman', size=16)
    plt.yticks(fontproperties='Times New Roman')
    plt.locator_params()
    plt.plot(prc_list, features[:-1], 'o', color='k')
    plt.show()


plot_roc(test_set)
plot_prc(test_set)
category_roc(test_set)
category_pr(test_set)
relative_importance(categoryABC[:-1])
single_plot(test_set, categoryABC)
