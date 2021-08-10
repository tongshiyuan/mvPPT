import numpy as np
import pandas as pd
import random as rd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import Binarizer
from sklearn.metrics import accuracy_score,precision_recall_curve, average_precision_score, precision_score, f1_score, matthews_corrcoef, multilabel_confusion_matrix, log_loss, roc_curve, auc, recall_score,classification_report


def predict(model, data, features):
    model_pred = model.predict_proba(data[features[:-1]])
    return model_pred[:,1]

def pltRoc(y, pred_y, tag, n): 
    fpr, tpr, thresholds = roc_curve(y, pred_y, pos_label=1, drop_intermediate=False)
    plt.plot(fpr, tpr, color=n, lw=2, label=tag+'(%0.3f)' % auc(fpr, tpr))

def pltPrc(y, pred_y, tag, n): 
    precision, recall, _ = precision_recall_curve(y, pred_y)
    average_precision = average_precision_score(y, pred_y)
    plt.step(recall, precision,lw=2, color=(n),
         label=(tag + '({0:0.3f})'.format(average_precision)), where='post')

def Roc(dfTest, compareList,cl, file):
    plt.figure(figsize=(10,10.5))
    for index,i in enumerate(compareList):
        pltRoc(dfTest['tags'], dfTest[i],i, cl[index])
    plt.plot([0, 1], [0, 1], color='silver', lw=2, linestyle='--')
    plt.xlim([0.0, 1.00])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False positive rate', size = 16)
    plt.ylabel('True positive rate', size = 16)
    plt.legend(loc="lower right", fontsize=16)
    plt.xticks(size = 16)
    plt.yticks(size = 16)
    plt.locator_params()
    plt.savefig(file)

def Prc(dfTest, compareList,cl, file):
    plt.figure(figsize=(10,10.5))
    for index,i in enumerate(compareList):
        pltPrc(dfTest['tags'], dfTest[i],i, cl[index])
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall', size=16)
    plt.ylabel('Precision', size=16)
    plt.legend(loc="lower left", fontsize=16)
    plt.xticks(size = 16)
    plt.yticks(size = 16)
    plt.locator_params()
    plt.savefig(file)

def auprc_auroc(y, pred_y):
    fpr, tpr, thresholds = roc_curve(y, pred_y, pos_label=1, drop_intermediate=False)
    return (auc(fpr, tpr), average_precision_score(y, pred_y))
    

def au(dfTest, compareList):
    audict = dict.fromkeys(compareList,[])
    for index,i in enumerate(compareList):
        rst = auprc_auroc(dfTest['tags'], dfTest[i])
        audict[i]=rst
    return audict

def lower_sample_data(df, num=[500,4500]):
    TP_data = df[df['tags'] == 1]
    TN_data = df[df['tags'] == 0]
    indexTP = np.random.randint(len(TP_data), size=num[0])
    indexTN = np.random.randint(len(TN_data), size=num[1])
    #下采样后数据样本
    TP_data = TP_data.iloc[list(indexTP)]  # 下采样
    TN_data = TN_data.iloc[list(indexTN)]  # 下采样
    return pd.concat([TP_data, TN_data])
