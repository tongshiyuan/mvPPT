import scipy.stats as ss
import numpy as np
from sklearn.preprocessing import Binarizer
from sklearn.metrics import accuracy_score, classification_report
from sklearn.metrics import precision_score, recall_score, f1_score, matthews_corrcoef, multilabel_confusion_matrix


def NRIcalculate(gold, prd_o, prd_n, cut_o, cut_n):
    '''
    :param gold: true tags of sample, type:list
    :param prd_o: old model predict value, type:list
    :param prd_n: new model predict value, type:list
    :param cut_o: threshold of old model, type:float
    :param cut_n: threshold of new model, type:float
    :return: no return, but print nri value, z value, and p value
    '''
    for i in gold:
        if i not in [0, 1]:
            print('only 0/1 can be identified')
            return
        elif not (len(gold) == len(prd_o) == len(prd_n)):
            print('wrong length of data')
            return
    B1 = 0
    C1 = 0
    B2 = 0
    C2 = 0
    N1 = 0
    N2 = 0
    for g, o, n in zip(gold, prd_o, prd_n):
        o_t = 1 if (o - cut_o > 0) else 0
        n_t = 1 if (n - cut_n > 0) else 0
        if g == 1:
            N1 += 1
            if o_t == 0 and n_t == 1:
                B1 += 1
            elif o_t == 1 and n_t == 0:
                C1 += 1
        else:
            N2 += 1
            if o_t == 0 and n_t == 1:
                B2 += 1
            elif o_t == 1 and n_t == 0:
                C2 += 1
    Additive_NRI = (B1 - C1) / N1 + (C2 - B2) / N2
    Z = Additive_NRI / np.sqrt((B1 + C1) / N1 ** 2 + (B2 + C2) / N2 ** 2)
    p = ss.norm.sf(abs(Z)) * 2

    print('Additive NRI:', Additive_NRI)
    print('Z:', Z)
    print('p:', p)


def metrics_of_model(data, model, features, th=0.5):
    """
    :param data: test set
    :param model: model trained by our study
    :param features: feature list
    :param th: threshold in out study is 0.5, and other tools use the threshold provide in their study
    :return:
    """
    y_test = data['tags']
    y_pred = Binarizer(threshold=th).transform(np.array(
        model.predict_proba(data[features].iloc[:, :len(features) - 1])[:, 1]).reshape(-1, 1))

    print(accuracy_score(y_test, y_pred))
    print(precision_score(y_test, y_pred))
    print(recall_score(y_test, y_pred))
    print(f1_score(y_test, y_pred))
    print(matthews_corrcoef(y_test, y_pred))
    mcm = multilabel_confusion_matrix(y_test, y_pred)
    tn = mcm[:, 0, 0]
    tp = mcm[:, 1, 1]
    fn = mcm[:, 1, 0]
    fp = mcm[:, 0, 1]
    TNR = tn / (tn + fp)
    FPR = fp / (tn + fp)
    print(TNR[-1])
    print(FPR[-1])
    print(classification_report(y_test, y_pred, digits=3))
