import pickle
import warnings
import numpy as np
import pandas as pd
import lightgbm as lgb
from sklearn.metrics import roc_auc_score
from bayes_opt import BayesianOptimization
from sklearn.model_selection import train_test_split

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


def train_model(train_data):
    np.random.seed(2020)
    X_train = train_data
    traind, val = train_test_split(X_train, test_size=0.2, random_state=2020)
    train_y = traind['tags']
    val_y = val['tags']
    train_x = traind.drop(['tags'], axis=1)
    val_x = val.drop(['tags'], axis=1)

    # BayesianOptimization
    def LGB_bayesian(
            num_leaves,  # int
            min_data_in_leaf,  # int
            learning_rate,
            min_sum_hessian_in_leaf,  # int
            feature_fraction,
            lambda_l1,
            lambda_l2,
            min_gain_to_split,
            max_depth,
            bagging_freq,
            bagging_fraction,
            max_bin):
        # LightGBM expects next three parameters need to be integer. So we make them integer
        num_leaves = int(num_leaves)
        min_data_in_leaf = int(min_data_in_leaf)
        max_depth = int(max_depth)
        bagging_freq = int(bagging_freq)
        max_bin = int(max_bin)

        assert type(num_leaves) == int
        assert type(min_data_in_leaf) == int
        assert type(max_depth) == int

        param = {
            'num_leaves': num_leaves,
            'min_data_in_leaf': min_data_in_leaf,
            'learning_rate': learning_rate,
            'min_sum_hessian_in_leaf': min_sum_hessian_in_leaf,
            'feature_fraction': feature_fraction,
            'lambda_l1': lambda_l1,
            'lambda_l2': lambda_l2,
            'min_gain_to_split': min_gain_to_split,
            'max_depth': max_depth,
            'save_binary': True,
            #         'max_bin': 63,
            #         'bagging_fraction': 0.4,
            #         'bagging_freq': 5,
            'max_bin': max_bin,
            'bagging_fraction': bagging_fraction,
            'bagging_freq': bagging_freq,
            'seed': 2020,
            # 'feature_fraction_seed': 2019,
            # 'bagging_seed': 2019,
            # 'drop_seed': 2019,
            # 'data_random_seed': 2019,
            'objective': 'binary',
            'boosting_type': 'gbdt',
            'verbose': -1,
            'metric': 'auc',
            # "tree_learner": "serial",
            # 'is_unbalance': True,
            # 'boost_from_average': False,
        }
        lgtrain = lgb.Dataset(train_x, label=train_y)
        lgval = lgb.Dataset(val_x, label=val_y)
        model = lgb.train(param, lgtrain, 20000, valid_sets=[lgval], early_stopping_rounds=100, verbose_eval=3000)
        pred_val_y = model.predict(val_x, num_iteration=model.best_iteration)
        score = roc_auc_score(val_y, pred_val_y)
        return score

    bounds_LGB = {
        'num_leaves': (5, 20),
        'min_data_in_leaf': (5, 100),
        'learning_rate': (0.005, 0.3),
        'min_sum_hessian_in_leaf': (0.00001, 20),
        'feature_fraction': (0.001, 0.5),
        'lambda_l1': (0, 10),
        'lambda_l2': (0, 10),
        'min_gain_to_split': (0, 1.0),
        'max_depth': (3, 200),
        'bagging_fraction': (0.5, 1),
        'max_bin': (5, 256),
        'bagging_freq': (0, 90),
    }

    LGB_BO = BayesianOptimization(LGB_bayesian, bounds_LGB, random_state=2020)

    init_points = 10
    n_iter = 200
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        LGB_BO.maximize(init_points=init_points, n_iter=n_iter, acq='ucb', xi=0.0, alpha=1e-6)

    print(LGB_BO.max['target'])
    print(LGB_BO.max['params'])

    param = {
        'num_leaves': int(round(LGB_BO.max['params']['num_leaves'])),  # remember to int here
        'max_bin': int(round(LGB_BO.max['params']['max_bin'])),
        'min_data_in_leaf': int(round(LGB_BO.max['params']['min_data_in_leaf'])),  # remember to int here
        'learning_rate': LGB_BO.max['params']['learning_rate'],
        'min_sum_hessian_in_leaf': LGB_BO.max['params']['min_sum_hessian_in_leaf'],
        'bagging_fraction': LGB_BO.max['params']['bagging_fraction'],
        'bagging_freq': int(round(LGB_BO.max['params']['bagging_freq'])),
        'feature_fraction': LGB_BO.max['params']['feature_fraction'],
        'lambda_l1': LGB_BO.max['params']['lambda_l1'],
        'lambda_l2': LGB_BO.max['params']['lambda_l2'],
        'min_gain_to_split': LGB_BO.max['params']['min_gain_to_split'],
        'max_depth': int(round(LGB_BO.max['params']['max_depth'])),  # remember to int here
        'objective': 'binary',
        'boosting_type': 'gbdt',
        'metric': 'auc',
    }

    model = lgb.LGBMClassifier(
        **param
    )
    model.fit(train_data.iloc[:, :-1], train_data.iloc[:, -1])
    return model, param


train1s = pd.read_table('../example/train1s.txt')
train2s = pd.read_table('../example/train2s.txt')
train3s = pd.read_table('../example/train3s.txt')

modelA, paramA = train_model(train3s[categoryA])
pickle.dump(modelA, open("../example/modelA_demo.m", 'wb'))
modelB, paramB = train_model(train3s[categoryB])
pickle.dump(modelB, open("../example/modelB_demo.m", 'wb'))
modelC, paramC = train_model(train3s[categoryC])
pickle.dump(modelC, open("../example/modelC_demo.m", 'wb'))
modelAB, paramAB = train_model(train3s[categoryAB])
pickle.dump(modelAB, open("../example/modelAB_demo.m", 'wb'))
modelBC, paramBC = train_model(train3s[categoryBC])
pickle.dump(modelBC, open("../example/modelBC_demo.m", 'wb'))
modelAC, paramAC = train_model(train3s[categoryAC])
pickle.dump(modelAC, open("../example/modelAC_demo.m", 'wb'))
modelABC, paramABC = train_model(train3s[categoryABC])
pickle.dump(modelABC, open("../example/modelABC_demo.m", 'wb'))

model1, param1 = train_model(train1s[categoryABC])
model2, param2 = train_model(train2s[categoryABC])
pickle.dump(model1, open("../example/model1_demo.m", 'wb'))
pickle.dump(model2, open("../example/model2_demo.m", 'wb'))

singleFeatureParam = {}
for i in categoryABC[:-1]:
    model, param = train_model(train3s[[i, 'tags']])
    pickle.dump(model, open("../example/model_%s.m" % i, 'wb'))
    singleFeatureParam[i] = param
