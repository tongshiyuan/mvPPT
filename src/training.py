import pickle
import warnings
import pandas as pd

import lightgbm as lgb  # 2.3.1

from bayes_opt import BayesianOptimization  # 1.1.0

warnings.filterwarnings("ignore")


def getTraning(file, features):
    trainingData = file[features]
    X = trainingData.drop('tags', axis=1)
    y = trainingData['tags']
    return X, y


def trainingModel(train, features, out, onehotFeature=[]):
    X, y = getTraning(train, features)
    opt_params = bayes_parameter_opt_lgb(X, y, onehotFeature, init_round=15, opt_round=100, n_folds=5, random_seed=1,
                                         n_estimators=100, learning_rate=0.05)
    model = training(X, y, opt_params)
    pickle.dump(model, open(out, 'wb'))


def bayes_parameter_opt_lgb(X, y, categorical_feats=[], init_round=15, opt_round=500, n_folds=5, random_seed=1,
                            n_estimators=10000, learning_rate=0.01, output_process=False):
    # prepare data
    train_data = lgb.Dataset(data=X, label=y, free_raw_data=False, categorical_feature=categorical_feats)

    # parameters
    def lgb_eval(num_leaves, feature_fraction, bagging_fraction, max_depth, lambda_l1, lambda_l2, min_split_gain,
                 min_child_weight):
        '''Step 1: parameters to be tuned,
        Note: values for parameters should make sense,
        e.g.: 'num_leaves' needs to be a integer and 'feature_fraction' should between 0 and 1'''
        params = {'application': 'binary',
                  'num_iterations': n_estimators,
                  'learning_rate': learning_rate,
                  'early_stopping_round': 100,
                  'metric': 'auc',
                  #                   'nthreads':16,
                  'boost_from_average': False}
        params["num_leaves"] = int(round(num_leaves))
        params['feature_fraction'] = max(min(feature_fraction, 1), 0)
        params['bagging_fraction'] = max(min(bagging_fraction, 1), 0)
        params['max_depth'] = int(round(max_depth))
        params['lambda_l1'] = max(lambda_l1, 0)
        params['lambda_l2'] = max(lambda_l2, 0)
        params['min_split_gain'] = min_split_gain
        params['min_child_weight'] = min_child_weight
        cv_result = lgb.cv(params, train_data, nfold=n_folds, seed=random_seed, stratified=True, verbose_eval=200,
                           metrics=['auc'])
        return max(cv_result['auc-mean'])

    # range
    # Step 2: Set the range for each parameter
    # Gentle reminder: try to make the range as narrow as possible
    lgbBO = BayesianOptimization(lgb_eval, {'num_leaves': (24, 45),
                                            'feature_fraction': (0.1, 0.9),
                                            'bagging_fraction': (0.8, 1),
                                            'max_depth': (5, 8.99),
                                            'lambda_l1': (0, 5),
                                            'lambda_l2': (0, 3),
                                            'min_split_gain': (0.001, 0.1),
                                            'min_child_weight': (5, 50)}, random_state=1)
    # optimize
    # Step 3: Bayesian Optimization: Maximize
    lgbBO.maximize(init_points=init_round, n_iter=opt_round)

    # output optimization process
    if output_process == True: lgbBO.points_to_csv("bayes_opt_result.csv")

    # return best parameters
    # Step 4: Get the parameters
    return lgbBO


def training(X, y, opt_params, n_estimators=100, learning_rate=0.01):
    param = {
        'num_leaves': int(round(opt_params.max['params']['num_leaves'])),
        'feature_fraction': opt_params.max['params']['feature_fraction'],
        'bagging_fraction': opt_params.max['params']['bagging_fraction'],
        'max_depth': int(round(opt_params.max['params']['max_depth'])),
        'lambda_l1': opt_params.max['params']['lambda_l1'],
        'lambda_l2': opt_params.max['params']['lambda_l2'],
        'min_split_gain': opt_params.max['params']['min_split_gain'],
        'min_child_weight': opt_params.max['params']['min_child_weight'],
        'application': 'binary',
        'num_iterations': n_estimators,
        'learning_rate': learning_rate,
        'metric': 'auc',
        #         'nthreads':16,
        'boost_from_average': False
    }

    model = lgb.LGBMClassifier(
        **param
    )
    model.fit(X, y)
    return model


train = pd.read_table('', low_memory=False)
trainFneatures = []
trainingModel(train, trainFneatures, 'model.m')
