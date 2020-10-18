�clightgbm.sklearn
LGBMClassifier
q )�q}q(X   boosting_typeqX   gbdtqX	   objectiveqX   binaryqX
   num_leavesqKX	   max_depthqK�X   learning_rateq	cnumpy.core.multiarray
scalar
q
cnumpy
dtype
qX   f8qK K�qRq(KX   <qNNNJ����J����K tqbC�/Ϡ�<�?q�qRqX   n_estimatorsqKdX   subsample_for_binqJ@ X   min_split_gainqG        X   min_child_weightqG?PbM���X   min_child_samplesqKX	   subsampleqG?�      X   subsample_freqqK X   colsample_bytreeqG?�      X	   reg_alphaqG        X
   reg_lambdaqG        X   random_stateqNX   n_jobsqJ����X   silentq �X   importance_typeq!X   splitq"X   _Boosterq#clightgbm.basic
Booster
q$)�q%}q&(X   handleq'X�3  tree
version=v3
num_class=1
num_tree_per_iteration=1
label_index=0
max_feature_idx=0
objective=binary sigmoid:1
feature_names=Gnomad_exomes_ASJ_wtf
feature_infos=[0:1]
tree_sizes=533 443 442 549 551 551 552 447 551 448 450 451 449 450 452 343 342 343 345 344 343 346 347 347

Tree=0
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=7435.27 64.8158 45.5995
threshold=0.99980113353686251 0.99990079365079365 0.99940393401183547
decision_type=2 2 2
left_child=2 -2 -1
right_child=1 -3 -4
leaf_value=0.25608218423178197 1.1837685906029789 1.5047180470970976 0.9426137150543914
leaf_weight=335.72276160120964 29.031098708510399 2876.9981005042791 30.977317616343498
leaf_count=2070 179 17739 191
internal_value=0 0.538701 -4.19915
internal_weight=0 2906.03 366.7
internal_count=20179 17918 2261
shrinkage=1


Tree=1
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=2628.96 24.6833
threshold=0.99990079365079365 0.9982121572596041
decision_type=2 2
left_child=1 -1
right_child=-2 -3
leaf_value=-0.51739881049905323 0.11660259002900916 -0.27944998012215444
leaf_weight=461.14828437566757 2637.7837287336588 118.61670273542404
leaf_count=1875 17739 565
internal_value=0 -1.88013
internal_weight=0 579.765
internal_count=20179 2440
shrinkage=0.253694


Tree=2
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=1490.35 22.9296
threshold=0.99990079365079365 0.99254269610309087
decision_type=2 2
left_child=1 -1
right_child=-2 -3
leaf_value=-0.40780444389120807 0.095149078116396377 -0.23629031204876877
leaf_weight=371.12812414765358 2444.2619986832142 218.0327672958374
leaf_count=1510 17739 930
internal_value=0 -1.388
internal_weight=0 589.161
internal_count=20179 2440
shrinkage=0.253694


Tree=3
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=927.682 9.36436 3.4466
threshold=0.99980113353686251 0.99990079365079365 0.99254269610309087
decision_type=2 2 2
left_child=2 -2 -1
right_child=1 -3 -4
leaf_value=-0.33979181114914386 -0.0074010894873951863 0.07664785990098158 -0.19936230426926679
leaf_weight=338.21019694209099 40.106593117117882 2290.5147879570723 179.40309892594814
leaf_count=1510 179 17739 751
internal_value=0 0.292097 -1.17968
internal_weight=0 2330.62 517.613
internal_count=20179 17918 2261
shrinkage=0.253694


Tree=4
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=591.341 9.08314 6.50911
threshold=0.99980113353686251 0.9982121572596041 0.99990079365079365
decision_type=2 2 2
left_child=1 -1 -2
right_child=2 -3 -4
leaf_value=-0.27955438586317827 -0.0057916393846594765 -0.10823391432845045 0.060972576064311404
leaf_weight=376.66668854653835 40.201832011342049 93.407046690583229 2170.0158790871501
leaf_count=1875 179 386 17739
internal_value=0 -0.996929 0.231005
internal_weight=0 470.074 2210.22
internal_count=20179 2261 17918
shrinkage=0.253694


Tree=5
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=382.746 13.9129 4.42841
threshold=0.99980113353686251 0.99254269610309087 0.99990079365079365
decision_type=2 2 2
left_child=1 -1 -2
right_child=2 -3 -4
leaf_value=-0.26767742447091519 -0.0045341004223153599 -0.11733611533589891 0.047945908823872313
leaf_weight=255.80577775835991 40.275823131203651 165.1482729613781 2076.5629931017756
leaf_count=1510 179 751 17739
internal_value=0 -0.857377 0.180322
internal_weight=0 420.954 2116.84
internal_count=20179 2261 17918
shrinkage=0.253694


Tree=6
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=249.642 15.4696 2.92473
threshold=0.99980113353686251 0.98569298860888022 0.99990079365079365
decision_type=2 2 2
left_child=1 -1 -2
right_child=2 -3 -4
leaf_value=-0.25263446526215178 -0.0035507900675803487 -0.10201954959764924 0.037321684190694482
leaf_weight=190.60750715434551 40.333420976996422 186.41950935125351 2004.6818205267191
leaf_count=1325 179 936 17739
internal_value=0 -0.741014 0.139048
internal_weight=0 377.027 2045.02
internal_count=20179 2261 17918
shrinkage=0.253694


Tree=7
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=169.987 11.1152
threshold=0.9982121572596041 0.99990079365079365
decision_type=2 2
left_child=-1 -2
right_child=1 -3
leaf_value=-0.19615789054037383 -0.043069247506608295 0.028802754918259385
leaf_weight=250.21898817270994 129.93249490857124 1949.7521076053381
leaf_count=1875 565 17739
internal_value=0 0.0905477
internal_weight=0 2079.68
internal_count=20179 18304
shrinkage=0.253694


Tree=8
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=112.307 6.42877 6.38337
threshold=0.9982121572596041 0.97183842176180868 0.99990079365079365
decision_type=2 2 2
left_child=1 -1 -2
right_child=2 -3 -4
leaf_value=-0.22959537631516913 -0.032785315724712938 -0.068748713309142701 0.022071285492765144
leaf_weight=119.13858236372471 129.65149366855621 99.307088911533356 1907.9913998991251
leaf_count=1139 565 736 17739
internal_value=0 -0.67847 0.0679734
internal_weight=0 218.446 2037.64
internal_count=20179 1875 18304
shrinkage=0.253694


Tree=9
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=79.0666 8.20159
threshold=0.99254269610309087 0.99990079365079365
decision_type=2 2
left_child=-1 -2
right_child=1 -3
leaf_value=-0.18389863713783922 -0.038492633770623154 0.016816609004074497
leaf_weight=138.51884593814611 185.17426434159279 1876.3697953149676
leaf_count=1510 930 17739
internal_value=0 0.0414402
internal_weight=0 2061.54
internal_count=20179 18669
shrinkage=0.253694


Tree=10
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=56.5695 6.95357
threshold=0.97183842176180868 0.99990079365079365
decision_type=2 2
left_child=-1 -2
right_child=1 -3
leaf_value=-0.19776391584840433 -0.035898567340432896 0.012755096963640257
leaf_weight=85.527907192707062 217.27691228687763 1852.5005107820034
leaf_count=1139 1301 17739
internal_value=0 0.0249388
internal_weight=0 2069.78
internal_count=20179 19040
shrinkage=0.253694


Tree=11
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=40.0174 3.67656
threshold=0.97183842176180868 0.99990079365079365
decision_type=2 2
left_child=-1 -2
right_child=1 -3
leaf_value=-0.18097106480750982 -0.027453007236968637 0.0096403517415896891
leaf_weight=72.286990508437157 214.65250280499458 1834.5264839529991
leaf_count=1139 1301 17739
internal_value=0 0.0175236
internal_weight=0 2049.18
internal_count=20179 19040
shrinkage=0.253694


Tree=12
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=28.0973 1.7591
threshold=0.97183842176180868 0.99990079365079365
decision_type=2 2
left_child=-1 -2
right_child=1 -3
leaf_value=-0.16417641607746813 -0.020934210748810619 0.0072663465110196215
leaf_weight=61.704010363668203 212.6226559355855 1821.0167551785707
leaf_count=1139 1301 17739
internal_value=0 0.0118967
internal_weight=0 2033.64
internal_count=20179 19040
shrinkage=0.253694


Tree=13
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=19.8279 0.959322
threshold=0.94991531433199572 0.99990079365079365
decision_type=2 2
left_child=-1 -2
right_child=1 -3
leaf_value=-0.16116274639848122 -0.01802985753016919 0.0054654815991239695
leaf_weight=44.30275684222579 220.04430298507214 1810.8767268955708
leaf_count=947 1493 17739
internal_value=0 0.00641329
internal_weight=0 2030.92
internal_count=20179 19232
shrinkage=0.253694


Tree=14
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=13.9451 0.159226
threshold=0.94991531433199572 0.99990079365079365
decision_type=2 2
left_child=-1 -2
right_child=1 -3
leaf_value=-0.14545017759916695 -0.013703802545201373 0.0041043909119972567
leaf_weight=38.266905508935452 218.54686445742846 1803.2742829173803
leaf_count=947 1493 17739
internal_value=0 0.00352294
internal_weight=0 2021.82
internal_count=20179 19232
shrinkage=0.253694


Tree=15
num_leaves=2
num_cat=0
split_feature=0
split_gain=9.7258
threshold=0.97183842176180868
decision_type=2
left_child=-1
right_child=-2
leaf_value=-0.11818090022019932 0.00061883839490289188
leaf_weight=42.196434009820223 2006.255480453372
leaf_count=1139 19040
internal_value=0
internal_weight=0
internal_count=20179
shrinkage=0.253694


Tree=16
num_leaves=2
num_cat=0
split_feature=0
split_gain=7.1761
threshold=0.9982121572596041
decision_type=2
left_child=-1
right_child=-2
leaf_value=-0.063789263744402511 0.0022570253992911138
leaf_weight=119.16151709109545 1924.104017406702
leaf_count=1875 18304
internal_value=0
internal_weight=0
internal_count=20179
shrinkage=0.253694


Tree=17
num_leaves=2
num_cat=0
split_feature=0
split_gain=5.58707
threshold=0.94991531433199572
decision_type=2
left_child=-1
right_child=-2
leaf_value=-0.10907001474775206 0.0004660517733198908
leaf_weight=28.243536002933979 2006.0411375388503
leaf_count=947 19232
internal_value=0
internal_weight=0
internal_count=20179
shrinkage=0.253694


Tree=18
num_leaves=2
num_cat=0
split_feature=0
split_gain=3.67054
threshold=0.94991531433199572
decision_type=2
left_child=-1
right_child=-2
leaf_value=-0.095495326339983921 0.00034834075058677816
leaf_weight=25.486886428669095 2005.4336205199361
leaf_count=947 19232
internal_value=0
internal_weight=0
internal_count=20179
shrinkage=0.253694


Tree=19
num_leaves=2
num_cat=0
split_feature=0
split_gain=2.30305
threshold=0.94991531433199572
decision_type=2
left_child=-1
right_child=-2
leaf_value=-0.082983315354841602 0.00026035713778579678
leaf_weight=23.282849848270416 2004.979699447751
leaf_count=947 19232
internal_value=0
internal_weight=0
internal_count=20179
shrinkage=0.253694


Tree=20
num_leaves=2
num_cat=0
split_feature=0
split_gain=1.74777
threshold=0.9982121572596041
decision_type=2
left_child=-1
right_child=-2
leaf_value=-0.038624327519308575 0.0014192031148284896
leaf_weight=106.6096428912133 1919.5459945499897
leaf_count=1875 18304
internal_value=0
internal_weight=0
internal_count=20179
shrinkage=0.253694


Tree=21
num_leaves=2
num_cat=0
split_feature=0
split_gain=0.948164
threshold=0.94991531433199572
decision_type=2
left_child=-1
right_child=-2
leaf_value=-0.066195171786071372 0.00026024484590016098
leaf_weight=20.736467184498906 2000.2457632794976
leaf_count=947 19232
internal_value=0
internal_weight=0
internal_count=20179
shrinkage=0.253694


Tree=22
num_leaves=2
num_cat=0
split_feature=0
split_gain=0.462266
threshold=0.99254269610309087
decision_type=2
left_child=-1
right_child=-2
leaf_value=-0.038077324107392452 0.00063355720966959492
leaf_weight=54.966307735070586 1964.4044191390276
leaf_count=1510 18669
internal_value=0
internal_weight=0
internal_count=20179
shrinkage=0.253694


Tree=23
num_leaves=2
num_cat=0
split_feature=0
split_gain=0.0998015
threshold=0.94991531433199572
decision_type=2
left_child=-1
right_child=-2
leaf_value=-0.051042117298845807 0.00020555916442422531
leaf_weight=18.766288246959448 1997.9275843575597
leaf_count=947 19232
internal_value=0
internal_weight=0
internal_count=20179
shrinkage=0.253694


end of trees

feature importances:
Gnomad_exomes_ASJ_wtf=45

parameters:
[boosting: gbdt]
[objective: binary]
[metric: auc]
[tree_learner: serial]
[device_type: cpu]
[data: ]
[valid: ]
[num_iterations: 100]
[learning_rate: 0.253694]
[num_leaves: 7]
[num_threads: -1]
[max_depth: 181]
[min_data_in_leaf: 97]
[min_sum_hessian_in_leaf: 16.9667]
[bagging_fraction: 0.775331]
[pos_bagging_fraction: 1]
[neg_bagging_fraction: 1]
[bagging_freq: 89]
[bagging_seed: 3]
[feature_fraction: 0.355777]
[feature_fraction_bynode: 1]
[feature_fraction_seed: 2]
[early_stopping_round: 0]
[first_metric_only: 0]
[max_delta_step: 0]
[lambda_l1: 9.95341]
[lambda_l2: 7.0806]
[min_gain_to_split: 0.947787]
[drop_rate: 0.1]
[max_drop: 50]
[skip_drop: 0.5]
[xgboost_dart_mode: 0]
[uniform_drop: 0]
[drop_seed: 4]
[top_rate: 0.2]
[other_rate: 0.1]
[min_data_per_group: 100]
[max_cat_threshold: 32]
[cat_l2: 10]
[cat_smooth: 10]
[max_cat_to_onehot: 4]
[top_k: 20]
[monotone_constraints: ]
[feature_contri: ]
[forcedsplits_filename: ]
[forcedbins_filename: ]
[refit_decay_rate: 0.9]
[cegb_tradeoff: 1]
[cegb_penalty_split: 0]
[cegb_penalty_feature_lazy: ]
[cegb_penalty_feature_coupled: ]
[verbosity: -1]
[max_bin: 15]
[max_bin_by_feature: ]
[min_data_in_bin: 3]
[bin_construct_sample_cnt: 200000]
[histogram_pool_size: -1]
[data_random_seed: 1]
[output_model: LightGBM_model.txt]
[snapshot_freq: -1]
[input_model: ]
[output_result: LightGBM_predict_result.txt]
[initscore_filename: ]
[valid_data_initscores: ]
[pre_partition: 0]
[enable_bundle: 1]
[max_conflict_rate: 0]
[is_enable_sparse: 1]
[sparse_threshold: 0.8]
[use_missing: 1]
[zero_as_missing: 0]
[two_round: 0]
[save_binary: 0]
[header: 0]
[label_column: ]
[weight_column: ]
[group_column: ]
[ignore_column: ]
[categorical_feature: ]
[predict_raw_score: 0]
[predict_leaf_index: 0]
[predict_contrib: 0]
[num_iteration_predict: -1]
[pred_early_stop: 0]
[pred_early_stop_freq: 10]
[pred_early_stop_margin: 10]
[convert_model_language: ]
[convert_model: gbdt_prediction.cpp]
[num_class: 1]
[is_unbalance: 0]
[scale_pos_weight: 1]
[sigmoid: 1]
[boost_from_average: 1]
[reg_sqrt: 0]
[alpha: 0.9]
[fair_c: 1]
[poisson_max_delta_step: 0.7]
[tweedie_variance_power: 1.5]
[max_position: 20]
[lambdamart_norm: 1]
[label_gain: ]
[metric_freq: 1]
[is_provide_training_metric: 0]
[eval_at: ]
[multi_error_top_k: 1]
[num_machines: 1]
[local_listen_port: 12400]
[time_out: 120]
[machine_list_filename: ]
[machines: ]
[gpu_platform_id: -1]
[gpu_device_id: -1]
[gpu_use_dp: 0]

end of parameters

pandas_categorical:[]
q(X   networkq)�X   _Booster__need_reload_eval_infoq*�X   _train_data_nameq+X   trainingq,X   _Booster__attrq-}q.X   _Booster__set_objective_to_noneq/�X   best_iterationq0K X
   best_scoreq1ccollections
defaultdict
q2ccollections
OrderedDict
q3�q4Rq5X   name_valid_setsq6]q7X   _Booster__num_datasetq8K X   _Booster__init_predictorq9NX   _Booster__num_classq:KX   _Booster__inner_predict_bufferq;]q<X   _Booster__is_predicted_cur_iterq=]q>X   _Booster__num_inner_evalq?KX   _Booster__name_inner_evalq@]qAX   aucqBaX"   _Booster__higher_better_inner_evalqC]qD�aX   pandas_categoricalqE]qFX   paramsqG}qH(hhhG?�      h	h
hC�/Ϡ�<�?qI�qJRqKhK�hKhG?PbM���hG        hJ����hKhNhG        hG        hG?�      hJ@ hK X   max_binqLKX   min_data_in_leafqMKaX   min_sum_hessian_in_leafqNh
hC��aw�0@qO�qPRqQX   bagging_fractionqRh
hC��L����?qS�qTRqUX   bagging_freqqVKYX   feature_fractionqWh
hCyYE-��?qX�qYRqZX	   lambda_l1q[h
hC���%�#@q\�q]Rq^X	   lambda_l2q_h
hCF�Z��R@q`�qaRqbX   min_gain_to_splitqch
hC�ET�?qd�qeRqfX   verboseqgJ����hhX   metricqh]qiX   aucqjauubX   _evals_resultqkNX   _best_scoreqlh5X   _best_iterationqmNX   _other_paramsqn}qo(hLKhMKahNh
hC��aw�0@qp�qqRqrhRh
hC��L����?qs�qtRquhVKYhWh
hCyYE-��?qv�qwRqxh[h
hC���%�#@qy�qzRq{h_h
hCF�Z��R@q|�q}Rq~hch
hC�ET�?q�q�Rq�hhhjuX
   _objectiveq�hX   class_weightq�NX   _class_weightq�NX
   _class_mapq�}q�(h
hX   i8q�K K�q�Rq�(KhNNNJ����J����K tq�bC        q��q�Rq�h
h�C        q��q�Rq�h
h�C       q��q�Rq�h
h�C       q��q�Rq�uX   _n_featuresq�KX   _classesq�cnumpy.core.multiarray
_reconstruct
q�cnumpy
ndarray
q�K �q�Cbq��q�Rq�(KK�q�h��C               q�tq�bX
   _n_classesq�KhLKhMKahNhrhRhuhVKYhWhxh[h{h_h~hch�hhhjX   _leq�csklearn.preprocessing._label
LabelEncoder
q�)�q�}q�(X   classes_q�h�X   _sklearn_versionq�X   0.22.1q�ubX   _fobjq�Nub.