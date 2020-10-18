�clightgbm.sklearn
LGBMClassifier
q )�q}q(X   boosting_typeqX   gbdtqX	   objectiveqX   binaryqX
   num_leavesqKX	   max_depthqKX   learning_rateq	cnumpy.core.multiarray
scalar
q
cnumpy
dtype
qX   f8qK K�qRq(KX   <qNNNJ����J����K tqbC�-�	�u�?q�qRqX   n_estimatorsqKdX   subsample_for_binqJ@ X   min_split_gainqG        X   min_child_weightqG?PbM���X   min_child_samplesqKX	   subsampleqG?�      X   subsample_freqqK X   colsample_bytreeqG?�      X	   reg_alphaqG        X
   reg_lambdaqG        X   random_stateqNX   n_jobsqJ����X   silentq �X   importance_typeq!X   splitq"X   _Boosterq#clightgbm.basic
Booster
q$)�q%}q&(X   handleq'X�H  tree
version=v3
num_class=1
num_tree_per_iteration=1
label_index=0
max_feature_idx=0
objective=binary sigmoid:1
feature_names=Gnomad_exomes_FIN_hetf
feature_infos=[0:0.69152095018145832]
tree_sizes=539 556 559 647 554 556 555 561 649 643 451 451 454 457 453 453 457 453 326 327 324 324 345 347 346 348 325 329 348 349 330 326 349 349 436 436

Tree=0
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=6687.59 185.96 119.375
threshold=0.00018527096036708332 0.0020472520234926556 1.0000000180025095e-35
decision_type=2 2 2
left_child=2 -2 -1
right_child=1 -3 -4
leaf_value=1.4905255223818459 0.8797316486508433 0.27026619225120374 1.199024647373981
leaf_weight=2975.930894985795 71.199175044894218 285.93199454247952 74.767243042588234
leaf_count=18349 439 1763 461
internal_value=0 -4.03606 0.483374
internal_weight=0 357.131 3050.7
internal_count=21012 2202 18810
shrinkage=1


Tree=1
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=2503.85 57.2847 15.1066
threshold=9.262260937755171e-05 0.0060057285923407798 1.0000000180025095e-35
decision_type=2 2 2
left_child=2 -2 -1
right_child=1 -3 -4
leaf_value=0.10485746189295614 -0.29757713977981898 -0.52802490417252401 -0.013440389447971633
leaf_weight=2753.1695533245802 187.86598606407642 377.56327569484711 40.403189793229103
leaf_count=18349 898 1538 227
internal_value=0 -1.87368 0.417121
internal_weight=0 565.429 2793.57
internal_count=21012 2436 18576
shrinkage=0.245781


Tree=2
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=1495.46 27.8408 22.3785
threshold=9.9448437238843426e-05 0.00064680064818120684 1.0000000180025095e-35
decision_type=2 2 2
left_child=2 -2 -1
right_child=1 -3 -4
leaf_value=0.087046556247701815 -0.14213249735651554 -0.38187373915985012 -0.039600414582792162
leaf_weight=2573.5141557008028 72.213396281003952 485.22827272117138 64.537037715315819
leaf_count=18357 326 1984 345
internal_value=0 -1.45098 0.338498
internal_weight=0 557.442 2638.05
internal_count=21012 2310 18702
shrinkage=0.245781


Tree=3
num_leaves=5
num_cat=0
split_feature=0 0 0 0
split_gain=921.322 25.5287 15.0664 0.0488791
threshold=9.9448437238843426e-05 0.0060057285923407798 1.0000000180025095e-35 9.2408631163409786e-05
decision_type=2 2 2 2
left_child=2 -2 -1 -4
right_child=1 -3 3 -5
leaf_value=0.070338686637966608 -0.17462623299760652 -0.33533492637337137 -0 -0.040257544444523057
leaf_weight=2426.7207724303007 184.2209325581789 347.48047161102295 23.611472859978676 42.186449810862541
leaf_count=18357 773 1537 129 216
internal_value=0 -1.16666 0.272248 -0.124162
internal_weight=0 531.701 2492.52 65.7979
internal_count=21012 2310 18702 345
shrinkage=0.245781


Tree=4
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=609.3 44.9615 6.53654
threshold=0.00064680064818120684 1.0000000180025095e-35 0.01568289706692064
decision_type=2 2 2
left_child=1 -1 -2
right_child=2 -3 -4
leaf_value=0.059552454612091001 -0.17114954873375895 -0.071528456965659903 -0.30913232967215382
leaf_weight=2305.0539468675852 151.01113067567348 151.33949330449104 265.10030889511108
leaf_count=18311 667 701 1333
internal_value=0 0.205732 -1.08987
internal_weight=0 2456.39 416.111
internal_count=21012 19012 2000
shrinkage=0.245781


Tree=5
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=407.37 27.5279 10.4433
threshold=0.00064680064818120684 1.0000000180025095e-35 0.01568289706692064
decision_type=2 2 2
left_child=1 -1 -2
right_child=2 -3 -4
leaf_value=0.047232249342359159 -0.13802137560128649 -0.053709889853799715 -0.28065605670536004
leaf_weight=2209.3603670746088 144.47686973214149 154.62602652609348 226.38056200742722
leaf_count=18311 667 701 1333
internal_value=0 0.161723 -0.953111
internal_weight=0 2363.99 370.857
internal_count=21012 19012 2000
shrinkage=0.245781


Tree=6
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=271.73 22.0488 3.28667
threshold=9.9448437238843426e-05 0.01568289706692064 1.0000000180025095e-35
decision_type=2 2 2
left_child=2 -2 -1
right_child=1 -3 -4
leaf_value=0.033388997195035755 -0.10979052613165073 -0.25666594567928119 -0.005479928231712632
leaf_weight=2144.3209803774953 215.67732453346252 185.81974343955517 69.094742804765701
leaf_count=18391 984 1299 338
internal_value=0 -0.757331 0.12783
internal_weight=0 401.497 2213.42
internal_count=21012 2283 18729
shrinkage=0.245781


Tree=7
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=179.026 24.3396 2.21092
threshold=9.9448437238843426e-05 0.022900119404348076 1.0000000180025095e-35
decision_type=2 2 2
left_child=2 -2 -1
right_child=1 -3 -4
leaf_value=0.025943693456194593 -0.09013219470313244 -0.24547358053223445 -0.0042862463318159129
leaf_weight=2092.3831334412098 224.35770070552826 142.28460171818733 69.253212377429008
leaf_count=18391 1096 1187 338
internal_value=0 -0.648845 0.0984772
internal_weight=0 366.642 2161.64
internal_count=21012 2283 18729
shrinkage=0.245781


Tree=8
num_leaves=5
num_cat=0
split_feature=0 0 0 0
split_gain=133.551 8.12224 0.673178 0.0956559
threshold=0.0060057285923407798 1.0000000180025095e-35 0.0012012567095732176 0.00018527096036708332
decision_type=2 2 2 2
left_child=1 -1 3 -3
right_child=-2 2 -4 -5
leaf_value=0.019164960332056412 -0.20946824574606071 -0.013753299639041795 -0.053184243372713449 0
leaf_weight=2051.3090020418167 165.01457254588604 96.655569866299629 67.369343519210815 79.170107454061508
leaf_count=18380 1530 449 330 323
internal_value=0 0.0552784 -0.10544 -0.0179235
internal_weight=0 2294.5 243.195 175.826
internal_count=21012 19482 1102 772
shrinkage=0.245781


Tree=9
num_leaves=5
num_cat=0
split_feature=0 0 0 0
split_gain=94.4137 5.17016 0.367546 0.0219732
threshold=0.0060057285923407798 1.0000000180025095e-35 0.00064680064818120684 0.00018527096036708332
decision_type=2 2 2 2
left_child=1 -1 3 -3
right_child=-2 2 -4 -5
leaf_value=0.0147176990336254 -0.19126655208605081 -0.010624805749474166 -0.038071744218739323 0
leaf_weight=2022.1839734911919 140.55671190470457 97.099017307162285 91.644626155495644 53.452253371477127
leaf_count=18380 1530 449 439 214
internal_value=0 0.0409896 -0.0864826 -0
internal_weight=0 2264.38 242.196 150.551
internal_count=21012 19482 1102 663
shrinkage=0.245781


Tree=10
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=69.9937 6.81279
threshold=0.010201984979849312 9.2408631163409786e-05
decision_type=2 2
left_child=1 -1
right_child=-2 -3
leaf_value=0.01302021347237319 -0.18703314047368641 -0.02957176462068678
leaf_weight=2021.8183720707893 107.84902316331863 235.2427020445466
leaf_count=18480 1414 1118
internal_value=0 0.0315341
internal_weight=0 2257.06
internal_count=21012 19598
shrinkage=0.245781


Tree=11
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=51.7656 4.62411
threshold=0.01568289706692064 1.0000000180025095e-35
decision_type=2 2
left_child=1 -1
right_child=-2 -3
leaf_value=0.010180394628077186 -0.18464266395513723 -0.02310531219017101
leaf_weight=1978.0000956952572 80.703087605535984 270.0353767350316
leaf_count=18357 1299 1356
internal_value=0 0.0218635
internal_weight=0 2248.04
internal_count=21012 19713
shrinkage=0.245781


Tree=12
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=39.0486 1.33543
threshold=0.022900119404348076 9.9448437238843426e-05
decision_type=2 2
left_child=1 -1
right_child=-2 -3
leaf_value=0.0089200367203196955 -0.17668217300133554 -0.010142025485529385
leaf_weight=2033.7037523165345 63.405565861612558 201.57403553277254
leaf_count=18688 1215 1109
internal_value=0 0.026175
internal_weight=0 2235.28
internal_count=21012 19797
shrinkage=0.245781


Tree=13
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=28.4522 0.765062
threshold=0.022900119404348076 9.9448437238843426e-05
decision_type=2 2
left_child=1 -1
right_child=-2 -3
leaf_value=0.0067898256008959877 -0.16249617585160886 -0.0077879700603626732
leaf_weight=2020.2456248253584 54.100402090698481 200.96213372051716
leaf_count=18688 1215 1109
internal_value=0 0.0191499
internal_weight=0 2221.21
internal_count=21012 19797
shrinkage=0.245781


Tree=14
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=17.6093 1.57151
threshold=0.022900119404348076 0.00064680064818120684
decision_type=2 2
left_child=1 -1
right_child=-2 -3
leaf_value=0.0032997191752069852 -0.13901718649165135 -0.025640442947893232
leaf_weight=2092.0518027096987 46.368837088346481 119.0804528594017
leaf_count=19025 1208 779
internal_value=0 0.0036533
internal_weight=0 2211.13
internal_count=21012 19804
shrinkage=0.245781


Tree=15
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=13.0147 0.96372
threshold=0.036016846889301853 0.00064680064818120684
decision_type=2 2
left_child=1 -1
right_child=-2 -3
leaf_value=0.0024992307472445534 -0.13094686765035068 -0.02074341764072677
leaf_weight=2087.1041375175118 37.262592494487762 121.0509730912745
leaf_count=19025 1104 883
internal_value=0 0.00162955
internal_weight=0 2208.16
internal_count=21012 19908
shrinkage=0.245781


Tree=16
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=12.2416 0.122293
threshold=0.0060057285923407798 0.0012012567095732176
decision_type=2 2
left_child=1 -1
right_child=-2 -3
leaf_value=0.0033345193174224028 -0.099990835265471312 -0.0083605647108574051
leaf_weight=2110.7820084840059 64.487413834780455 63.854858100414276
leaf_count=19130 1540 342
internal_value=0 0.0089667
internal_weight=0 2174.64
internal_count=21012 19472
shrinkage=0.245781


Tree=17
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=8.58558 0.308097
threshold=0.01568289706692064 0.0012012567095732176
decision_type=2 2
left_child=1 -1
right_child=-2 -3
leaf_value=0.002525503617131112 -0.10428274550855235 -0.014059046326554435
leaf_weight=2105.8269284814596 38.955283092334867 83.683950759470463
leaf_count=19130 1320 562
internal_value=0 0.00441966
internal_weight=0 2189.51
internal_count=21012 19692
shrinkage=0.245781


Tree=18
num_leaves=2
num_cat=0
split_feature=0
split_gain=6.87987
threshold=0.022900119404348076
decision_type=2
left_child=-1
right_child=-2
leaf_value=-0 -0.10525364204158225
leaf_weight=2189.9141000881791 29.708473324775696
leaf_count=19817 1195
internal_value=0
internal_weight=0
internal_count=21012
shrinkage=0.245781


Tree=19
num_leaves=2
num_cat=0
split_feature=0
split_gain=4.97866
threshold=0.022900119404348076
decision_type=2
left_child=-1
right_child=-2
leaf_value=-0 -0.093348964945362226
leaf_weight=2189.9141000881791 26.877448696643114
leaf_count=19817 1195
internal_value=0
internal_weight=0
internal_count=21012
shrinkage=0.245781


Tree=20
num_leaves=2
num_cat=0
split_feature=0
split_gain=2.9207
threshold=0.036016846889301853
decision_type=2
left_child=-1
right_child=-2
leaf_value=0 -0.077086226960848714
leaf_weight=2190.162201764062 22.149050725623965
leaf_count=19925 1087
internal_value=0
internal_weight=0
internal_count=21012
shrinkage=0.245781


Tree=21
num_leaves=2
num_cat=0
split_feature=0
split_gain=2.0958
threshold=0.036016846889301853
decision_type=2
left_child=-1
right_child=-2
leaf_value=0 -0.067691846237868947
leaf_weight=2190.162201764062 20.569274136796594
leaf_count=19925 1087
internal_value=0
internal_weight=0
internal_count=21012
shrinkage=0.245781


Tree=22
num_leaves=2
num_cat=0
split_feature=0
split_gain=3.02597
threshold=0.022900119404348076
decision_type=2
left_child=-1
right_child=-2
leaf_value=0.0006599677985893147 -0.077735758791806636
leaf_weight=2184.8436221107841 22.33814281411469
leaf_count=19785 1227
internal_value=0
internal_weight=0
internal_count=21012
shrinkage=0.245781


Tree=23
num_leaves=2
num_cat=0
split_feature=0
split_gain=2.16693
threshold=0.022900119404348076
decision_type=2
left_child=-1
right_child=-2
leaf_value=0.00049872327409338207 -0.068287296166245293
leaf_weight=2183.8960262052715 20.725297832861543
leaf_count=19785 1227
internal_value=0
internal_weight=0
internal_count=21012
shrinkage=0.245781


Tree=24
num_leaves=2
num_cat=0
split_feature=0
split_gain=1.1877
threshold=0.0012012567095732176
decision_type=2
left_child=-1
right_child=-2
leaf_value=0.00040338489533014847 -0.026758715927005527
leaf_weight=2099.310624435544 105.36520297825336
leaf_count=19142 1870
internal_value=0
internal_weight=0
internal_count=21012
shrinkage=0.245781


Tree=25
num_leaves=2
num_cat=0
split_feature=0
split_gain=0.66034
threshold=0.0012012567095732176
decision_type=2
left_child=-1
right_child=-2
leaf_value=0.00030479293165020904 -0.021114464906706803
leaf_weight=2098.7116609737277 103.54231200926006
leaf_count=19142 1870
internal_value=0
internal_weight=0
internal_count=21012
shrinkage=0.245781


Tree=26
num_leaves=2
num_cat=0
split_feature=0
split_gain=0.12492
threshold=0.01568289706692064
decision_type=2
left_child=-1
right_child=-2
leaf_value=-0 -0.023869155503570015
leaf_weight=2179.0616495087743 23.33276576641947
leaf_count=19697 1315
internal_value=0
internal_weight=0
internal_count=21012
shrinkage=0.245781


Tree=27
num_leaves=2
num_cat=0
split_feature=0
split_gain=0.049424
threshold=9.262260937755171e-05
decision_type=2
left_child=-1
right_child=-2
leaf_value=0 -0.0075072446718502965
leaf_weight=1962.2005877718329 239.66778599843383
leaf_count=18566 2446
internal_value=0
internal_weight=0
internal_count=21012
shrinkage=0.245781


Tree=28
num_leaves=2
num_cat=0
split_feature=0
split_gain=0.361669
threshold=0.0012012567095732176
decision_type=2
left_child=-1
right_child=-2
leaf_value=0.00016425199120878432 -0.017109617597311711
leaf_weight=2097.9066784158349 101.4950684197247
leaf_count=19138 1874
internal_value=0
internal_weight=0
internal_count=21012
shrinkage=0.245781


Tree=29
num_leaves=2
num_cat=0
split_feature=0
split_gain=0.151083
threshold=0.0012012567095732176
decision_type=2
left_child=-1
right_child=-2
leaf_value=0.00012408800643224645 -0.013424959418738051
leaf_weight=2097.6630469784141 100.35297907330096
leaf_count=19138 1874
internal_value=0
internal_weight=0
internal_count=21012
shrinkage=0.245781


Tree=30
num_leaves=2
num_cat=0
split_feature=0
split_gain=0.272631
threshold=1.0000000180025095e-35
decision_type=2
left_child=-1
right_child=-2
leaf_value=0 -0.0098280968981363341
leaf_weight=1915.9693440422416 285.51252831891179
leaf_count=18325 2687
internal_value=0
internal_weight=0
internal_count=21012
shrinkage=0.245781


Tree=31
num_leaves=2
num_cat=0
split_feature=0
split_gain=0.163023
threshold=0.01568289706692064
decision_type=2
left_child=-1
right_child=-2
leaf_value=-0 -0.02602746764482124
leaf_weight=2179.5113751962781 21.513823501765728
leaf_count=19698 1314
internal_value=0
internal_weight=0
internal_count=21012
shrinkage=0.245781


Tree=32
num_leaves=2
num_cat=0
split_feature=0
split_gain=0.666883
threshold=1.0000000180025095e-35
decision_type=2
left_child=-1
right_child=-2
leaf_value=0.0031212124555000816 -0.010755303595560895
leaf_weight=1917.1194484233856 275.71292837336659
leaf_count=18336 2676
internal_value=0
internal_weight=0
internal_count=21012
shrinkage=0.245781


Tree=33
num_leaves=2
num_cat=0
split_feature=0
split_gain=0.316465
threshold=1.0000000180025095e-35
decision_type=2
left_child=-1
right_child=-2
leaf_value=0.002364356168740056 -0.0082170153210416076
leaf_weight=1912.5588684082031 275.20776311401278
leaf_count=18336 2676
internal_value=0
internal_weight=0
internal_count=21012
shrinkage=0.245781


Tree=34
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=0.611404 0.455523
threshold=0.01568289706692064 0.00064680064818120684
decision_type=2 2
left_child=1 -1
right_child=-2 -3
leaf_value=0 -0.038589630350509621 0.019053719863461341
leaf_weight=2065.7145905792713 21.046428102999926 101.59867696464062
leaf_count=19017 1325 670
internal_value=0 0.00401187
internal_weight=0 2167.31
internal_count=21012 19687
shrinkage=0.245781


Tree=35
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=0.381407 0.198159
threshold=0.01568289706692064 0.00064680064818120684
decision_type=2 2
left_child=1 -1
right_child=-2 -3
leaf_value=0 -0.033006188503945148 0.014615861576488635
leaf_weight=2065.7145905792713 20.280867731198668 102.62249217927456
leaf_count=19017 1325 670
internal_value=0 0.00311664
internal_weight=0 2168.34
internal_count=21012 19687
shrinkage=0.245781


end of trees

feature importances:
Gnomad_exomes_FIN_hetf=69

parameters:
[boosting: gbdt]
[objective: binary]
[metric: auc]
[tree_learner: serial]
[device_type: cpu]
[data: ]
[valid: ]
[num_iterations: 100]
[learning_rate: 0.245781]
[num_leaves: 7]
[num_threads: -1]
[max_depth: 30]
[min_data_in_leaf: 60]
[min_sum_hessian_in_leaf: 19.8892]
[bagging_fraction: 0.807079]
[pos_bagging_fraction: 1]
[neg_bagging_fraction: 1]
[bagging_freq: 2]
[bagging_seed: 3]
[feature_fraction: 0.150362]
[feature_fraction_bynode: 1]
[feature_fraction_seed: 2]
[early_stopping_round: 0]
[first_metric_only: 0]
[max_delta_step: 0]
[lambda_l1: 6.68541]
[lambda_l2: 9.52552]
[min_gain_to_split: 0.177279]
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
[max_bin: 25]
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
hC�-�	�u�?qI�qJRqKhKhKhG?PbM���hG        hJ����hKhNhG        hG        hG?�      hJ@ hK X   max_binqLKX   min_data_in_leafqMK<X   min_sum_hessian_in_leafqNh
hC�*��3@qO�qPRqQX   bagging_fractionqRh
hC�e~w���?qS�qTRqUX   bagging_freqqVKX   feature_fractionqWh
hCg�Q4?�?qX�qYRqZX	   lambda_l1q[h
hC��l�۽@q\�q]Rq^X	   lambda_l2q_h
hC,[E�#@q`�qaRqbX   min_gain_to_splitqch
hC������?qd�qeRqfX   verboseqgJ����hhX   metricqh]qiX   aucqjauubX   _evals_resultqkNX   _best_scoreqlh5X   _best_iterationqmNX   _other_paramsqn}qo(hLKhMK<hNh
hC�*��3@qp�qqRqrhRh
hC�e~w���?qs�qtRquhVKhWh
hCg�Q4?�?qv�qwRqxh[h
hC��l�۽@qy�qzRq{h_h
hC,[E�#@q|�q}Rq~hch
hC������?q�q�Rq�hhhjuX
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
   _n_classesq�KhLKhMK<hNhrhRhuhVKhWhxh[h{h_h~hch�hhhjX   _leq�csklearn.preprocessing._label
LabelEncoder
q�)�q�}q�(X   classes_q�h�X   _sklearn_versionq�X   0.22.1q�ubX   _fobjq�Nub.