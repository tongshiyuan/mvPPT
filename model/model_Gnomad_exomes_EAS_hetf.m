�clightgbm.sklearn
LGBMClassifier
q )�q}q(X   boosting_typeqX   gbdtqX	   objectiveqX   binaryqX
   num_leavesqKX	   max_depthqK�X   learning_rateq	cnumpy.core.multiarray
scalar
q
cnumpy
dtype
qX   f8qK K�qRq(KX   <qNNNJ����J����K tqbC����� �?q�qRqX   n_estimatorsqKdX   subsample_for_binqJ@ X   min_split_gainqG        X   min_child_weightqG?PbM���X   min_child_samplesqKX	   subsampleqG?�      X   subsample_freqqK X   colsample_bytreeqG?�      X	   reg_alphaqG        X
   reg_lambdaqG        X   random_stateqNX   n_jobsqJ����X   silentq �X   importance_typeq!X   splitq"X   _Boosterq#clightgbm.basic
Booster
q$)�q%}q&(X   handleq'X.=  tree
version=v3
num_class=1
num_tree_per_iteration=1
label_index=0
max_feature_idx=0
objective=binary sigmoid:1
feature_names=Gnomad_exomes_EAS_hetf
feature_infos=[0:0.53498694516971279]
tree_sizes=543 555 557 747 1070 1197 768 876 554 453 672 559 455 454 345 860 564 458 346 346 458

Tree=0
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=4507.15 479.41 127.554
threshold=0.0001112532684090029 0.0057222583290612019 1.0000000180025095e-35
decision_type=2 2 2
left_child=2 -2 -1
right_child=1 -3 -4
leaf_value=1.5008376513608568 0.75440789759764382 -0.0038970943827114635 1.1749406837809881
leaf_weight=2735.086305141449 173.70003752410412 170.61852425336838 104.60926629602909
leaf_count=16864 1071 1052 645
internal_value=0 -3.37503 0.41696
internal_weight=0 344.319 2839.7
internal_count=19632 2123 17509
shrinkage=1


Tree=1
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=1412.72 88.2149 29.5729
threshold=0.0001087961706967472 0.004364403181297616 1.0000000180025095e-35
decision_type=2 2 2
left_child=2 -2 -1
right_child=1 -3 -4
leaf_value=0.10663261500613316 -0.28835980628190366 -0.55510077102120292 -0.075648613267169554
leaf_weight=2513.867639541626 267.40034361183643 268.43775069713593 76.440634369850159
leaf_count=16864 1267 1077 424
internal_value=0 -1.43552 0.337446
internal_weight=0 535.838 2590.31
internal_count=19632 2344 17288
shrinkage=0.298879


Tree=2
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=720.124 67.9349 28.1008
threshold=0.00010979358937797895 0.0062572941900099687 1.0000000180025095e-35
decision_type=2 2 2
left_child=2 -2 -1
right_child=1 -3 -4
leaf_value=0.081371663691088864 -0.1968947341473547 -0.43242278842029863 -0.066221133637076143
leaf_weight=2345.295398235321 265.18826606869698 241.16458714008331 116.25834921002388
leaf_count=16864 1124 1042 602
internal_value=0 -1.05382 0.247651
internal_weight=0 506.353 2461.55
internal_count=19632 2166 17466
shrinkage=0.298879


Tree=3
num_leaves=6
num_cat=0
split_feature=0 0 0 0 0
split_gain=403.286 45.6795 27.988 2.33367 2.06031
threshold=0.0003262820173877495 1.0000000180025095e-35 0.0094793538131164954 0.00031791938993145231 0.00010873701983971457
decision_type=2 2 2 2 2
left_child=1 -1 -2 4 -3
right_child=2 3 -4 -5 -6
leaf_value=0.061048237203763275 -0.17432370935468564 -0 -0.37379996024888945 0.02760569648166953 -0.096905178113384893
leaf_weight=2219.5485410690308 149.65045221149921 43.944612100720406 195.8110032081604 16.941320791840553 194.00535632669926
leaf_count=16864 620 227 992 69 860
internal_value=0 0.15726 -0.985745 -0.239737 -0.272048
internal_weight=0 2474.44 345.461 254.891 237.95
internal_count=19632 18020 1612 1156 1087
shrinkage=0.298879


Tree=4
num_leaves=9
num_cat=0
split_feature=0 0 0 0 0 0 0 0
split_gain=236.023 30.2065 23.9993 0.4311 0.660187 3.109 2.93084 0.0733196
threshold=0.00032661949216767532 0.0062572941900099687 1.0000000180025095e-35 0.00031791938993145231 0.00021967167113694721 0.00020378814219052256 0.0001112532684090029 0.00021747403967942915
decision_type=2 2 2 2 2 2 2 2
left_child=2 -2 -1 4 5 6 -4 -7
right_child=1 -3 3 -5 -6 7 -8 -9
leaf_value=0.045117381361974213 -0.11915014238432534 -0.32795135599219305 -0.029492986294269184 0.0025876865715796678 -0.13800114169928768 0.073338609440785962 -0.13549998024350682 -0
leaf_weight=2127.2075943946838 135.69587253034115 169.78613856434822 132.81279969215393 22.862252220511436 19.852383136749268 16.626370877027512 45.660481214523315 27.545181602239609
leaf_count=16864 546 1042 645 93 80 67 184 111
internal_value=0 -0.810057 0.114434 -0.168398 -0.191152 -0.14501 -0.213439 0.0914788
internal_weight=0 305.482 2392.57 265.359 242.497 222.645 178.473 44.1716
internal_count=19632 1588 18044 1180 1087 1007 829 178
shrinkage=0.298879


Tree=5
num_leaves=10
num_cat=0
split_feature=0 0 0 0 0 0 0 0 0
split_gain=166.919 23.0566 3.2139 8.46882 0.611148 0.2291 0.135015 1.54343 0.115605
threshold=0.0034655634440324819 1.0000000180025095e-35 0.00022181556104367196 0.00011976770272650648 0.00011126564707594003 0.00076140751843204997 0.00043630022123747727 0.00023363255293883581 0.00081244921308450103
decision_type=2 2 2 2 2 2 2 2 2
left_child=1 -1 3 4 -3 6 7 -4 -7
right_child=-2 2 5 -5 -6 8 -8 -9 -10
leaf_value=0.034443342286838735 -0.29192803691822394 -0.033743239400544772 -0.15051020485827754 0.073549779391648085 -0.12255604959774559 0.03655120266925084 -0.15825605315622743 -0.0092360209557355195 -0.047805160739166457
leaf_weight=2054.9385000541806 154.92699933052063 145.93708914518356 11.496628940105436 63.301422536373138 24.491730242967606 4.4993345439434043 33.745009079575539 50.868509560823441 57.741459980607033
leaf_count=16821 1123 700 46 255 98 18 135 205 231
internal_value=0 0.0715027 -0.151507 -0.0544137 -0.180916 -0.275812 -0.3544 -0.187808 -0.105424
internal_weight=0 2447.02 392.081 233.73 170.429 158.351 96.1101 62.3651 62.2408
internal_count=19632 18509 1688 1053 798 635 386 251 249
shrinkage=0.298879


Tree=6
num_leaves=6
num_cat=0
split_feature=0 0 0 0 0
split_gain=109.567 13.3581 1.24837 1.0704 5.57889
threshold=0.004364403181297616 0.00010873701983971457 0.023590827662452001 0.00022181556104367196 0.00020378814219052256
decision_type=2 2 2 2 2
left_child=1 -1 -2 4 -3
right_child=2 3 -4 -5 -6
leaf_value=0.024242959023758219 -0.10672723135063146 -0.047728102394283142 -0.29983531813389513 -0.064784986093184255 0.073780410899677934
leaf_weight=2054.2905444726348 27.478311374783516 141.12570008635521 93.871363699436188 164.09219112992287 44.433132588863373
leaf_count=17070 214 623 884 660 181
internal_value=0 0.048726 -0.900216 -0.134693 -0.0479946
internal_weight=0 2403.94 121.35 349.651 185.559
internal_count=19632 18534 1098 1464 804
shrinkage=0.298879


Tree=7
num_leaves=7
num_cat=0
split_feature=0 0 0 0 0 0
split_gain=74.5375 6.56215 1.6369 0.26468 2.49259 0.108258
threshold=0.004364403181297616 1.0000000180025095e-35 0.027593311865070965 0.00022181556104367196 0.00011976770272650648 0.00076140751843204997
decision_type=2 2 2 2 2 2
left_child=1 -1 -2 4 -3 -5
right_child=2 3 -4 5 -6 -7
leaf_value=0.018082455117355317 -0.092604560421327656 -0.032354539344602327 -0.28309455940985784 -0.066826497022309014 0.044663397859553131 -0.007605418382322782
leaf_weight=1969.669174388051 27.66561534255743 173.29659017920494 72.106325775384903 95.62268628180027 62.425208851695061 68.005638748407364
leaf_count=16821 238 798 860 386 255 274
internal_value=0 0.0346162 -0.821206 -0.0878502 -0.030856 -0.155478
internal_weight=0 2369.02 99.7719 399.35 235.722 163.628
internal_count=19632 18534 1098 1713 1053 660
shrinkage=0.298879


Tree=8
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=51.9356 4.51677 0.213748
threshold=0.017396607385458052 0.00022181556104367196 0.015605486216445768
decision_type=2 2 2
left_child=1 -1 -3
right_child=-2 2 -4
leaf_value=0.010653494717889201 -0.25486566045845799 -0.047168241925496415 0.046013828437015118
leaf_weight=2180.8345576971769 63.392152078449726 180.93031158298254 1.7275648191571225
leaf_count=17874 935 804 19
internal_value=0 0.020958 -0.140242
internal_weight=0 2363.49 182.658
internal_count=19632 18697 823
shrinkage=0.298879


Tree=9
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=36.7251 2.62576
threshold=0.027593311865070965 0.00010873701983971457
decision_type=2 2
left_child=1 -1
right_child=-2 -3
leaf_value=0.0094751961268959104 -0.24939274229664057 -0.022019263386398473
leaf_weight=1978.2227356731892 45.538949817419052 373.56004170328379
leaf_count=17070 860 1702
internal_value=0 0.0141012
internal_weight=0 2351.78
internal_count=19632 18772
shrinkage=0.298879


Tree=10
num_leaves=5
num_cat=0
split_feature=0 0 0 0
split_gain=23.6531 6.80825 0.546719 0.314192
threshold=0.0052757306238804978 0.00011130279961011503 0.00022634686760885647 0.00032929049924350476
decision_type=2 2 2 2
left_child=1 -1 -3 -4
right_child=-2 2 3 -5
leaf_value=0.0090748379968199552 -0.18552142976767258 -0.076416233847360376 0.020693329684878713 -0.036208944156318056
leaf_weight=2065.6219162344933 55.670606888830662 95.987474679946899 31.542057946324348 121.0673925280571
leaf_count=17544 1072 392 127 497
internal_value=0 0.00942928 -0.153874 -0.0704592
internal_weight=0 2314.22 248.597 152.609
internal_count=19632 18560 1016 624
shrinkage=0.298879


Tree=11
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=16.5845 3.59692 0.282355
threshold=0.021443344871708517 0.00011126564707594003 0.015605486216445768
decision_type=2 2 2
left_child=1 -1 -3
right_child=-2 2 -4
leaf_value=0.0064476367717107548 -0.19545735172888754 -0.038973731451723809 0.051852195130198599
leaf_weight=2050.0929449647665 33.136464335024357 261.5113213583827 3.2178430482745162
leaf_count=17535 903 1146 48
internal_value=0 0.00482166 -0.115931
internal_weight=0 2314.82 264.729
internal_count=19632 18729 1194
shrinkage=0.298879


Tree=12
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=11.5388 1.28218
threshold=0.028475727435572449 0.00011130279961011503
decision_type=2 2
left_child=1 -1
right_child=-2 -3
leaf_value=0.0045431849155619687 -0.18379040970042707 -0.024343908641226642
leaf_weight=2042.9004647433758 25.182014876976609 263.99236802011728
leaf_count=17544 847 1241
internal_value=0 0.00321384
internal_weight=0 2306.89
internal_count=19632 18785
shrinkage=0.298879


Tree=13
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=7.93564 0.167079
threshold=0.028475727435572449 0.00011130279961011503
decision_type=2 2
left_child=1 -1
right_child=-2 -3
leaf_value=0.0032031531167244771 -0.16550786430423639 -0.017350472662136764
leaf_weight=2036.2792864441872 21.171695116907358 263.09909399226308
leaf_count=17544 847 1241
internal_value=0 0.001989
internal_weight=0 2299.38
internal_count=19632 18785
shrinkage=0.298879


Tree=14
num_leaves=2
num_cat=0
split_feature=0
split_gain=5.38599
threshold=0.028475727435572449
decision_type=2
left_child=-1
right_child=-2
leaf_value=0.00034189741483351609 -0.14794399797554025
leaf_weight=2294.0424034446478 18.083484753966331
leaf_count=18785 847
internal_value=0
internal_weight=0
internal_count=19632
shrinkage=0.298879


Tree=15
num_leaves=7
num_cat=0
split_feature=0 0 0 0 0 0
split_gain=3.93689 1.34859 1.10883 0.247979 1.05585 0.135525
threshold=0.028475727435572449 0.00020378814219052256 0.0003262820173877495 0.00022634686760885647 0.00021967167113694721 0.00017439845838332387
decision_type=2 2 2 2 2 2
left_child=1 5 3 4 -3 -1
right_child=-2 2 -4 -5 -6 -7
leaf_value=-0 -0.13355046205721091 0.058906180434448036 0.0050947995723480861 0.12444579083173543 -0.04826163809818735 -0.083013879510662594
leaf_weight=2067.2363989129663 15.987085215747355 42.465488716959953 143.52051306888461 22.374001145362854 11.484055399894713 6.4939820468425742
leaf_count=17672 863 172 761 90 48 26
internal_value=0 0.00764574 0.103612 0.230588 0.0842335 -0.00212632
internal_weight=0 2293.57 219.844 76.3235 53.9495 2073.73
internal_count=19632 18769 1071 310 220 17698
shrinkage=0.298879


Tree=16
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=2.53627 0.412349 0.184502
threshold=0.028475727435572449 0.00020378814219052256 0.0003262820173877495
decision_type=2 2 2
left_child=1 -1 -3
right_child=-2 2 -4
leaf_value=-0.00055826273439269782 -0.11809303259490023 0.052427685538823939 0.0036467464718251907
leaf_weight=2073.7028707042336 14.054707610979674 75.898399114608765 143.71292016282678
leaf_count=17698 863 310 761
internal_value=0 0.00550873 0.0797374
internal_weight=0 2293.31 219.611
internal_count=19632 18769 1071
shrinkage=0.298879


Tree=17
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=1.6906 0.0738898
threshold=0.0094793538131164954 0.00020378814219052256
decision_type=2 2
left_child=1 -1
right_child=-2 -3
leaf_value=-0.0003918862546384879 -0.087326262527875786 0.021225103963487358
leaf_weight=2074.5118878558278 22.573978129774332 209.43106587231159
leaf_count=17698 1027 907
internal_value=0 0.00485486
internal_weight=0 2283.94
internal_count=19632 18605
shrinkage=0.298879


Tree=18
num_leaves=2
num_cat=0
split_feature=0
split_gain=1.00459
threshold=0.044726099728915884
decision_type=2
left_child=-1
right_child=-2
leaf_value=0.00083869876117123089 -0.097421423282912339
leaf_weight=2295.3443910442293 10.448054549284278
leaf_count=18849 783
internal_value=0
internal_weight=0
internal_count=19632
shrinkage=0.298879


Tree=19
num_leaves=2
num_cat=0
split_feature=0
split_gain=0.48797
threshold=0.044726099728915884
decision_type=2
left_child=-1
right_child=-2
leaf_value=0.00058929132690993557 -0.086095506289732879
leaf_weight=2294.1539305113256 9.5020462526008469
leaf_count=18849 783
internal_value=0
internal_weight=0
internal_count=19632
shrinkage=0.298879


Tree=20
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=0.318028 0.525392
threshold=0.017396607385458052 0.00076140751843204997
decision_type=2 2
left_child=1 -1
right_child=-2 -3
leaf_value=-0.0020308612029076275 -0.073990865403233491 0.036572650950107324
leaf_weight=2202.0371161550283 13.601520801894365 85.330743901431561
leaf_count=18241 942 449
internal_value=0 -0.00087812
internal_weight=0 2287.37
internal_count=19632 18690
shrinkage=0.298879


end of trees

feature importances:
Gnomad_exomes_EAS_hetf=74

parameters:
[boosting: gbdt]
[objective: binary]
[metric: auc]
[tree_learner: serial]
[device_type: cpu]
[data: ]
[valid: ]
[num_iterations: 100]
[learning_rate: 0.298879]
[num_leaves: 20]
[num_threads: -1]
[max_depth: 184]
[min_data_in_leaf: 17]
[min_sum_hessian_in_leaf: 0.0494719]
[bagging_fraction: 0.754542]
[pos_bagging_fraction: 1]
[neg_bagging_fraction: 1]
[bagging_freq: 5]
[bagging_seed: 3]
[feature_fraction: 0.327177]
[feature_fraction_bynode: 1]
[feature_fraction_seed: 2]
[early_stopping_round: 0]
[first_metric_only: 0]
[max_delta_step: 0]
[lambda_l1: 1.58251]
[lambda_l2: 7.98636]
[min_gain_to_split: 0.972157]
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
[max_bin: 252]
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
hC����� �?qI�qJRqKhK�hKhG?PbM���hG        hJ����hKhNhG        hG        hG?�      hJ@ hK X   max_binqLK�X   min_data_in_leafqMKX   min_sum_hessian_in_leafqNh
hC�%��aT�?qO�qPRqQX   bagging_fractionqRh
hC��<4%�?qS�qTRqUX   bagging_freqqVKX   feature_fractionqWh
hC$��u��?qX�qYRqZX	   lambda_l1q[h
hCj����Q�?q\�q]Rq^X	   lambda_l2q_h
hC�qS��@q`�qaRqbX   min_gain_to_splitqch
hC��ݷ��?qd�qeRqfX   verboseqgJ����hhX   metricqh]qiX   aucqjauubX   _evals_resultqkNX   _best_scoreqlh5X   _best_iterationqmNX   _other_paramsqn}qo(hLK�hMKhNh
hC�%��aT�?qp�qqRqrhRh
hC��<4%�?qs�qtRquhVKhWh
hC$��u��?qv�qwRqxh[h
hCj����Q�?qy�qzRq{h_h
hC�qS��@q|�q}Rq~hch
hC��ݷ��?q�q�Rq�hhhjuX
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
   _n_classesq�KhLK�hMKhNhrhRhuhVKhWhxh[h{h_h~hch�hhhjX   _leq�csklearn.preprocessing._label
LabelEncoder
q�)�q�}q�(X   classes_q�h�X   _sklearn_versionq�X   0.22.1q�ubX   _fobjq�Nub.