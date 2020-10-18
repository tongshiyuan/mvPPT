�clightgbm.sklearn
LGBMClassifier
q )�q}q(X   boosting_typeqX   gbdtqX	   objectiveqX   binaryqX
   num_leavesqKX	   max_depthqK X   learning_rateq	cnumpy.core.multiarray
scalar
q
cnumpy
dtype
qX   f8qK K�qRq(KX   <qNNNJ����J����K tqbCҞ�l���?q�qRqX   n_estimatorsqKdX   subsample_for_binqJ@ X   min_split_gainqG        X   min_child_weightqG?PbM���X   min_child_samplesqKX	   subsampleqG?�      X   subsample_freqqK X   colsample_bytreeqG?�      X	   reg_alphaqG        X
   reg_lambdaqG        X   random_stateqNX   n_jobsqJ����X   silentq �X   importance_typeq!X   splitq"X   _Boosterq#clightgbm.basic
Booster
q$)�q%}q&(X   handleq'X�p  tree
version=v3
num_class=1
num_tree_per_iteration=1
label_index=0
max_feature_idx=0
objective=binary sigmoid:1
feature_names=Gnomad_exomes_ASJ_AF
feature_infos=[0:1]
tree_sizes=544 543 659 552 554 649 554 662 650 649 652 547 543 548 549 548 657 547 547 550 554 555 664 440 550 440 548 546 546 451 450 453 452 456 335 454 457 344 341 344 335 344 346 343 344 342 348 345 345 348 345 345 347 346

Tree=0
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=8006.66 122.402 105.372
threshold=9.9591700000000011e-05 1.0000000180025095e-35 0.0011922500000000004
decision_type=2 2 2
left_child=1 -1 -2
right_child=2 -3 -4
leaf_value=1.4368333302233305 0.98720929688976478 1.1988795473896785 0.76433702800112546
leaf_weight=3040.9670434892178 70.226065590977669 34.545385614037514 312.85468943417072
leaf_count=18750 433 213 1929
internal_value=0 0.537806 -4.29386
internal_weight=0 3075.51 383.081
internal_count=21325 18963 2362
shrinkage=1


Tree=1
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=4622.08 113.861 2.98807
threshold=1.0000000180025095e-35 0.0003976535000000001 0.01050795
decision_type=2 2 2
left_child=-1 -2 -3
right_child=1 2 -4
leaf_value=0.066799242272665552 -0.18881721888234507 -0.32721316412313595 -0.40145639170112069
leaf_weight=2909.2180542647839 80.632879391312599 181.86506080627441 279.19806241989136
leaf_count=18750 429 858 1288
internal_value=0 -2.66382 -2.86376
internal_weight=0 541.696 461.063
internal_count=21325 2575 2146
shrinkage=0.13089


Tree=2
num_leaves=5
num_cat=0
split_feature=0 0 0 0
split_gain=3198.01 69.126 2.80826 2.44217
threshold=1.0000000180025095e-35 0.0003976535000000001 0.0053693950000000008 9.9591700000000011e-05
decision_type=2 2 2 2
left_child=-1 3 -3 -2
right_child=1 2 -4 -5
leaf_value=0.060729050699203008 -0.098764231540566114 -0.24465412694644006 -0.30700133017999276 -0.18500341694048358
leaf_weight=2790.0458313524723 41.683231338858604 148.88583359122276 364.32053428888321 46.232703566551208
leaf_count=18750 213 637 1509 216
internal_value=0 -2.07386 -2.22155 -1.15166
internal_weight=0 601.122 513.206 87.9159
internal_count=21325 2575 2146 429
shrinkage=0.13089


Tree=3
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=2357.83 52.2111 7.06416
threshold=1.0000000180025095e-35 0.0011922500000000004 9.9591700000000011e-05
decision_type=2 2 2
left_child=-1 2 -2
right_child=1 -3 -4
leaf_value=0.055054499044746755 -0.082919132798531275 -0.24876034572343964 -0.17140873050735037
leaf_weight=2682.8468777239323 43.562277987599373 480.79910795390606 101.19902144372463
leaf_count=18750 213 1929 433
internal_value=0 -1.7313 -1.1327
internal_weight=0 625.56 144.761
internal_count=21325 2575 646
shrinkage=0.13089


Tree=4
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=1792.69 47.1114 8.18094
threshold=1.0000000180025095e-35 0.0025796250000000008 9.9591700000000011e-05
decision_type=2 2 2
left_child=-1 2 -2
right_child=1 -3 -4
leaf_value=0.049763734099770994 -0.070298404420347446 -0.22259508961240135 -0.15392629687212017
leaf_weight=2586.8651457130909 45.066386237740517 425.53888818621635 158.35493142902851
leaf_count=18750 213 1715 647
internal_value=0 -1.49877 -1.05122
internal_weight=0 628.96 203.421
internal_count=21325 2575 860
shrinkage=0.13089


Tree=5
num_leaves=5
num_cat=0
split_feature=0 0 0 0
split_gain=1389.83 44.5895 2.68844 2.55126
threshold=1.0000000180025095e-35 0.0011922500000000004 9.9591700000000011e-05 0.01050795
decision_type=2 2 2 2
left_child=-1 2 -2 -3
right_child=1 3 -4 -5
leaf_value=0.044847828983332769 -0.060027804392238565 -0.16185041358770688 -0.12183171150294084 -0.20957399240216953
leaf_weight=2501.2624450027943 46.276774227619171 157.09710785746574 107.03542570769787 308.49429202079773
leaf_count=18750 213 641 433 1288
internal_value=0 -1.3281 -0.808967 -1.48977
internal_weight=0 618.904 153.312 465.591
internal_count=21325 2575 646 1929
shrinkage=0.13089


Tree=6
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=1095.02 24.1322 23.53
threshold=9.9591700000000011e-05 1.0000000180025095e-35 0.0053693950000000008
decision_type=2 2 2
left_child=1 -1 -2
right_child=2 -3 -4
leaf_value=0.040298452866221843 -0.12326573046374881 -0.051534388826997152 -0.1883717635928665
leaf_weight=2425.1709692180157 208.63316358625889 47.255743369460106 344.38617545366287
leaf_count=18750 853 213 1509
internal_value=0 0.293426 -1.26011
internal_weight=0 2472.43 553.019
internal_count=21325 18963 2362
shrinkage=0.13089


Tree=7
num_leaves=5
num_cat=0
split_feature=0 0 0 0
split_gain=869.349 24.9011 18.7787 0.638915
threshold=9.9591700000000011e-05 0.0053693950000000008 1.0000000180025095e-35 0.0003976535000000001
decision_type=2 2 2 2
left_child=2 3 -1 -2
right_child=1 -3 -4 -5
leaf_value=0.036106490340896656 -0.070951681860645932 -0.17532446673686269 -0.044424487704043301 -0.11951337688086247
leaf_weight=2357.7268235385418 53.975097298622131 323.24247480928898 48.051439270377159 151.88535799086094
leaf_count=18750 216 1509 213 637
internal_value=0 -1.1504 0.262556 -0.831832
internal_weight=0 529.103 2405.78 205.86
internal_count=21325 2362 18963 853
shrinkage=0.13089


Tree=8
num_leaves=5
num_cat=0
split_feature=0 0 0 0
split_gain=695.762 26.4628 14.5703 1.85994
threshold=9.9591700000000011e-05 0.01050795 1.0000000180025095e-35 0.0003976535000000001
decision_type=2 2 2 2
left_child=2 3 -1 -2
right_child=1 -3 -4 -5
leaf_value=0.032261142750042757 -0.06212426441389101 -0.17069422475446658 -0.038416836536540808 -0.11290627618260332
leaf_weight=2298.0948910117149 53.989420294761658 254.55020773410797 48.70127959549427 193.95621874928474
leaf_count=18750 216 1288 213 858
internal_value=0 -1.05966 0.234291 -0.790525
internal_weight=0 502.496 2346.8 247.946
internal_count=21325 2362 18963 1074
shrinkage=0.13089


Tree=9
num_leaves=5
num_cat=0
split_feature=0 0 0 0
split_gain=559.16 27.4185 11.2622 2.0696
threshold=9.9591700000000011e-05 0.01050795 1.0000000180025095e-35 0.0011922500000000004
decision_type=2 2 2 2
left_child=2 3 -1 -2
right_child=1 -3 -4 -5
leaf_value=0.028749583157472433 -0.067675562851941123 -0.16238015447866574 -0.033303636190972701 -0.1096873373440652
leaf_weight=2245.4830817878246 106.15381360054016 234.10474729537964 49.234505921602249 134.98040702939034
leaf_count=18750 433 1288 213 641
internal_value=0 -0.98078 0.208523 -0.71065
internal_weight=0 475.239 2294.72 241.134
internal_count=21325 2362 18963 1074
shrinkage=0.13089


Tree=10
num_leaves=5
num_cat=0
split_feature=0 0 0 0
split_gain=450.817 28.162 8.66391 2.02218
threshold=9.9591700000000011e-05 0.01050795 1.0000000180025095e-35 0.0011922500000000004
decision_type=2 2 2 2
left_child=2 3 -1 -2
right_child=1 -3 -4 -5
leaf_value=0.025556936640350789 -0.059614198579411519 -0.15542643632263728 -0.028927055217719488 -0.1003176795162482
leaf_weight=2199.1514135152102 105.19944210350513 214.06398403644562 49.674062713980675 128.96910388767719
leaf_count=18750 433 1288 213 641
internal_value=0 -0.910882 0.185131 -0.64029
internal_weight=0 448.233 2248.83 234.169
internal_count=21325 2362 18963 1074
shrinkage=0.13089


Tree=11
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=366.216 19.4159 16.9789
threshold=0.0003976535000000001 1.0000000180025095e-35 0.01050795
decision_type=2 2 2
left_child=1 -1 -2
right_child=2 -3 -4
leaf_value=0.022666537774444731 -0.08479374303212156 -0.035068139240828397 -0.14947161378600282
leaf_weight=2158.4164816886187 173.99194985628128 103.40492968261242 194.88180601596832
leaf_count=18750 858 429 1288
internal_value=0 0.152036 -0.919789
internal_weight=0 2261.82 368.874
internal_count=21325 19179 2146
shrinkage=0.13089


Tree=12
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=298 17.9262 14.823
threshold=0.0003976535000000001 0.01050795 1.0000000180025095e-35
decision_type=2 2 2
left_child=2 -2 -1
right_child=1 -3 -4
leaf_value=0.020060262560475186 -0.076762712386462623 -0.14425694647187781 -0.030542624954949443
leaf_weight=2122.6529963314533 167.79056780040264 176.83249604701996 103.60474348068237
leaf_count=18750 858 1288 429
internal_value=0 -0.862359 0.134325
internal_weight=0 344.623 2226.26
internal_count=21325 2146 19179
shrinkage=0.13089


Tree=13
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=242.618 18.6239 11.2717
threshold=0.0003976535000000001 0.01050795 1.0000000180025095e-35
decision_type=2 2 2
left_child=2 -2 -1
right_child=1 -3 -4
leaf_value=0.017719008517400551 -0.069420613642358658 -0.13959092390904848 -0.026620479453806772
leaf_weight=2091.2932232022285 161.98002295196056 160.06459599733353 103.73212864995003
leaf_count=18750 858 1288 429
internal_value=0 -0.808575 0.118445
internal_weight=0 322.045 2195.03
internal_count=21325 2146 19179
shrinkage=0.13089


Tree=14
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=200.529 15.9906 9.82177
threshold=0.0011922500000000004 1.0000000180025095e-35 0.01050795
decision_type=2 2 2
left_child=1 -1 -2
right_child=2 -3 -4
leaf_value=0.015623162767826411 -0.071860699555013161 -0.028666928506826318 -0.13532797756484888
leaf_weight=2063.8248883187771 109.24754539132118 151.16220445930958 144.63925379514694
leaf_count=18750 641 646 1288
internal_value=0 0.0953457 -0.839637
internal_weight=0 2214.99 253.887
internal_count=21325 19396 1929
shrinkage=0.13089


Tree=15
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=165.105 12.1111 10.3021
threshold=0.0011922500000000004 1.0000000180025095e-35 0.01050795
decision_type=2 2 2
left_child=1 -1 -2
right_child=2 -3 -4
leaf_value=0.013753003838970287 -0.065456299523249495 -0.025063894522140792 -0.13135472039219004
leaf_weight=2039.7879648953676 104.84128470718861 150.71637351810932 130.55808162689209
leaf_count=18750 641 646 1288
internal_value=0 0.0837684 -0.794404
internal_weight=0 2190.5 235.399
internal_count=21325 19396 1929
shrinkage=0.13089


Tree=16
num_leaves=5
num_cat=0
split_feature=0 0 0 0
split_gain=136.309 14.5188 4.3726 0.136272
threshold=0.0025796250000000008 1.0000000180025095e-35 0.01050795 9.9591700000000011e-05
decision_type=2 2 2 2
left_child=1 -1 -2 -3
right_child=2 3 -4 -5
leaf_value=0.012089100024626171 -0.066628273910447358 -0.0066241436094967758 -0.12758121029992966 -0.032393382920739143
leaf_weight=2018.7714602798223 64.942745372653008 51.602518200874329 117.78304916620255 134.57430674135685
leaf_count=18750 427 213 1288 647
internal_value=0 0.0665754 -0.827791 -0.202336
internal_weight=0 2204.95 182.726 186.177
internal_count=21325 19610 1715 860
shrinkage=0.13089


Tree=17
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=113.28 11.0849 4.73269
threshold=0.0025796250000000008 1.0000000180025095e-35 0.01050795
decision_type=2 2 2
left_child=1 -1 -2
right_child=2 -3 -4
leaf_value=0.010612638572518077 -0.060863797561325379 -0.023433156188948215 -0.12393517204908359
leaf_weight=2000.4092250019312 62.258247300982475 184.67508155107498 106.2508493065834
leaf_count=18750 427 860 1288
internal_value=0 0.0582137 -0.788413
internal_weight=0 2185.08 168.509
internal_count=21325 19610 1715
shrinkage=0.13089


Tree=18
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=94.083 8.35614 4.97986
threshold=0.0025796250000000008 9.9591700000000011e-05 0.01050795
decision_type=2 2 2
left_child=1 -1 -2
right_child=2 -3 -4
leaf_value=0.0088970647047500583 -0.055465394997646883 -0.025823187431573653 -0.12035807851071022
leaf_weight=2036.2410068660975 59.83521132171154 131.86057107150555 95.883156657218933
leaf_count=18963 427 647 1288
internal_value=0 0.0509179 -0.749659
internal_weight=0 2168.1 155.718
internal_count=21325 19610 1715
shrinkage=0.13089


Tree=19
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=78.0365 6.54972 5.11662
threshold=0.0025796250000000008 1.0000000180025095e-35 0.01050795
decision_type=2 2 2
left_child=1 -1 -2
right_child=2 -3 -4
leaf_value=0.0082019391629894745 -0.050420443501466024 -0.018608996000748765 -0.11680260107141865
leaf_weight=1970.9979649633169 57.657026320695877 182.35884165763855 86.593843996524811
leaf_count=18750 427 860 1288
internal_value=0 0.0444636 -0.711371
internal_weight=0 2153.36 144.251
internal_count=21325 19610 1715
shrinkage=0.13089


Tree=20
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=65.1784 7.9691 1.4246
threshold=0.0053693950000000008 9.9591700000000011e-05 0.04402385000000001
decision_type=2 2 2
left_child=1 -1 -2
right_child=2 -3 -4
leaf_value=0.0068515272398920624 -0.06599848612582912 -0.024486101527776222 -0.12656662398514071
leaf_weight=2010.6606224030256 53.632948517799377 158.2627262622118 51.730003096163273
leaf_count=18963 658 853 851
internal_value=0 0.0339687 -0.765035
internal_weight=0 2168.92 105.363
internal_count=21325 19816 1509
shrinkage=0.13089


Tree=21
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=54.7063 5.95844 1.69458
threshold=0.0053693950000000008 9.9591700000000011e-05 0.04402385000000001
decision_type=2 2 2
left_child=1 -1 -2
right_child=2 -3 -4
leaf_value=0.0059888933412892593 -0.06098200572094635 -0.021587231537420305 -0.12406613692161805
leaf_weight=2000.3821887820959 50.884312272071838 156.53964526951313 46.292951822280884
leaf_count=18963 658 853 851
internal_value=0 0.02958 -0.731232
internal_weight=0 2156.92 97.1773
internal_count=21325 19816 1509
shrinkage=0.13089


Tree=22
num_leaves=5
num_cat=0
split_feature=0 0 0 0
split_gain=45.8612 4.56948 1.91028 0.0182958
threshold=0.0053693950000000008 1.0000000180025095e-35 0.04402385000000001 0.0003976535000000001
decision_type=2 2 2 2
left_child=1 -1 -2 -3
right_child=2 3 -4 -5
leaf_value=0.00556226357586066 -0.056170438185957378 -0.0054932526816375291 -0.12157172679766572 -0.024252226281742208
leaf_weight=1939.5867828279734 48.438122946768999 102.4775842577219 41.448063153773546 104.36961236596107
leaf_count=18750 658 429 851 637
internal_value=0 0.0257276 -0.69746 -0.122882
internal_weight=0 2146.43 89.8862 206.847
internal_count=21325 19816 1509 1066
shrinkage=0.13089


Tree=23
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=38.8902 5.7017
threshold=0.01050795 0.0003976535000000001
decision_type=2 2
left_child=1 -1
right_child=-2 -3
leaf_value=0.0042703630227069423 -0.10308613967071399 -0.025572608040684074
leaf_weight=2033.7917536497116 59.746310696005821 126.60679809749126
leaf_count=19179 1288 858
internal_value=0 0.0183367
internal_weight=0 2160.4
internal_count=21325 20037
shrinkage=0.13089


Tree=24
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=32.8012 4.24835 0.407883
threshold=0.01050795 0.0011922500000000004 1.0000000180025095e-35
decision_type=2 2 2
left_child=1 2 -1
right_child=-2 -3 -4
leaf_value=0.0043143926298171344 -0.099244772804821885 -0.028620488894483563 -0.006929954247553694
leaf_weight=1925.0204320997 54.416697561740875 81.905148580670357 145.22261695563793
leaf_count=18750 1288 641 646
internal_value=0 0.0158739 0.0261591
internal_weight=0 2152.15 2070.24
internal_count=21325 20037 19396
shrinkage=0.13089


Tree=25
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=27.5917 3.17816
threshold=0.01050795 0.0003976535000000001
decision_type=2 2
left_child=1 -1
right_child=-2 -3
leaf_value=0.0032461001584045551 -0.095311779381167028 -0.020260488318784076
leaf_weight=2021.0935054421425 49.692862194031477 122.90897027403116
leaf_count=19179 1288 858
internal_value=0 0.0136217
internal_weight=0 2144
internal_count=21325 20037
shrinkage=0.13089


Tree=26
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=23.1523 2.48747 0.266918
threshold=0.01050795 1.0000000180025095e-35 0.04402385000000001
decision_type=2 2 2
left_child=1 -1 -2
right_child=2 -3 -4
leaf_value=0.0033485408085107741 -0.035900327190669407 -0.012928229706923247 -0.11242294149629284
leaf_weight=1913.8687290251255 17.270923379808664 223.83370804786682 28.238933559507132
leaf_count=18750 437 1287 851
internal_value=0 0.0117566 -0.697458
internal_weight=0 2137.7 45.5099
internal_count=21325 20037 1288
shrinkage=0.13089


Tree=27
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=19.7007 2.0238 0.274785
threshold=0.01050795 0.0011922500000000004 0.04402385000000001
decision_type=2 2 2
left_child=1 -1 -2
right_child=2 -3 -4
leaf_value=0.0023374824330875542 -0.032860062812921413 -0.021866207679274243 -0.10968086716532251
leaf_weight=2053.2934209406376 16.710434205830097 78.447846926748753 25.421625388786197
leaf_count=19396 437 641 851
internal_value=0 0.0101436 -0.669323
internal_weight=0 2131.74 42.1321
internal_count=21325 20037 1288
shrinkage=0.13089


Tree=28
num_leaves=4
num_cat=0
split_feature=0 0 0
split_gain=16.7471 1.41546 0.263438
threshold=0.01050795 0.0011922500000000004 0.04402385000000001
decision_type=2 2 2
left_child=1 -1 -2
right_child=2 -3 -4
leaf_value=0.002035553117815515 -0.030017236684714045 -0.019416014737067059 -0.10685003974833024
leaf_weight=2049.9182535707951 16.211968369781971 77.242374464869499 22.927413567900658
leaf_count=19396 437 641 851
internal_value=0 0.0087131 -0.641105
internal_weight=0 2127.16 39.1394
internal_count=21325 20037 1288
shrinkage=0.13089


Tree=29
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=14.2563 2.26737
threshold=0.04402385000000001 0.0011922500000000004
decision_type=2 2
left_child=1 -1
right_child=-2 -3
leaf_value=0.0017722297397573463 -0.1039270291991902 -0.021777523763976007
leaf_weight=2046.9820948541164 20.720163905993104 91.948586139827967
leaf_count=19396 851 1078
internal_value=0 0.00490076
internal_weight=0 2138.93
internal_count=21325 20474
shrinkage=0.13089


Tree=30
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=12.2308 1.63457
threshold=0.04402385000000001 0.0025796250000000008
decision_type=2 2
left_child=1 -1
right_child=-2 -3
leaf_value=0.0013856204870998597 -0.10091225687432451 -0.023832057314608676
leaf_weight=2073.8016383647919 18.767303956672549 61.078140139579773
leaf_count=19610 851 864
internal_value=0 0.0041432
internal_weight=0 2134.88
internal_count=21325 20474
shrinkage=0.13089


Tree=31
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=10.4765 1.30191
threshold=0.04402385000000001 1.0000000180025095e-35
decision_type=2 2
left_child=1 -1
right_child=-2 -3
leaf_value=0.0019488886929256715 -0.097809763777693653 -0.010724854538772081
leaf_weight=1897.896034643054 17.039601102471352 233.86431689932942
leaf_count=18750 851 1724
internal_value=0 0.00347924
internal_weight=0 2131.76
internal_count=21325 20474
shrinkage=0.13089


Tree=32
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=8.95808 0.912267
threshold=0.04402385000000001 0.0025796250000000008
decision_type=2 2
left_child=1 -1
right_child=-2 -3
leaf_value=0.001091849573601314 -0.094626947531172442 -0.020097532132370106
leaf_weight=2068.5039180070162 15.510904181748627 59.413728449493647
leaf_count=19610 851 864
internal_value=0 0.00292558
internal_weight=0 2127.92
internal_count=21325 20474
shrinkage=0.13089


Tree=33
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=7.64342 0.648756
threshold=0.04402385000000001 1.0000000180025095e-35
decision_type=2 2
left_child=1 -1
right_child=-2 -3
leaf_value=0.0015556050440260511 -0.091374199169083944 -0.0088622516447988547
leaf_weight=1893.4474792331457 14.157899381592868 231.96995370462537
leaf_count=18750 851 1724
internal_value=0 0.00242476
internal_weight=0 2125.42
internal_count=21325 20474
shrinkage=0.13089


Tree=34
num_leaves=2
num_cat=0
split_feature=0
split_gain=6.51045
threshold=0.01050795
decision_type=2
left_child=-1
right_child=-2
leaf_value=0.00051834688991247958 -0.06351612278428391
leaf_weight=2107.7610267102718 27.528344225138426
leaf_count=20037 1288
internal_value=0
internal_weight=0
internal_count=21325
shrinkage=0.13089


Tree=35
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=5.78967 0.247874
threshold=0.04402385000000001 0.0011922500000000004
decision_type=2 2
left_child=1 -1
right_child=-2 -3
leaf_value=0.00091255324634250052 -0.085662958029800423 -0.012702551462742125
leaf_weight=2034.5041081905365 12.185477577149866 86.27077791839838
leaf_count=19396 851 1078
internal_value=0 0.00191893
internal_weight=0 2120.77
internal_count=21325 20474
shrinkage=0.13089


Tree=36
num_leaves=3
num_cat=0
split_feature=0 0
split_gain=4.90337 0.0926784
threshold=0.04402385000000001 1.0000000180025095e-35
decision_type=2 2
left_child=1 -1
right_child=-2 -3
leaf_value=0.0011684427230035121 -0.082290525809403103 -0.0068165486607681125
leaf_weight=1889.0849780291319 11.211826302111147 229.54791922494769
leaf_count=18750 851 1724
internal_value=0 0.00156011
internal_weight=0 2118.63
internal_count=21325 20474
shrinkage=0.13089


Tree=37
num_leaves=2
num_cat=0
split_feature=0
split_gain=4.13796
threshold=0.04402385000000001
decision_type=2
left_child=-1
right_child=-2
leaf_value=0.00016467471883775653 -0.078904217644919214
leaf_weight=2116.3086513914168 10.347959816455839
leaf_count=20474 851
internal_value=0
internal_weight=0
internal_count=21325
shrinkage=0.13089


Tree=38
num_leaves=2
num_cat=0
split_feature=0
split_gain=3.47834
threshold=0.04402385000000001
decision_type=2
left_child=-1
right_child=-2
leaf_value=0.0001431647191983432 -0.0755211975253922
leaf_weight=2116.0833415165544 9.5807343795895559
leaf_count=20474 851
internal_value=0
internal_weight=0
internal_count=21325
shrinkage=0.13089


Tree=39
num_leaves=2
num_cat=0
split_feature=0
split_gain=2.90896
threshold=0.04402385000000001
decision_type=2
left_child=-1
right_child=-2
leaf_value=0.00012446112354654837 -0.072158435436264476
leaf_weight=2115.8876501619816 8.8985683778300864
leaf_count=20474 851
internal_value=0
internal_weight=0
internal_count=21325
shrinkage=0.13089


Tree=40
num_leaves=2
num_cat=0
split_feature=0
split_gain=2.48071
threshold=0.01050795
decision_type=2
left_child=-1
right_child=-2
leaf_value=0.00029115405598437345 -0.046787770887242085
leaf_weight=2102.2301584482193 21.77847608923912
leaf_count=20037 1288
internal_value=0
internal_weight=0
internal_count=21325
shrinkage=0.13089


Tree=41
num_leaves=2
num_cat=0
split_feature=0
split_gain=2.12579
threshold=0.04402385000000001
decision_type=2
left_child=-1
right_child=-2
leaf_value=0.00010858902169945831 -0.066617029967861549
leaf_weight=2114.7367005199194 7.9194040866568676
leaf_count=20474 851
internal_value=0
internal_weight=0
internal_count=21325
shrinkage=0.13089


Tree=42
num_leaves=2
num_cat=0
split_feature=0
split_gain=1.77506
threshold=0.0053693950000000008
decision_type=2
left_child=-1
right_child=-2
leaf_value=0.00042414227351096818 -0.03117067834830867
leaf_weight=2081.9380534440279 40.068188180215657
leaf_count=19816 1509
internal_value=0
internal_weight=0
internal_count=21325
shrinkage=0.13089


Tree=43
num_leaves=2
num_cat=0
split_feature=0
split_gain=1.57647
threshold=0.04402385000000001
decision_type=2
left_child=-1
right_child=-2
leaf_value=0.00010197234197438947 -0.061842009835006914
leaf_weight=2113.1346718631685 7.194230060093104
leaf_count=20474 851
internal_value=0
internal_weight=0
internal_count=21325
shrinkage=0.13089


Tree=44
num_leaves=2
num_cat=0
split_feature=0
split_gain=1.26886
threshold=0.04402385000000001
decision_type=2
left_child=-1
right_child=-2
leaf_value=8.8652322040586707e-05 -0.058723612693134065
leaf_weight=2112.9951883777976 6.7697230624034992
leaf_count=20474 851
internal_value=0
internal_weight=0
internal_count=21325
shrinkage=0.13089


Tree=45
num_leaves=2
num_cat=0
split_feature=0
split_gain=1.0038
threshold=0.04402385000000001
decision_type=2
left_child=-1
right_child=-2
leaf_value=7.7069570483583134e-05 -0.055694387033665901
leaf_weight=2112.8739228565246 6.389471571892499
leaf_count=20474 851
internal_value=0
internal_weight=0
internal_count=21325
shrinkage=0.13089


Tree=46
num_leaves=2
num_cat=0
split_feature=0
split_gain=0.841681
threshold=0.0053693950000000008
decision_type=2
left_child=-1
right_child=-2
leaf_value=0.00033389667956715919 -0.025206151487094128
leaf_weight=2080.9726273268461 37.844197936356068
leaf_count=19816 1509
internal_value=0
internal_weight=0
internal_count=21325
shrinkage=0.13089


Tree=47
num_leaves=2
num_cat=0
split_feature=0
split_gain=0.679399
threshold=0.04402385000000001
decision_type=2
left_child=-1
right_child=-2
leaf_value=7.3101036785521857e-05 -0.051417389826547671
leaf_weight=2111.6197035070509 5.8998517799191168
leaf_count=20474 851
internal_value=0
internal_weight=0
internal_count=21325
shrinkage=0.13089


Tree=48
num_leaves=2
num_cat=0
split_feature=0
split_gain=0.495814
threshold=0.04402385000000001
decision_type=2
left_child=-1
right_child=-2
leaf_value=6.3553295190337268e-05 -0.048639751006127774
leaf_weight=2111.5197754185647 5.6080875270999959
leaf_count=20474 851
internal_value=0
internal_weight=0
internal_count=21325
shrinkage=0.13089


Tree=49
num_leaves=2
num_cat=0
split_feature=0
split_gain=0.349318
threshold=0.0053693950000000008
decision_type=2
left_child=-1
right_child=-2
leaf_value=0.00027248288990187474 -0.021102302606919959
leaf_weight=2080.3162171542645 36.461820169817656
leaf_count=19816 1509
internal_value=0
internal_weight=0
internal_count=21325
shrinkage=0.13089


Tree=50
num_leaves=2
num_cat=0
split_feature=0
split_gain=0.273556
threshold=0.04402385000000001
decision_type=2
left_child=-1
right_child=-2
leaf_value=6.0469665200889638e-05 -0.044803241179460079
leaf_weight=2110.4903337787837 5.2349733952432862
leaf_count=20474 851
internal_value=0
internal_weight=0
internal_count=21325
shrinkage=0.13089


Tree=51
num_leaves=2
num_cat=0
split_feature=0
split_gain=0.145868
threshold=0.04402385000000001
decision_type=2
left_child=-1
right_child=-2
leaf_value=5.2563984547940494e-05 -0.042299375401756018
leaf_weight=2110.4076942373067 5.0083223329856983
leaf_count=20474 851
internal_value=0
internal_weight=0
internal_count=21325
shrinkage=0.13089


Tree=52
num_leaves=2
num_cat=0
split_feature=0
split_gain=0.070346
threshold=0.0025796250000000008
decision_type=2
left_child=-1
right_child=-2
leaf_value=0.0003259007314232656 -0.014373578368592269
leaf_weight=2056.6183921545744 58.520661424379796
leaf_count=19610 1715
internal_value=0
internal_weight=0
internal_count=21325
shrinkage=0.13089


Tree=53
num_leaves=2
num_cat=0
split_feature=0
split_gain=0.00057606
threshold=0.04402385000000001
decision_type=2
left_child=-1
right_child=-2
leaf_value=5.175836497515176e-05 -0.039091491520440565
leaf_weight=2109.2513902299106 4.7354676621034733
leaf_count=20474 851
internal_value=0
internal_weight=0
internal_count=21325
shrinkage=0.13089


end of trees

feature importances:
Gnomad_exomes_ASJ_AF=125

parameters:
[boosting: gbdt]
[objective: binary]
[metric: auc]
[tree_learner: serial]
[device_type: cpu]
[data: ]
[valid: ]
[num_iterations: 100]
[learning_rate: 0.13089]
[num_leaves: 6]
[num_threads: -1]
[max_depth: 32]
[min_data_in_leaf: 6]
[min_sum_hessian_in_leaf: 1.26004]
[bagging_fraction: 0.819095]
[pos_bagging_fraction: 1]
[neg_bagging_fraction: 1]
[bagging_freq: 87]
[bagging_seed: 3]
[feature_fraction: 0.177558]
[feature_fraction_bynode: 1]
[feature_fraction_seed: 2]
[early_stopping_round: 0]
[first_metric_only: 0]
[max_delta_step: 0]
[lambda_l1: 1.4729]
[lambda_l2: 2.92951]
[min_gain_to_split: 0.683451]
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
[max_bin: 13]
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
hCҞ�l���?qI�qJRqKhK hKhG?PbM���hG        hJ����hKhNhG        hG        hG?�      hJ@ hK X   max_binqLKX   min_data_in_leafqMKX   min_sum_hessian_in_leafqNh
hCln4�)�?qO�qPRqQX   bagging_fractionqRh
hCT�@�6�?qS�qTRqUX   bagging_freqqVKWX   feature_fractionqWh
hC��{�8��?qX�qYRqZX	   lambda_l1q[h
hC5r�5���?q\�q]Rq^X	   lambda_l2q_h
hC���ˢo@q`�qaRqbX   min_gain_to_splitqch
hC�-�����?qd�qeRqfX   verboseqgJ����hhX   metricqh]qiX   aucqjauubX   _evals_resultqkNX   _best_scoreqlh5X   _best_iterationqmNX   _other_paramsqn}qo(hLKhMKhNh
hCln4�)�?qp�qqRqrhRh
hCT�@�6�?qs�qtRquhVKWhWh
hC��{�8��?qv�qwRqxh[h
hC5r�5���?qy�qzRq{h_h
hC���ˢo@q|�q}Rq~hch
hC�-�����?q�q�Rq�hhhjuX
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
   _n_classesq�KhLKhMKhNhrhRhuhVKWhWhxh[h{h_h~hch�hhhjX   _leq�csklearn.preprocessing._label
LabelEncoder
q�)�q�}q�(X   classes_q�h�X   _sklearn_versionq�X   0.22.1q�ubX   _fobjq�Nub.