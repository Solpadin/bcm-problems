ID D:\Users\Dima\My_document,FEMAP
SOL SESTATICS
TIME 10000
CEND
  ECHO = NONE
  DISPLACEMENT(PLOT) = ALL
  SPCFORCE(PLOT) = ALL
  OLOAD(PLOT) = ALL
  FORCE(PLOT,CORNER) = ALL
  STRESS(PLOT,CORNER) = ALL
  LOAD = 1
BEGIN BULK
$ ***************************************************************************
$   Written by : FEMAP
$   Version    : 7.00
$   Translator : MSC/NASTRAN
$   From Model : D:\Users\Dima\My_documents\PAVT_2009\models\Pipe\Solid\Pipe_model.MOD
$   Date       : Sat Feb 14 17:01:59 2009
$ ***************************************************************************
$
PARAM,POST,-1
PARAM,OGEOM,NO
PARAM,AUTOSPC,YES
PARAM,GRDPNT,0
CORD2C         1       0      0.      0.      0.      0.      0.      1.+FEMAPC1
+FEMAPC1      1.      0.      1.
CORD2S         2       0      0.      0.      0.      0.      0.      1.+FEMAPC2
+FEMAPC2      1.      0.      1.
$ FEMAP Load Set 1 : Source Displn -.001
PLOAD4         1     151      1.                             193     204
PLOAD4         1     157      1.                             194     208
PLOAD4         1     163      1.                             141     212
PLOAD4         1     169      1.                             195     216
PLOAD4         1     175      1.                             196     226
PLOAD4         1     181      1.                             251     203
PLOAD4         1     187      1.                             204     207
PLOAD4         1     193      1.                             208     211
PLOAD4         1     199      1.                             212     215
PLOAD4         1     205      1.                             216     197
PLOAD4         1     211      1.                             192     202
PLOAD4         1     217      1.                             203     206
PLOAD4         1     223      1.                             207     210
PLOAD4         1     229      1.                             211     214
PLOAD4         1     235      1.                             215     225
PLOAD4         1     241      1.                             252     201
PLOAD4         1     247      1.                             202     205
PLOAD4         1     253      1.                             206     209
PLOAD4         1     259      1.                             210     213
PLOAD4         1     265      1.                             214     198
PLOAD4         1     271      1.                             191     200
PLOAD4         1     277      1.                             201     168
PLOAD4         1     283      1.                             205     169
PLOAD4         1     289      1.                             209     170
PLOAD4         1     295      1.                             213     199
$ FEMAP Load Set 2 : Absorbing Uptake 1.0
PLOAD4         2     476      2.                             531     436
PLOAD4         2     477      2.                             475     432
PLOAD4         2     478      2.                             423     428
PLOAD4         2     479      2.                             474     424
PLOAD4         2     480      2.                             473     507
PLOAD4         2     506      2.                             422     437
PLOAD4         2     507      2.                             436     433
PLOAD4         2     508      2.                             432     429
PLOAD4         2     509      2.                             428     425
PLOAD4         2     510      2.                             424     414
PLOAD4         2     536      2.                             532     438
PLOAD4         2     537      2.                             437     434
PLOAD4         2     538      2.                             433     430
PLOAD4         2     539      2.                             429     426
PLOAD4         2     540      2.                             425     415
PLOAD4         2     566      2.                             421     439
PLOAD4         2     567      2.                             438     435
PLOAD4         2     568      2.                             434     431
PLOAD4         2     569      2.                             430     427
PLOAD4         2     570      2.                             426     416
PLOAD4         2     596      2.                             420     446
PLOAD4         2     597      2.                             439     447
PLOAD4         2     598      2.                             435     418
PLOAD4         2     599      2.                             431     448
PLOAD4         2     600      2.                             427     417
$ FEMAP Load Set 3 : Boundary Skews
PLOAD4         3     151      3.                             193     142
PLOAD4         3     152      3.                             133     143
PLOAD4         3     153      3.                             134     144
PLOAD4         3     154      3.                             135     145
PLOAD4         3     155      3.                             136     146
PLOAD4         3     156      3.                             137     590
PLOAD4         3     157      3.                             194     147
PLOAD4         3     158      3.                             142     148
PLOAD4         3     159      3.                             143     149
PLOAD4         3     160      3.                             144     150
PLOAD4         3     161      3.                             145     151
PLOAD4         3     162      3.                             146     276
PLOAD4         3     163      3.                             141     152
PLOAD4         3     164      3.                             147     153
PLOAD4         3     165      3.                             148     154
PLOAD4         3     166      3.                             149     155
PLOAD4         3     167      3.                             150     156
PLOAD4         3     168      3.                             151     591
PLOAD4         3     169      3.                             195     157
PLOAD4         3     170      3.                             152     158
PLOAD4         3     171      3.                             153     159
PLOAD4         3     172      3.                             154     160
PLOAD4         3     173      3.                             155     161
PLOAD4         3     174      3.                             156     592
PLOAD4         3     175      3.                             196     218
PLOAD4         3     176      3.                             157     140
PLOAD4         3     177      3.                             158     139
PLOAD4         3     178      3.                             159     138
PLOAD4         3     179      3.                             160     219
PLOAD4         3     180      3.                             161     275
PLOAD4         3       1      3.                               1       9
PLOAD4         3       2      3.                               2      10
PLOAD4         3       3      3.                              80      11
PLOAD4         3       4      3.                               3      12
PLOAD4         3       5      3.                               4     123
PLOAD4         3       6      3.                              28      13
PLOAD4         3       7      3.                               9      14
PLOAD4         3       8      3.                              10      15
PLOAD4         3       9      3.                              11      16
PLOAD4         3      10      3.                              12     122
PLOAD4         3      11      3.                              29      17
PLOAD4         3      12      3.                              13      18
PLOAD4         3      13      3.                              14      19
PLOAD4         3      14      3.                              15      20
PLOAD4         3      15      3.                              16     684
PLOAD4         3      16      3.                               8      21
PLOAD4         3      17      3.                              17      22
PLOAD4         3      18      3.                              18      23
PLOAD4         3      19      3.                              19      24
PLOAD4         3      20      3.                              20       5
PLOAD4         3      21      3.                               7      57
PLOAD4         3      22      3.                              21       6
PLOAD4         3      23      3.                              22      58
PLOAD4         3      24      3.                              23      59
PLOAD4         3      25      3.                              24     613
PLOAD4         3       1      3.                               1      45
PLOAD4         3       6      3.                              28      40
PLOAD4         3      11      3.                              29      35
PLOAD4         3      16      3.                               8      30
PLOAD4         3      21      3.                               7      25
PLOAD4         3      26      3.                              81      46
PLOAD4         3      31      3.                              45      41
PLOAD4         3      36      3.                              40      36
PLOAD4         3      41      3.                              35      31
PLOAD4         3      46      3.                              30      55
PLOAD4         3      51      3.                              82      47
PLOAD4         3      56      3.                              46      42
PLOAD4         3      61      3.                              41      37
PLOAD4         3      66      3.                              36      32
PLOAD4         3      71      3.                              31      54
PLOAD4         3      76      3.                              83      48
PLOAD4         3      81      3.                              47      43
PLOAD4         3      86      3.                              42      38
PLOAD4         3      91      3.                              37      33
PLOAD4         3      96      3.                              32      53
PLOAD4         3     101      3.                              84      49
PLOAD4         3     106      3.                              48      44
PLOAD4         3     111      3.                              43      39
PLOAD4         3     116      3.                              38      34
PLOAD4         3     121      3.                              33      26
PLOAD4         3     126      3.                              85     109
PLOAD4         3     131      3.                              49     557
PLOAD4         3     136      3.                              44      27
PLOAD4         3     141      3.                              39     470
PLOAD4         3     146      3.                              34     110
PLOAD4         3     571      3.                             554     468
PLOAD4         3     572      3.                             121     463
PLOAD4         3     573      3.                             571     458
PLOAD4         3     574      3.                             112     453
PLOAD4         3     575      3.                             555     441
PLOAD4         3     576      3.                             442     467
PLOAD4         3     577      3.                             468     462
PLOAD4         3     578      3.                             463     457
PLOAD4         3     579      3.                             458     452
PLOAD4         3     580      3.                             453     504
PLOAD4         3     581      3.                             443     466
PLOAD4         3     582      3.                             467     461
PLOAD4         3     583      3.                             462     456
PLOAD4         3     584      3.                             457     451
PLOAD4         3     585      3.                             452     440
PLOAD4         3     586      3.                             444     465
PLOAD4         3     587      3.                             466     460
PLOAD4         3     588      3.                             461     455
PLOAD4         3     589      3.                             456     450
PLOAD4         3     590      3.                             451     505
PLOAD4         3     591      3.                             445     464
PLOAD4         3     592      3.                             465     459
PLOAD4         3     593      3.                             460     454
PLOAD4         3     594      3.                             455     449
PLOAD4         3     595      3.                             450     506
PLOAD4         3     596      3.                             533     446
PLOAD4         3     597      3.                             464     447
PLOAD4         3     598      3.                             459     418
PLOAD4         3     599      3.                             454     448
PLOAD4         3     600      3.                             449     417
PLOAD4         3     601      3.                             554     584
PLOAD4         3     607      3.                             121     580
PLOAD4         3     613      3.                             571     576
PLOAD4         3     619      3.                             112     572
PLOAD4         3     625      3.                             555     629
PLOAD4         3     631      3.                             570     585
PLOAD4         3     637      3.                             584     581
PLOAD4         3     643      3.                             580     577
PLOAD4         3     649      3.                             576     573
PLOAD4         3     655      3.                             572     628
PLOAD4         3     661      3.                             659     586
PLOAD4         3     667      3.                             585     582
PLOAD4         3     673      3.                             581     578
PLOAD4         3     679      3.                             577     574
PLOAD4         3     685      3.                             573     627
PLOAD4         3     691      3.                             569     587
PLOAD4         3     697      3.                             586     583
PLOAD4         3     703      3.                             582     579
PLOAD4         3     709      3.                             578     575
PLOAD4         3     715      3.                             574     626
PLOAD4         3     721      3.                             660     568
PLOAD4         3     727      3.                             587     596
PLOAD4         3     733      3.                             583     597
PLOAD4         3     739      3.                             579     332
PLOAD4         3     745      3.                             575     567
PLOAD4         3     271      3.                             253     186
PLOAD4         3     272      3.                             167     187
PLOAD4         3     273      3.                             166     188
PLOAD4         3     274      3.                             165     189
PLOAD4         3     275      3.                             164     190
PLOAD4         3     276      3.                             254     278
PLOAD4         3     277      3.                             200     181
PLOAD4         3     278      3.                             186     182
PLOAD4         3     279      3.                             187     183
PLOAD4         3     280      3.                             188     184
PLOAD4         3     281      3.                             189     185
PLOAD4         3     282      3.                             190     290
PLOAD4         3     283      3.                             168     176
PLOAD4         3     284      3.                             181     177
PLOAD4         3     285      3.                             182     178
PLOAD4         3     286      3.                             183     179
PLOAD4         3     287      3.                             184     180
PLOAD4         3     288      3.                             185     291
PLOAD4         3     289      3.                             169     171
PLOAD4         3     290      3.                             176     172
PLOAD4         3     291      3.                             177     173
PLOAD4         3     292      3.                             178     174
PLOAD4         3     293      3.                             179     175
PLOAD4         3     294      3.                             180     292
PLOAD4         3     295      3.                             170     224
PLOAD4         3     296      3.                             171     223
PLOAD4         3     297      3.                             172     222
PLOAD4         3     298      3.                             173     162
PLOAD4         3     299      3.                             174     163
PLOAD4         3     300      3.                             175     293
PLOAD4         3     175      3.                             217     227
PLOAD4         3     176      3.                             218     228
PLOAD4         3     177      3.                             140     229
PLOAD4         3     178      3.                             139     230
PLOAD4         3     179      3.                             138     231
PLOAD4         3     180      3.                             219     345
PLOAD4         3     205      3.                             226     232
PLOAD4         3     206      3.                             227     233
PLOAD4         3     207      3.                             228     234
PLOAD4         3     208      3.                             229     235
PLOAD4         3     209      3.                             230     236
PLOAD4         3     210      3.                             231     279
PLOAD4         3     235      3.                             197     237
PLOAD4         3     236      3.                             232     238
PLOAD4         3     237      3.                             233     239
PLOAD4         3     238      3.                             234     240
PLOAD4         3     239      3.                             235     241
PLOAD4         3     240      3.                             236     220
PLOAD4         3     265      3.                             225     242
PLOAD4         3     266      3.                             237     243
PLOAD4         3     267      3.                             238     244
PLOAD4         3     268      3.                             239     245
PLOAD4         3     269      3.                             240     246
PLOAD4         3     270      3.                             241     221
PLOAD4         3     295      3.                             198     224
PLOAD4         3     296      3.                             242     223
PLOAD4         3     297      3.                             243     222
PLOAD4         3     298      3.                             244     162
PLOAD4         3     299      3.                             245     163
PLOAD4         3     300      3.                             246     293
PLOAD4         3     151      3.                             193     274
PLOAD4         3     152      3.                             133     270
PLOAD4         3     153      3.                             134     266
PLOAD4         3     154      3.                             135     262
PLOAD4         3     155      3.                             136     258
PLOAD4         3     156      3.                             137     277
PLOAD4         3     181      3.                             251     273
PLOAD4         3     182      3.                             274     269
PLOAD4         3     183      3.                             270     265
PLOAD4         3     184      3.                             266     261
PLOAD4         3     185      3.                             262     257
PLOAD4         3     186      3.                             258     250
PLOAD4         3     211      3.                             192     272
PLOAD4         3     212      3.                             273     268
PLOAD4         3     213      3.                             269     264
PLOAD4         3     214      3.                             265     260
PLOAD4         3     215      3.                             261     256
PLOAD4         3     216      3.                             257     249
PLOAD4         3     241      3.                             252     271
PLOAD4         3     242      3.                             272     267
PLOAD4         3     243      3.                             268     263
PLOAD4         3     244      3.                             264     259
PLOAD4         3     245      3.                             260     255
PLOAD4         3     246      3.                             256     248
PLOAD4         3     271      3.                             191     167
PLOAD4         3     272      3.                             271     166
PLOAD4         3     273      3.                             267     165
PLOAD4         3     274      3.                             263     164
PLOAD4         3     275      3.                             259     254
PLOAD4         3     276      3.                             255     247
PLOAD4         3      21      3.                              56      75
PLOAD4         3      22      3.                              57      70
PLOAD4         3      23      3.                               6      65
PLOAD4         3      24      3.                              58      60
PLOAD4         3      25      3.                              59     631
PLOAD4         3      46      3.                              25      76
PLOAD4         3      47      3.                              75      71
PLOAD4         3      48      3.                              70      66
PLOAD4         3      49      3.                              65      61
PLOAD4         3      50      3.                              60      50
PLOAD4         3      71      3.                              55      77
PLOAD4         3      72      3.                              76      72
PLOAD4         3      73      3.                              71      67
PLOAD4         3      74      3.                              66      62
PLOAD4         3      75      3.                              61     614
PLOAD4         3      96      3.                              54      78
PLOAD4         3      97      3.                              77      73
PLOAD4         3      98      3.                              72      68
PLOAD4         3      99      3.                              67      63
PLOAD4         3     100      3.                              62     630
PLOAD4         3     121      3.                              53      79
PLOAD4         3     122      3.                              78      74
PLOAD4         3     123      3.                              73      69
PLOAD4         3     124      3.                              68      64
PLOAD4         3     125      3.                              63      51
PLOAD4         3     146      3.                              26      52
PLOAD4         3     147      3.                              79     502
PLOAD4         3     148      3.                              74     556
PLOAD4         3     149      3.                              69     111
PLOAD4         3     150      3.                              64     503
PLOAD4         3       1      3.                               1     102
PLOAD4         3       2      3.                               2     103
PLOAD4         3       3      3.                              80     104
PLOAD4         3       4      3.                               3     105
PLOAD4         3       5      3.                               4     120
PLOAD4         3      26      3.                              81      98
PLOAD4         3      27      3.                             102      99
PLOAD4         3      28      3.                             103     100
PLOAD4         3      29      3.                             104     101
PLOAD4         3      30      3.                             105     656
PLOAD4         3      51      3.                              82      94
PLOAD4         3      52      3.                              98      95
PLOAD4         3      53      3.                              99      96
PLOAD4         3      54      3.                             100      97
PLOAD4         3      55      3.                             101     657
PLOAD4         3      76      3.                              83      90
PLOAD4         3      77      3.                              94      91
PLOAD4         3      78      3.                              95      92
PLOAD4         3      79      3.                              96      93
PLOAD4         3      80      3.                              97     611
PLOAD4         3     101      3.                              84      86
PLOAD4         3     102      3.                              90      87
PLOAD4         3     103      3.                              91      88
PLOAD4         3     104      3.                              92      89
PLOAD4         3     105      3.                              93     658
PLOAD4         3     126      3.                              85     108
PLOAD4         3     127      3.                              86     107
PLOAD4         3     128      3.                              87     106
PLOAD4         3     129      3.                              88     528
PLOAD4         3     130      3.                              89     554
PLOAD4         3     421      3.                             247     304
PLOAD4         3     422      3.                             369     308
PLOAD4         3     423      3.                             370     312
PLOAD4         3     424      3.                             371     316
PLOAD4         3     425      3.                             303     320
PLOAD4         3     426      3.                             302     301
PLOAD4         3     427      3.                             278     305
PLOAD4         3     428      3.                             304     309
PLOAD4         3     429      3.                             308     313
PLOAD4         3     430      3.                             312     317
PLOAD4         3     431      3.                             316     321
PLOAD4         3     432      3.                             320     300
PLOAD4         3     433      3.                             290     306
PLOAD4         3     434      3.                             305     310
PLOAD4         3     435      3.                             309     314
PLOAD4         3     436      3.                             313     318
PLOAD4         3     437      3.                             317     322
PLOAD4         3     438      3.                             321     299
PLOAD4         3     439      3.                             291     307
PLOAD4         3     440      3.                             306     311
PLOAD4         3     441      3.                             310     315
PLOAD4         3     442      3.                             314     319
PLOAD4         3     443      3.                             318     323
PLOAD4         3     444      3.                             322     395
PLOAD4         3     445      3.                             292     344
PLOAD4         3     446      3.                             307     294
PLOAD4         3     447      3.                             311     295
PLOAD4         3     448      3.                             315     296
PLOAD4         3     449      3.                             319     297
PLOAD4         3     450      3.                             323     298
PLOAD4         3     325      3.                             275     346
PLOAD4         3     326      3.                             333     347
PLOAD4         3     327      3.                             593     348
PLOAD4         3     328      3.                             594     349
PLOAD4         3     329      3.                             662     350
PLOAD4         3     330      3.                             595     396
PLOAD4         3     355      3.                             345     351
PLOAD4         3     356      3.                             346     352
PLOAD4         3     357      3.                             347     353
PLOAD4         3     358      3.                             348     354
PLOAD4         3     359      3.                             349     355
PLOAD4         3     360      3.                             350     341
PLOAD4         3     385      3.                             279     356
PLOAD4         3     386      3.                             351     357
PLOAD4         3     387      3.                             352     358
PLOAD4         3     388      3.                             353     359
PLOAD4         3     389      3.                             354     360
PLOAD4         3     390      3.                             355     342
PLOAD4         3     415      3.                             220     361
PLOAD4         3     416      3.                             356     362
PLOAD4         3     417      3.                             357     363
PLOAD4         3     418      3.                             358     364
PLOAD4         3     419      3.                             359     365
PLOAD4         3     420      3.                             360     343
PLOAD4         3     445      3.                             221     344
PLOAD4         3     446      3.                             361     294
PLOAD4         3     447      3.                             362     295
PLOAD4         3     448      3.                             363     296
PLOAD4         3     449      3.                             364     297
PLOAD4         3     450      3.                             365     298
PLOAD4         3     301      3.                             589     391
PLOAD4         3     302      3.                             330     387
PLOAD4         3     303      3.                             331     383
PLOAD4         3     304      3.                             368     379
PLOAD4         3     305      3.                             635     375
PLOAD4         3     306      3.                             588     392
PLOAD4         3     331      3.                             277     390
PLOAD4         3     332      3.                             391     386
PLOAD4         3     333      3.                             387     382
PLOAD4         3     334      3.                             383     378
PLOAD4         3     335      3.                             379     374
PLOAD4         3     336      3.                             375     367
PLOAD4         3     361      3.                             250     389
PLOAD4         3     362      3.                             390     385
PLOAD4         3     363      3.                             386     381
PLOAD4         3     364      3.                             382     377
PLOAD4         3     365      3.                             378     373
PLOAD4         3     366      3.                             374     393
PLOAD4         3     391      3.                             249     388
PLOAD4         3     392      3.                             389     384
PLOAD4         3     393      3.                             385     380
PLOAD4         3     394      3.                             381     376
PLOAD4         3     395      3.                             377     372
PLOAD4         3     396      3.                             373     394
PLOAD4         3     421      3.                             248     369
PLOAD4         3     422      3.                             388     370
PLOAD4         3     423      3.                             384     371
PLOAD4         3     424      3.                             380     303
PLOAD4         3     425      3.                             376     302
PLOAD4         3     426      3.                             372     366
PLOAD4         3     306      3.                             332     392
PLOAD4         3     312      3.                             597     400
PLOAD4         3     318      3.                             596     399
PLOAD4         3     324      3.                             568     398
PLOAD4         3     330      3.                             661     397
PLOAD4         3     336      3.                             400     367
PLOAD4         3     342      3.                             399     404
PLOAD4         3     348      3.                             398     403
PLOAD4         3     354      3.                             397     402
PLOAD4         3     360      3.                             396     401
PLOAD4         3     366      3.                             404     393
PLOAD4         3     372      3.                             403     408
PLOAD4         3     378      3.                             402     407
PLOAD4         3     384      3.                             401     406
PLOAD4         3     390      3.                             341     405
PLOAD4         3     396      3.                             408     394
PLOAD4         3     402      3.                             407     412
PLOAD4         3     408      3.                             406     411
PLOAD4         3     414      3.                             405     410
PLOAD4         3     420      3.                             342     409
PLOAD4         3     426      3.                             412     366
PLOAD4         3     432      3.                             411     301
PLOAD4         3     438      3.                             410     300
PLOAD4         3     444      3.                             409     299
PLOAD4         3     450      3.                             343     395
PLOAD4         3     451      3.                             469     479
PLOAD4         3     452      3.                             109     480
PLOAD4         3     453      3.                             557     481
PLOAD4         3     454      3.                              27     482
PLOAD4         3     455      3.                             470     471
PLOAD4         3     456      3.                             529     483
PLOAD4         3     457      3.                             479     484
PLOAD4         3     458      3.                             480     485
PLOAD4         3     459      3.                             481     486
PLOAD4         3     460      3.                             482     472
PLOAD4         3     461      3.                             478     487
PLOAD4         3     462      3.                             483     488
PLOAD4         3     463      3.                             484     489
PLOAD4         3     464      3.                             485     490
PLOAD4         3     465      3.                             486     501
PLOAD4         3     466      3.                             530     491
PLOAD4         3     467      3.                             487     492
PLOAD4         3     468      3.                             488     493
PLOAD4         3     469      3.                             489     494
PLOAD4         3     470      3.                             490     500
PLOAD4         3     471      3.                             477     495
PLOAD4         3     472      3.                             491     496
PLOAD4         3     473      3.                             492     497
PLOAD4         3     474      3.                             493     498
PLOAD4         3     475      3.                             494     499
PLOAD4         3     476      3.                             476     475
PLOAD4         3     477      3.                             495     423
PLOAD4         3     478      3.                             496     474
PLOAD4         3     479      3.                             497     473
PLOAD4         3     480      3.                             498     413
PLOAD4         3     455      3.                             471      52
PLOAD4         3     460      3.                             472     512
PLOAD4         3     465      3.                             501     511
PLOAD4         3     470      3.                             500     510
PLOAD4         3     475      3.                             499     509
PLOAD4         3     480      3.                             413     508
PLOAD4         3     485      3.                             512     502
PLOAD4         3     490      3.                             511     517
PLOAD4         3     495      3.                             510     516
PLOAD4         3     500      3.                             509     515
PLOAD4         3     505      3.                             508     514
PLOAD4         3     510      3.                             507     513
PLOAD4         3     515      3.                             517     556
PLOAD4         3     520      3.                             516     522
PLOAD4         3     525      3.                             515     521
PLOAD4         3     530      3.                             514     520
PLOAD4         3     535      3.                             513     519
PLOAD4         3     540      3.                             414     518
PLOAD4         3     545      3.                             522     111
PLOAD4         3     550      3.                             521     527
PLOAD4         3     555      3.                             520     526
PLOAD4         3     560      3.                             519     525
PLOAD4         3     565      3.                             518     524
PLOAD4         3     570      3.                             415     523
PLOAD4         3     575      3.                             527     503
PLOAD4         3     580      3.                             526     441
PLOAD4         3     585      3.                             525     504
PLOAD4         3     590      3.                             524     440
PLOAD4         3     595      3.                             523     505
PLOAD4         3     600      3.                             416     506
PLOAD4         3     451      3.                             469     537
PLOAD4         3     456      3.                             529     541
PLOAD4         3     461      3.                             478     545
PLOAD4         3     466      3.                             530     549
PLOAD4         3     471      3.                             477     553
PLOAD4         3     476      3.                             476     422
PLOAD4         3     481      3.                             108     536
PLOAD4         3     486      3.                             537     540
PLOAD4         3     491      3.                             541     544
PLOAD4         3     496      3.                             545     548
PLOAD4         3     501      3.                             549     552
PLOAD4         3     506      3.                             553     532
PLOAD4         3     511      3.                             107     535
PLOAD4         3     516      3.                             536     539
PLOAD4         3     521      3.                             540     543
PLOAD4         3     526      3.                             544     547
PLOAD4         3     531      3.                             548     551
PLOAD4         3     536      3.                             552     421
PLOAD4         3     541      3.                             106     534
PLOAD4         3     546      3.                             535     538
PLOAD4         3     551      3.                             539     542
PLOAD4         3     556      3.                             543     546
PLOAD4         3     561      3.                             547     550
PLOAD4         3     566      3.                             551     420
PLOAD4         3     571      3.                             528     442
PLOAD4         3     576      3.                             534     443
PLOAD4         3     581      3.                             538     444
PLOAD4         3     586      3.                             542     445
PLOAD4         3     591      3.                             546     533
PLOAD4         3     596      3.                             550     419
PLOAD4         3     625      3.                             503     639
PLOAD4         3     626      3.                              51     643
PLOAD4         3     627      3.                             630     647
PLOAD4         3     628      3.                             614     651
PLOAD4         3     629      3.                              50     655
PLOAD4         3     630      3.                             631     632
PLOAD4         3     655      3.                             629     638
PLOAD4         3     656      3.                             639     642
PLOAD4         3     657      3.                             643     646
PLOAD4         3     658      3.                             647     650
PLOAD4         3     659      3.                             651     654
PLOAD4         3     660      3.                             655     633
PLOAD4         3     685      3.                             628     637
PLOAD4         3     686      3.                             638     641
PLOAD4         3     687      3.                             642     645
PLOAD4         3     688      3.                             646     649
PLOAD4         3     689      3.                             650     653
PLOAD4         3     690      3.                             654     634
PLOAD4         3     715      3.                             627     636
PLOAD4         3     716      3.                             637     640
PLOAD4         3     717      3.                             641     644
PLOAD4         3     718      3.                             645     648
PLOAD4         3     719      3.                             649     652
PLOAD4         3     720      3.                             653     688
PLOAD4         3     745      3.                             626     588
PLOAD4         3     746      3.                             636     635
PLOAD4         3     747      3.                             640     368
PLOAD4         3     748      3.                             644     331
PLOAD4         3     749      3.                             648     330
PLOAD4         3     750      3.                             652     589
PLOAD4         3     601      3.                             554     668
PLOAD4         3     602      3.                             658     667
PLOAD4         3     603      3.                             611     666
PLOAD4         3     604      3.                             657     665
PLOAD4         3     605      3.                             656     664
PLOAD4         3     606      3.                             120     685
PLOAD4         3     631      3.                             570     673
PLOAD4         3     632      3.                             668     672
PLOAD4         3     633      3.                             667     671
PLOAD4         3     634      3.                             666     670
PLOAD4         3     635      3.                             665     669
PLOAD4         3     636      3.                             664     663
PLOAD4         3     661      3.                             659     678
PLOAD4         3     662      3.                             673     677
PLOAD4         3     663      3.                             672     676
PLOAD4         3     664      3.                             671     675
PLOAD4         3     665      3.                             670     674
PLOAD4         3     666      3.                             669     686
PLOAD4         3     691      3.                             569     683
PLOAD4         3     692      3.                             678     682
PLOAD4         3     693      3.                             677     681
PLOAD4         3     694      3.                             676     680
PLOAD4         3     695      3.                             675     679
PLOAD4         3     696      3.                             674     687
PLOAD4         3     721      3.                             660     595
PLOAD4         3     722      3.                             683     662
PLOAD4         3     723      3.                             682     594
PLOAD4         3     724      3.                             681     593
PLOAD4         3     725      3.                             680     333
PLOAD4         3     726      3.                             679     275
PLOAD4         3     606      3.                             123     685
PLOAD4         3     612      3.                             122     692
PLOAD4         3     618      3.                             684     691
PLOAD4         3     624      3.                               5     690
PLOAD4         3     630      3.                             613     689
PLOAD4         3     636      3.                             692     663
PLOAD4         3     642      3.                             691     696
PLOAD4         3     648      3.                             690     695
PLOAD4         3     654      3.                             689     694
PLOAD4         3     660      3.                             632     693
PLOAD4         3     666      3.                             696     686
PLOAD4         3     672      3.                             695     700
PLOAD4         3     678      3.                             694     699
PLOAD4         3     684      3.                             693     698
PLOAD4         3     690      3.                             633     697
PLOAD4         3     696      3.                             700     687
PLOAD4         3     702      3.                             699     704
PLOAD4         3     708      3.                             698     703
PLOAD4         3     714      3.                             697     702
PLOAD4         3     720      3.                             634     701
PLOAD4         3     726      3.                             704     275
PLOAD4         3     732      3.                             703     592
PLOAD4         3     738      3.                             702     591
PLOAD4         3     744      3.                             701     276
PLOAD4         3     750      3.                             688     590
$ FEMAP Property 1 : Acoustic space
PSOLID         1       1       0        
$ FEMAP Property 2 : Rigid body
PSOLID         2       1       0        
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           1       0     1.2     -2.     0.5       0        
GRID           2       0     1.2    -1.8     0.5       0        
GRID           3       0     1.2    -1.4     0.5       0        
GRID           4       0     1.2    -1.2     0.5       0        
GRID           5       0     1.2     -1.    -0.3       0        
GRID           6       0     1.2    -1.6    -0.5       0        
GRID           7       0     1.2     -2.    -0.3       0        
GRID           8       0     1.2     -2.    -0.1       0        
GRID           9       0     1.2    -1.8     0.3       0        
GRID          10       0     1.2    -1.6     0.3       0        
GRID          11       0     1.2    -1.4     0.3       0        
GRID          12       0     1.2    -1.2     0.3       0        
GRID          13       0     1.2    -1.8     0.1       0        
GRID          14       0     1.2    -1.6     0.1       0        
GRID          15       0     1.2    -1.4     0.1       0        
GRID          16       0     1.2    -1.2     0.1       0        
GRID          17       0     1.2    -1.8    -0.1       0        
GRID          18       0     1.2    -1.6    -0.1       0        
GRID          19       0     1.2    -1.4    -0.1       0        
GRID          20       0     1.2    -1.2    -0.1       0        
GRID          21       0     1.2    -1.8    -0.3       0        
GRID          22       0     1.2    -1.6    -0.3       0        
GRID          23       0     1.2    -1.4    -0.3       0        
GRID          24       0     1.2    -1.2    -0.3       0        
GRID          25       0     1.4     -2.    -0.5       0        
GRID          26       0     2.2     -2.    -0.5       0        
GRID          27       0     2.4     -2.    -0.1       0        
GRID          28       0     1.2     -2.     0.3       0        
GRID          29       0     1.2     -2.     0.1       0        
GRID          30       0     1.4     -2.    -0.3       0        
GRID          31       0     1.6     -2.    -0.3       0        
GRID          32       0     1.8     -2.    -0.3       0        
GRID          33       0      2.     -2.    -0.3       0        
GRID          34       0     2.2     -2.    -0.3       0        
GRID          35       0     1.4     -2.    -0.1       0        
GRID          36       0     1.6     -2.    -0.1       0        
GRID          37       0     1.8     -2.    -0.1       0        
GRID          38       0      2.     -2.    -0.1       0        
GRID          39       0     2.2     -2.    -0.1       0        
GRID          40       0     1.4     -2.     0.1       0        
GRID          41       0     1.6     -2.     0.1       0        
GRID          42       0     1.8     -2.     0.1       0        
GRID          43       0      2.     -2.     0.1       0        
GRID          44       0     2.2     -2.     0.1       0        
GRID          45       0     1.4     -2.     0.3       0        
GRID          46       0     1.6     -2.     0.3       0        
GRID          47       0     1.8     -2.     0.3       0        
GRID          48       0      2.     -2.     0.3       0        
GRID          49       0     2.2     -2.     0.3       0        
GRID          50       0     1.6     -1.    -0.5       0        
GRID          51       0     2.2     -1.    -0.5       0        
GRID          52       0     2.4    -1.8    -0.5       0        
GRID          53       0      2.     -2.    -0.5       0        
GRID          54       0     1.8     -2.    -0.5       0        
GRID          55       0     1.6     -2.    -0.5       0        
GRID          56       0     1.2     -2.    -0.5       0        
GRID          57       0     1.2    -1.8    -0.5       0        
GRID          58       0     1.2    -1.4    -0.5       0        
GRID          59       0     1.2    -1.2    -0.5       0        
GRID          60       0     1.4    -1.2    -0.5       0        
GRID          61       0     1.6    -1.2    -0.5       0        
GRID          62       0     1.8    -1.2    -0.5       0        
GRID          63       0      2.    -1.2    -0.5       0        
GRID          64       0     2.2    -1.2    -0.5       0        
GRID          65       0     1.4    -1.4    -0.5       0        
GRID          66       0     1.6    -1.4    -0.5       0        
GRID          67       0     1.8    -1.4    -0.5       0        
GRID          68       0      2.    -1.4    -0.5       0        
GRID          69       0     2.2    -1.4    -0.5       0        
GRID          70       0     1.4    -1.6    -0.5       0        
GRID          71       0     1.6    -1.6    -0.5       0        
GRID          72       0     1.8    -1.6    -0.5       0        
GRID          73       0      2.    -1.6    -0.5       0        
GRID          74       0     2.2    -1.6    -0.5       0        
GRID          75       0     1.4    -1.8    -0.5       0        
GRID          76       0     1.6    -1.8    -0.5       0        
GRID          77       0     1.8    -1.8    -0.5       0        
GRID          78       0      2.    -1.8    -0.5       0        
GRID          79       0     2.2    -1.8    -0.5       0        
GRID          80       0     1.2    -1.6     0.5       0        
GRID          81       0     1.4     -2.     0.5       0        
GRID          82       0     1.6     -2.     0.5       0        
GRID          83       0     1.8     -2.     0.5       0        
GRID          84       0      2.     -2.     0.5       0        
GRID          85       0     2.2     -2.     0.5       0        
GRID          86       0     2.2    -1.8     0.5       0        
GRID          87       0     2.2    -1.6     0.5       0        
GRID          88       0     2.2    -1.4     0.5       0        
GRID          89       0     2.2    -1.2     0.5       0        
GRID          90       0      2.    -1.8     0.5       0        
GRID          91       0      2.    -1.6     0.5       0        
GRID          92       0      2.    -1.4     0.5       0        
GRID          93       0      2.    -1.2     0.5       0        
GRID          94       0     1.8    -1.8     0.5       0        
GRID          95       0     1.8    -1.6     0.5       0        
GRID          96       0     1.8    -1.4     0.5       0        
GRID          97       0     1.8    -1.2     0.5       0        
GRID          98       0     1.6    -1.8     0.5       0        
GRID          99       0     1.6    -1.6     0.5       0        
GRID         100       0     1.6    -1.4     0.5       0        
GRID         101       0     1.6    -1.2     0.5       0        
GRID         102       0     1.4    -1.8     0.5       0        
GRID         103       0     1.4    -1.6     0.5       0        
GRID         104       0     1.4    -1.4     0.5       0        
GRID         105       0     1.4    -1.2     0.5       0        
GRID         106       0     2.4    -1.4     0.5       0        
GRID         107       0     2.4    -1.6     0.5       0        
GRID         108       0     2.4    -1.8     0.5       0        
GRID         109       0     2.4     -2.     0.3       0        
GRID         110       0     2.4     -2.    -0.5       0        
GRID         111       0     2.4    -1.2    -0.5       0        
GRID         112       0     2.4     -1.    -0.1       0        
GRID         113       0     2.4    -1.6     0.3       0        
GRID         114       0     2.4    -1.2     0.1       0        
GRID         115       0     2.4    -1.4     0.1       0        
GRID         116       0     2.4    -1.8     0.1       0        
GRID         117       0     2.4    -1.6    -0.1       0        
GRID         118       0     2.4    -1.2    -0.3       0        
GRID         119       0     2.4    -1.6    -0.3       0        
GRID         120       0     1.4     -1.     0.5       0        
GRID         121       0     2.4     -1.     0.3       0        
GRID         122       0     1.2     -1.     0.1       0        
GRID         123       0     1.2     -1.     0.3       0        
GRID         124       0     1.4     -1.     0.3       0        
GRID         125       0     2.2     -1.     0.3       0        
GRID         126       0     1.4     -1.     0.1       0        
GRID         127       0     1.6     -1.     0.1       0        
GRID         128       0      2.     -1.     0.1       0        
GRID         129       0     2.2     -1.     0.1       0        
GRID         130       0     1.8     -1.    -0.1       0        
GRID         131       0     1.6     -1.    -0.3       0        
GRID         132       0     1.8     -1.    -0.3       0        
GRID         133       0     0.2      0.    -0.5       0        
GRID         134       0     0.4      0.    -0.5       0        
GRID         135       0     0.6      0.    -0.5       0        
GRID         136       0     0.8      0.    -0.5       0        
GRID         137       0      1.      0.    -0.5       0        
GRID         138       0     0.8      0.     0.5       0        
GRID         139       0     0.6      0.     0.5       0        
GRID         140       0     0.4      0.     0.5       0        
GRID         141       0      0.      0.    -0.1       0        
GRID         142       0     0.2      0.    -0.3       0        
GRID         143       0     0.4      0.    -0.3       0        
GRID         144       0     0.6      0.    -0.3       0        
GRID         145       0     0.8      0.    -0.3       0        
GRID         146       0      1.      0.    -0.3       0        
GRID         147       0     0.2      0.    -0.1       0        
GRID         148       0     0.4      0.    -0.1       0        
GRID         149       0     0.6      0.    -0.1       0        
GRID         150       0     0.8      0.    -0.1       0        
GRID         151       0      1.      0.    -0.1       0        
GRID         152       0     0.2      0.     0.1       0        
GRID         153       0     0.4      0.     0.1       0        
GRID         154       0     0.6      0.     0.1       0        
GRID         155       0     0.8      0.     0.1       0        
GRID         156       0      1.      0.     0.1       0        
GRID         157       0     0.2      0.     0.3       0        
GRID         158       0     0.4      0.     0.3       0        
GRID         159       0     0.6      0.     0.3       0        
GRID         160       0     0.8      0.     0.3       0        
GRID         161       0      1.      0.     0.3       0        
GRID         162       0     0.8      1.     0.5       0        
GRID         163       0      1.      1.     0.5       0        
GRID         164       0     0.8      1.    -0.5       0        
GRID         165       0     0.6      1.    -0.5       0        
GRID         166       0     0.4      1.    -0.5       0        
GRID         167       0     0.2      1.    -0.5       0        
GRID         168       0      0.      1.    -0.1       0        
GRID         169       0      0.      1.     0.1       0        
GRID         170       0      0.      1.     0.3       0        
GRID         171       0     0.2      1.     0.3       0        
GRID         172       0     0.4      1.     0.3       0        
GRID         173       0     0.6      1.     0.3       0        
GRID         174       0     0.8      1.     0.3       0        
GRID         175       0      1.      1.     0.3       0        
GRID         176       0     0.2      1.     0.1       0        
GRID         177       0     0.4      1.     0.1       0        
GRID         178       0     0.6      1.     0.1       0        
GRID         179       0     0.8      1.     0.1       0        
GRID         180       0      1.      1.     0.1       0        
GRID         181       0     0.2      1.    -0.1       0        
GRID         182       0     0.4      1.    -0.1       0        
GRID         183       0     0.6      1.    -0.1       0        
GRID         184       0     0.8      1.    -0.1       0        
GRID         185       0      1.      1.    -0.1       0        
GRID         186       0     0.2      1.    -0.3       0        
GRID         187       0     0.4      1.    -0.3       0        
GRID         188       0     0.6      1.    -0.3       0        
GRID         189       0     0.8      1.    -0.3       0        
GRID         190       0      1.      1.    -0.3       0        
GRID         191       0      0.     0.8    -0.5       0        
GRID         192       0      0.     0.4    -0.5       0        
GRID         193       0      0.      0.    -0.5       0        
GRID         194       0      0.      0.    -0.3       0        
GRID         195       0      0.      0.     0.1       0        
GRID         196       0      0.      0.     0.3       0        
GRID         197       0      0.     0.4     0.5       0        
GRID         198       0      0.     0.8     0.5       0        
GRID         199       0      0.      1.     0.5       0        
GRID         200       0      0.      1.    -0.3       0        
GRID         201       0      0.     0.8    -0.3       0        
GRID         202       0      0.     0.6    -0.3       0        
GRID         203       0      0.     0.4    -0.3       0        
GRID         204       0      0.     0.2    -0.3       0        
GRID         205       0      0.     0.8    -0.1       0        
GRID         206       0      0.     0.6    -0.1       0        
GRID         207       0      0.     0.4    -0.1       0        
GRID         208       0      0.     0.2    -0.1       0        
GRID         209       0      0.     0.8     0.1       0        
GRID         210       0      0.     0.6     0.1       0        
GRID         211       0      0.     0.4     0.1       0        
GRID         212       0      0.     0.2     0.1       0        
GRID         213       0      0.     0.8     0.3       0        
GRID         214       0      0.     0.6     0.3       0        
GRID         215       0      0.     0.4     0.3       0        
GRID         216       0      0.     0.2     0.3       0        
GRID         217       0      0.      0.     0.5       0        
GRID         218       0     0.2      0.     0.5       0        
GRID         219       0      1.      0.     0.5       0        
GRID         220       0     1.2     0.6     0.5       0        
GRID         221       0     1.2     0.8     0.5       0        
GRID         222       0     0.6      1.     0.5       0        
GRID         223       0     0.4      1.     0.5       0        
GRID         224       0     0.2      1.     0.5       0        
GRID         225       0      0.     0.6     0.5       0        
GRID         226       0      0.     0.2     0.5       0        
GRID         227       0     0.2     0.2     0.5       0        
GRID         228       0     0.4     0.2     0.5       0        
GRID         229       0     0.6     0.2     0.5       0        
GRID         230       0     0.8     0.2     0.5       0        
GRID         231       0      1.     0.2     0.5       0        
GRID         232       0     0.2     0.4     0.5       0        
GRID         233       0     0.4     0.4     0.5       0        
GRID         234       0     0.6     0.4     0.5       0        
GRID         235       0     0.8     0.4     0.5       0        
GRID         236       0      1.     0.4     0.5       0        
GRID         237       0     0.2     0.6     0.5       0        
GRID         238       0     0.4     0.6     0.5       0        
GRID         239       0     0.6     0.6     0.5       0        
GRID         240       0     0.8     0.6     0.5       0        
GRID         241       0      1.     0.6     0.5       0        
GRID         242       0     0.2     0.8     0.5       0        
GRID         243       0     0.4     0.8     0.5       0        
GRID         244       0     0.6     0.8     0.5       0        
GRID         245       0     0.8     0.8     0.5       0        
GRID         246       0      1.     0.8     0.5       0        
GRID         247       0     1.2      1.    -0.5       0        
GRID         248       0     1.2     0.8    -0.5       0        
GRID         249       0     1.2     0.6    -0.5       0        
GRID         250       0     1.2     0.4    -0.5       0        
GRID         251       0      0.     0.2    -0.5       0        
GRID         252       0      0.     0.6    -0.5       0        
GRID         253       0      0.      1.    -0.5       0        
GRID         254       0      1.      1.    -0.5       0        
GRID         255       0      1.     0.8    -0.5       0        
GRID         256       0      1.     0.6    -0.5       0        
GRID         257       0      1.     0.4    -0.5       0        
GRID         258       0      1.     0.2    -0.5       0        
GRID         259       0     0.8     0.8    -0.5       0        
GRID         260       0     0.8     0.6    -0.5       0        
GRID         261       0     0.8     0.4    -0.5       0        
GRID         262       0     0.8     0.2    -0.5       0        
GRID         263       0     0.6     0.8    -0.5       0        
GRID         264       0     0.6     0.6    -0.5       0        
GRID         265       0     0.6     0.4    -0.5       0        
GRID         266       0     0.6     0.2    -0.5       0        
GRID         267       0     0.4     0.8    -0.5       0        
GRID         268       0     0.4     0.6    -0.5       0        
GRID         269       0     0.4     0.4    -0.5       0        
GRID         270       0     0.4     0.2    -0.5       0        
GRID         271       0     0.2     0.8    -0.5       0        
GRID         272       0     0.2     0.6    -0.5       0        
GRID         273       0     0.2     0.4    -0.5       0        
GRID         274       0     0.2     0.2    -0.5       0        
GRID         275       0     1.2      0.     0.5       0        
GRID         276       0     1.2      0.    -0.1       0        
GRID         277       0     1.2     0.2    -0.5       0        
GRID         278       0     1.2      1.    -0.3       0        
GRID         279       0     1.2     0.4     0.5       0        
GRID         280       0     1.2     0.2     0.3       0        
GRID         281       0     1.2     0.2     0.1       0        
GRID         282       0     1.2     0.2    -0.1       0        
GRID         283       0     1.2     0.2    -0.3       0        
GRID         284       0     1.2     0.4     0.3       0        
GRID         285       0     1.2     0.4    -0.3       0        
GRID         286       0     1.2     0.6    -0.1       0        
GRID         287       0     1.2     0.6    -0.3       0        
GRID         288       0     1.2     0.8     0.3       0        
GRID         289       0     1.2     0.8    -0.3       0        
GRID         290       0     1.2      1.    -0.1       0        
GRID         291       0     1.2      1.     0.1       0        
GRID         292       0     1.2      1.     0.3       0        
GRID         293       0     1.2      1.     0.5       0        
GRID         294       0     1.6      1.     0.5       0        
GRID         295       0     1.8      1.     0.5       0        
GRID         296       0      2.      1.     0.5       0        
GRID         297       0     2.2      1.     0.5       0        
GRID         298       0     2.4      1.     0.5       0        
GRID         299       0     2.4      1.     0.1       0        
GRID         300       0     2.4      1.    -0.1       0        
GRID         301       0     2.4      1.    -0.3       0        
GRID         302       0     2.2      1.    -0.5       0        
GRID         303       0      2.      1.    -0.5       0        
GRID         304       0     1.4      1.    -0.3       0        
GRID         305       0     1.4      1.    -0.1       0        
GRID         306       0     1.4      1.     0.1       0        
GRID         307       0     1.4      1.     0.3       0        
GRID         308       0     1.6      1.    -0.3       0        
GRID         309       0     1.6      1.    -0.1       0        
GRID         310       0     1.6      1.     0.1       0        
GRID         311       0     1.6      1.     0.3       0        
GRID         312       0     1.8      1.    -0.3       0        
GRID         313       0     1.8      1.    -0.1       0        
GRID         314       0     1.8      1.     0.1       0        
GRID         315       0     1.8      1.     0.3       0        
GRID         316       0      2.      1.    -0.3       0        
GRID         317       0      2.      1.    -0.1       0        
GRID         318       0      2.      1.     0.1       0        
GRID         319       0      2.      1.     0.3       0        
GRID         320       0     2.2      1.    -0.3       0        
GRID         321       0     2.2      1.    -0.1       0        
GRID         322       0     2.2      1.     0.1       0        
GRID         323       0     2.2      1.     0.3       0        
GRID         324       0     1.2     0.4    -0.1       0        
GRID         325       0     1.2     0.4     0.1       0        
GRID         326       0     1.2     0.6     0.1       0        
GRID         327       0     1.2     0.6     0.3       0        
GRID         328       0     1.2     0.8    -0.1       0        
GRID         329       0     1.2     0.8     0.1       0        
GRID         330       0     1.4      0.    -0.5       0        
GRID         331       0     1.6      0.    -0.5       0        
GRID         332       0     2.4      0.    -0.3       0        
GRID         333       0     1.4      0.     0.5       0        
GRID         334       0     1.4      0.    -0.1       0        
GRID         335       0     2.2      0.    -0.1       0        
GRID         336       0     1.6      0.     0.1       0        
GRID         337       0     2.2      0.     0.1       0        
GRID         338       0     1.4      0.     0.3       0        
GRID         339       0     1.8      0.     0.3       0        
GRID         340       0      2.      0.     0.3       0        
GRID         341       0     2.4     0.4     0.5       0        
GRID         342       0     2.4     0.6     0.5       0        
GRID         343       0     2.4     0.8     0.5       0        
GRID         344       0     1.4      1.     0.5       0        
GRID         345       0     1.2     0.2     0.5       0        
GRID         346       0     1.4     0.2     0.5       0        
GRID         347       0     1.6     0.2     0.5       0        
GRID         348       0     1.8     0.2     0.5       0        
GRID         349       0      2.     0.2     0.5       0        
GRID         350       0     2.2     0.2     0.5       0        
GRID         351       0     1.4     0.4     0.5       0        
GRID         352       0     1.6     0.4     0.5       0        
GRID         353       0     1.8     0.4     0.5       0        
GRID         354       0      2.     0.4     0.5       0        
GRID         355       0     2.2     0.4     0.5       0        
GRID         356       0     1.4     0.6     0.5       0        
GRID         357       0     1.6     0.6     0.5       0        
GRID         358       0     1.8     0.6     0.5       0        
GRID         359       0      2.     0.6     0.5       0        
GRID         360       0     2.2     0.6     0.5       0        
GRID         361       0     1.4     0.8     0.5       0        
GRID         362       0     1.6     0.8     0.5       0        
GRID         363       0     1.8     0.8     0.5       0        
GRID         364       0      2.     0.8     0.5       0        
GRID         365       0     2.2     0.8     0.5       0        
GRID         366       0     2.4      1.    -0.5       0        
GRID         367       0     2.4     0.4    -0.5       0        
GRID         368       0     1.8      0.    -0.5       0        
GRID         369       0     1.4      1.    -0.5       0        
GRID         370       0     1.6      1.    -0.5       0        
GRID         371       0     1.8      1.    -0.5       0        
GRID         372       0     2.2     0.8    -0.5       0        
GRID         373       0     2.2     0.6    -0.5       0        
GRID         374       0     2.2     0.4    -0.5       0        
GRID         375       0     2.2     0.2    -0.5       0        
GRID         376       0      2.     0.8    -0.5       0        
GRID         377       0      2.     0.6    -0.5       0        
GRID         378       0      2.     0.4    -0.5       0        
GRID         379       0      2.     0.2    -0.5       0        
GRID         380       0     1.8     0.8    -0.5       0        
GRID         381       0     1.8     0.6    -0.5       0        
GRID         382       0     1.8     0.4    -0.5       0        
GRID         383       0     1.8     0.2    -0.5       0        
GRID         384       0     1.6     0.8    -0.5       0        
GRID         385       0     1.6     0.6    -0.5       0        
GRID         386       0     1.6     0.4    -0.5       0        
GRID         387       0     1.6     0.2    -0.5       0        
GRID         388       0     1.4     0.8    -0.5       0        
GRID         389       0     1.4     0.6    -0.5       0        
GRID         390       0     1.4     0.4    -0.5       0        
GRID         391       0     1.4     0.2    -0.5       0        
GRID         392       0     2.4     0.2    -0.5       0        
GRID         393       0     2.4     0.6    -0.5       0        
GRID         394       0     2.4     0.8    -0.5       0        
GRID         395       0     2.4      1.     0.3       0        
GRID         396       0     2.4     0.2     0.5       0        
GRID         397       0     2.4     0.2     0.3       0        
GRID         398       0     2.4     0.2     0.1       0        
GRID         399       0     2.4     0.2    -0.1       0        
GRID         400       0     2.4     0.2    -0.3       0        
GRID         401       0     2.4     0.4     0.3       0        
GRID         402       0     2.4     0.4     0.1       0        
GRID         403       0     2.4     0.4    -0.1       0        
GRID         404       0     2.4     0.4    -0.3       0        
GRID         405       0     2.4     0.6     0.3       0        
GRID         406       0     2.4     0.6     0.1       0        
GRID         407       0     2.4     0.6    -0.1       0        
GRID         408       0     2.4     0.6    -0.3       0        
GRID         409       0     2.4     0.8     0.3       0        
GRID         410       0     2.4     0.8     0.1       0        
GRID         411       0     2.4     0.8    -0.1       0        
GRID         412       0     2.4     0.8    -0.3       0        
GRID         413       0     3.6     -2.    -0.5       0        
GRID         414       0     3.6    -1.6    -0.5       0        
GRID         415       0     3.6    -1.4    -0.5       0        
GRID         416       0     3.6    -1.2    -0.5       0        
GRID         417       0     3.6     -1.    -0.5       0        
GRID         418       0     3.6     -1.    -0.1       0        
GRID         419       0     3.6     -1.     0.5       0        
GRID         420       0     3.6    -1.2     0.5       0        
GRID         421       0     3.6    -1.4     0.5       0        
GRID         422       0     3.6    -1.8     0.5       0        
GRID         423       0     3.6     -2.     0.1       0        
GRID         424       0     3.6    -1.8    -0.3       0        
GRID         425       0     3.6    -1.6    -0.3       0        
GRID         426       0     3.6    -1.4    -0.3       0        
GRID         427       0     3.6    -1.2    -0.3       0        
GRID         428       0     3.6    -1.8    -0.1       0        
GRID         429       0     3.6    -1.6    -0.1       0        
GRID         430       0     3.6    -1.4    -0.1       0        
GRID         431       0     3.6    -1.2    -0.1       0        
GRID         432       0     3.6    -1.8     0.1       0        
GRID         433       0     3.6    -1.6     0.1       0        
GRID         434       0     3.6    -1.4     0.1       0        
GRID         435       0     3.6    -1.2     0.1       0        
GRID         436       0     3.6    -1.8     0.3       0        
GRID         437       0     3.6    -1.6     0.3       0        
GRID         438       0     3.6    -1.4     0.3       0        
GRID         439       0     3.6    -1.2     0.3       0        
GRID         440       0      3.     -1.    -0.5       0        
GRID         441       0     2.6     -1.    -0.5       0        
GRID         442       0     2.6     -1.     0.5       0        
GRID         443       0     2.8     -1.     0.5       0        
GRID         444       0      3.     -1.     0.5       0        
GRID         445       0     3.2     -1.     0.5       0        
GRID         446       0     3.6     -1.     0.3       0        
GRID         447       0     3.6     -1.     0.1       0        
GRID         448       0     3.6     -1.    -0.3       0        
GRID         449       0     3.4     -1.    -0.3       0        
GRID         450       0     3.2     -1.    -0.3       0        
GRID         451       0      3.     -1.    -0.3       0        
GRID         452       0     2.8     -1.    -0.3       0        
GRID         453       0     2.6     -1.    -0.3       0        
GRID         454       0     3.4     -1.    -0.1       0        
GRID         455       0     3.2     -1.    -0.1       0        
GRID         456       0      3.     -1.    -0.1       0        
GRID         457       0     2.8     -1.    -0.1       0        
GRID         458       0     2.6     -1.    -0.1       0        
GRID         459       0     3.4     -1.     0.1       0        
GRID         460       0     3.2     -1.     0.1       0        
GRID         461       0      3.     -1.     0.1       0        
GRID         462       0     2.8     -1.     0.1       0        
GRID         463       0     2.6     -1.     0.1       0        
GRID         464       0     3.4     -1.     0.3       0        
GRID         465       0     3.2     -1.     0.3       0        
GRID         466       0      3.     -1.     0.3       0        
GRID         467       0     2.8     -1.     0.3       0        
GRID         468       0     2.6     -1.     0.3       0        
GRID         469       0     2.4     -2.     0.5       0        
GRID         470       0     2.4     -2.    -0.3       0        
GRID         471       0     2.6     -2.    -0.5       0        
GRID         472       0     2.8     -2.    -0.5       0        
GRID         473       0     3.6     -2.    -0.3       0        
GRID         474       0     3.6     -2.    -0.1       0        
GRID         475       0     3.6     -2.     0.3       0        
GRID         476       0     3.4     -2.     0.5       0        
GRID         477       0     3.2     -2.     0.5       0        
GRID         478       0     2.8     -2.     0.5       0        
GRID         479       0     2.6     -2.     0.3       0        
GRID         480       0     2.6     -2.     0.1       0        
GRID         481       0     2.6     -2.    -0.1       0        
GRID         482       0     2.6     -2.    -0.3       0        
GRID         483       0     2.8     -2.     0.3       0        
GRID         484       0     2.8     -2.     0.1       0        
GRID         485       0     2.8     -2.    -0.1       0        
GRID         486       0     2.8     -2.    -0.3       0        
GRID         487       0      3.     -2.     0.3       0        
GRID         488       0      3.     -2.     0.1       0        
GRID         489       0      3.     -2.    -0.1       0        
GRID         490       0      3.     -2.    -0.3       0        
GRID         491       0     3.2     -2.     0.3       0        
GRID         492       0     3.2     -2.     0.1       0        
GRID         493       0     3.2     -2.    -0.1       0        
GRID         494       0     3.2     -2.    -0.3       0        
GRID         495       0     3.4     -2.     0.3       0        
GRID         496       0     3.4     -2.     0.1       0        
GRID         497       0     3.4     -2.    -0.1       0        
GRID         498       0     3.4     -2.    -0.3       0        
GRID         499       0     3.4     -2.    -0.5       0        
GRID         500       0     3.2     -2.    -0.5       0        
GRID         501       0      3.     -2.    -0.5       0        
GRID         502       0     2.4    -1.6    -0.5       0        
GRID         503       0     2.4     -1.    -0.5       0        
GRID         504       0     2.8     -1.    -0.5       0        
GRID         505       0     3.2     -1.    -0.5       0        
GRID         506       0     3.4     -1.    -0.5       0        
GRID         507       0     3.6    -1.8    -0.5       0        
GRID         508       0     3.4    -1.8    -0.5       0        
GRID         509       0     3.2    -1.8    -0.5       0        
GRID         510       0      3.    -1.8    -0.5       0        
GRID         511       0     2.8    -1.8    -0.5       0        
GRID         512       0     2.6    -1.8    -0.5       0        
GRID         513       0     3.4    -1.6    -0.5       0        
GRID         514       0     3.2    -1.6    -0.5       0        
GRID         515       0      3.    -1.6    -0.5       0        
GRID         516       0     2.8    -1.6    -0.5       0        
GRID         517       0     2.6    -1.6    -0.5       0        
GRID         518       0     3.4    -1.4    -0.5       0        
GRID         519       0     3.2    -1.4    -0.5       0        
GRID         520       0      3.    -1.4    -0.5       0        
GRID         521       0     2.8    -1.4    -0.5       0        
GRID         522       0     2.6    -1.4    -0.5       0        
GRID         523       0     3.4    -1.2    -0.5       0        
GRID         524       0     3.2    -1.2    -0.5       0        
GRID         525       0      3.    -1.2    -0.5       0        
GRID         526       0     2.8    -1.2    -0.5       0        
GRID         527       0     2.6    -1.2    -0.5       0        
GRID         528       0     2.4    -1.2     0.5       0        
GRID         529       0     2.6     -2.     0.5       0        
GRID         530       0      3.     -2.     0.5       0        
GRID         531       0     3.6     -2.     0.5       0        
GRID         532       0     3.6    -1.6     0.5       0        
GRID         533       0     3.4     -1.     0.5       0        
GRID         534       0     2.6    -1.2     0.5       0        
GRID         535       0     2.6    -1.4     0.5       0        
GRID         536       0     2.6    -1.6     0.5       0        
GRID         537       0     2.6    -1.8     0.5       0        
GRID         538       0     2.8    -1.2     0.5       0        
GRID         539       0     2.8    -1.4     0.5       0        
GRID         540       0     2.8    -1.6     0.5       0        
GRID         541       0     2.8    -1.8     0.5       0        
GRID         542       0      3.    -1.2     0.5       0        
GRID         543       0      3.    -1.4     0.5       0        
GRID         544       0      3.    -1.6     0.5       0        
GRID         545       0      3.    -1.8     0.5       0        
GRID         546       0     3.2    -1.2     0.5       0        
GRID         547       0     3.2    -1.4     0.5       0        
GRID         548       0     3.2    -1.6     0.5       0        
GRID         549       0     3.2    -1.8     0.5       0        
GRID         550       0     3.4    -1.2     0.5       0        
GRID         551       0     3.4    -1.4     0.5       0        
GRID         552       0     3.4    -1.6     0.5       0        
GRID         553       0     3.4    -1.8     0.5       0        
GRID         554       0     2.4     -1.     0.5       0        
GRID         555       0     2.4     -1.    -0.3       0        
GRID         556       0     2.4    -1.4    -0.5       0        
GRID         557       0     2.4     -2.     0.1       0        
GRID         558       0     2.4    -1.8     0.3       0        
GRID         559       0     2.4    -1.4     0.3       0        
GRID         560       0     2.4    -1.2     0.3       0        
GRID         561       0     2.4    -1.6     0.1       0        
GRID         562       0     2.4    -1.8    -0.1       0        
GRID         563       0     2.4    -1.4    -0.1       0        
GRID         564       0     2.4    -1.2    -0.1       0        
GRID         565       0     2.4    -1.8    -0.3       0        
GRID         566       0     2.4    -1.4    -0.3       0        
GRID         567       0     2.4      0.    -0.5       0        
GRID         568       0     2.4      0.     0.3       0        
GRID         569       0     2.4    -0.4     0.5       0        
GRID         570       0     2.4    -0.8     0.5       0        
GRID         571       0     2.4     -1.     0.1       0        
GRID         572       0     2.4    -0.8    -0.3       0        
GRID         573       0     2.4    -0.6    -0.3       0        
GRID         574       0     2.4    -0.4    -0.3       0        
GRID         575       0     2.4    -0.2    -0.3       0        
GRID         576       0     2.4    -0.8    -0.1       0        
GRID         577       0     2.4    -0.6    -0.1       0        
GRID         578       0     2.4    -0.4    -0.1       0        
GRID         579       0     2.4    -0.2    -0.1       0        
GRID         580       0     2.4    -0.8     0.1       0        
GRID         581       0     2.4    -0.6     0.1       0        
GRID         582       0     2.4    -0.4     0.1       0        
GRID         583       0     2.4    -0.2     0.1       0        
GRID         584       0     2.4    -0.8     0.3       0        
GRID         585       0     2.4    -0.6     0.3       0        
GRID         586       0     2.4    -0.4     0.3       0        
GRID         587       0     2.4    -0.2     0.3       0        
GRID         588       0     2.2      0.    -0.5       0        
GRID         589       0     1.2      0.    -0.5       0        
GRID         590       0     1.2      0.    -0.3       0        
GRID         591       0     1.2      0.     0.1       0        
GRID         592       0     1.2      0.     0.3       0        
GRID         593       0     1.6      0.     0.5       0        
GRID         594       0     1.8      0.     0.5       0        
GRID         595       0     2.2      0.     0.5       0        
GRID         596       0     2.4      0.     0.1       0        
GRID         597       0     2.4      0.    -0.1       0        
GRID         598       0     2.2      0.    -0.3       0        
GRID         599       0      2.      0.    -0.3       0        
GRID         600       0     1.8      0.    -0.3       0        
GRID         601       0     1.6      0.    -0.3       0        
GRID         602       0     1.4      0.    -0.3       0        
GRID         603       0      2.      0.    -0.1       0        
GRID         604       0     1.8      0.    -0.1       0        
GRID         605       0     1.6      0.    -0.1       0        
GRID         606       0      2.      0.     0.1       0        
GRID         607       0     1.8      0.     0.1       0        
GRID         608       0     1.4      0.     0.1       0        
GRID         609       0     2.2      0.     0.3       0        
GRID         610       0     1.6      0.     0.3       0        
GRID         611       0      2.     -1.     0.5       0        
GRID         612       0     1.2     -1.     0.5       0        
GRID         613       0     1.2     -1.    -0.5       0        
GRID         614       0     1.8     -1.    -0.5       0        
GRID         615       0      2.     -1.     0.3       0        
GRID         616       0     1.8     -1.     0.3       0        
GRID         617       0     1.6     -1.     0.3       0        
GRID         618       0     1.8     -1.     0.1       0        
GRID         619       0     2.2     -1.    -0.1       0        
GRID         620       0      2.     -1.    -0.1       0        
GRID         621       0     1.6     -1.    -0.1       0        
GRID         622       0     1.4     -1.    -0.1       0        
GRID         623       0     2.2     -1.    -0.3       0        
GRID         624       0      2.     -1.    -0.3       0        
GRID         625       0     1.4     -1.    -0.3       0        
GRID         626       0     2.4    -0.2    -0.5       0        
GRID         627       0     2.4    -0.4    -0.5       0        
GRID         628       0     2.4    -0.6    -0.5       0        
GRID         629       0     2.4    -0.8    -0.5       0        
GRID         630       0      2.     -1.    -0.5       0        
GRID         631       0     1.4     -1.    -0.5       0        
GRID         632       0     1.2    -0.8    -0.5       0        
GRID         633       0     1.2    -0.6    -0.5       0        
GRID         634       0     1.2    -0.4    -0.5       0        
GRID         635       0      2.      0.    -0.5       0        
GRID         636       0     2.2    -0.2    -0.5       0        
GRID         637       0     2.2    -0.4    -0.5       0        
GRID         638       0     2.2    -0.6    -0.5       0        
GRID         639       0     2.2    -0.8    -0.5       0        
GRID         640       0      2.    -0.2    -0.5       0        
GRID         641       0      2.    -0.4    -0.5       0        
GRID         642       0      2.    -0.6    -0.5       0        
GRID         643       0      2.    -0.8    -0.5       0        
GRID         644       0     1.8    -0.2    -0.5       0        
GRID         645       0     1.8    -0.4    -0.5       0        
GRID         646       0     1.8    -0.6    -0.5       0        
GRID         647       0     1.8    -0.8    -0.5       0        
GRID         648       0     1.6    -0.2    -0.5       0        
GRID         649       0     1.6    -0.4    -0.5       0        
GRID         650       0     1.6    -0.6    -0.5       0        
GRID         651       0     1.6    -0.8    -0.5       0        
GRID         652       0     1.4    -0.2    -0.5       0        
GRID         653       0     1.4    -0.4    -0.5       0        
GRID         654       0     1.4    -0.6    -0.5       0        
GRID         655       0     1.4    -0.8    -0.5       0        
GRID         656       0     1.6     -1.     0.5       0        
GRID         657       0     1.8     -1.     0.5       0        
GRID         658       0     2.2     -1.     0.5       0        
GRID         659       0     2.4    -0.6     0.5       0        
GRID         660       0     2.4    -0.2     0.5       0        
GRID         661       0     2.4      0.     0.5       0        
GRID         662       0      2.      0.     0.5       0        
GRID         663       0     1.2    -0.6     0.5       0        
GRID         664       0     1.4    -0.8     0.5       0        
GRID         665       0     1.6    -0.8     0.5       0        
GRID         666       0     1.8    -0.8     0.5       0        
GRID         667       0      2.    -0.8     0.5       0        
GRID         668       0     2.2    -0.8     0.5       0        
GRID         669       0     1.4    -0.6     0.5       0        
GRID         670       0     1.6    -0.6     0.5       0        
GRID         671       0     1.8    -0.6     0.5       0        
GRID         672       0      2.    -0.6     0.5       0        
GRID         673       0     2.2    -0.6     0.5       0        
GRID         674       0     1.4    -0.4     0.5       0        
GRID         675       0     1.6    -0.4     0.5       0        
GRID         676       0     1.8    -0.4     0.5       0        
GRID         677       0      2.    -0.4     0.5       0        
GRID         678       0     2.2    -0.4     0.5       0        
GRID         679       0     1.4    -0.2     0.5       0        
GRID         680       0     1.6    -0.2     0.5       0        
GRID         681       0     1.8    -0.2     0.5       0        
GRID         682       0      2.    -0.2     0.5       0        
GRID         683       0     2.2    -0.2     0.5       0        
GRID         684       0     1.2     -1.    -0.1       0        
GRID         685       0     1.2    -0.8     0.5       0        
GRID         686       0     1.2    -0.4     0.5       0        
GRID         687       0     1.2    -0.2     0.5       0        
GRID         688       0     1.2    -0.2    -0.5       0        
GRID         689       0     1.2    -0.8    -0.3       0        
GRID         690       0     1.2    -0.8    -0.1       0        
GRID         691       0     1.2    -0.8     0.1       0        
GRID         692       0     1.2    -0.8     0.3       0        
GRID         693       0     1.2    -0.6    -0.3       0        
GRID         694       0     1.2    -0.6    -0.1       0        
GRID         695       0     1.2    -0.6     0.1       0        
GRID         696       0     1.2    -0.6     0.3       0        
GRID         697       0     1.2    -0.4    -0.3       0        
GRID         698       0     1.2    -0.4    -0.1       0        
GRID         699       0     1.2    -0.4     0.1       0        
GRID         700       0     1.2    -0.4     0.3       0        
GRID         701       0     1.2    -0.2    -0.3       0        
GRID         702       0     1.2    -0.2    -0.1       0        
GRID         703       0     1.2    -0.2     0.1       0        
GRID         704       0     1.2    -0.2     0.3       0        
GRID         705       0     1.4    -1.8     0.3       0        
GRID         706       0     1.4    -1.6     0.3       0        
GRID         707       0     1.4    -1.4     0.3       0        
GRID         708       0     1.4    -1.2     0.3       0        
GRID         709       0     1.4    -1.8     0.1       0        
GRID         710       0     1.4    -1.6     0.1       0        
GRID         711       0     1.4    -1.4     0.1       0        
GRID         712       0     1.4    -1.2     0.1       0        
GRID         713       0     1.4    -1.8    -0.1       0        
GRID         714       0     1.4    -1.6    -0.1       0        
GRID         715       0     1.4    -1.4    -0.1       0        
GRID         716       0     1.4    -1.2    -0.1       0        
GRID         717       0     1.4    -1.8    -0.3       0        
GRID         718       0     1.4    -1.6    -0.3       0        
GRID         719       0     1.4    -1.4    -0.3       0        
GRID         720       0     1.4    -1.2    -0.3       0        
GRID         721       0     1.6    -1.8     0.3       0        
GRID         722       0     1.6    -1.6     0.3       0        
GRID         723       0     1.6    -1.4     0.3       0        
GRID         724       0     1.6    -1.2     0.3       0        
GRID         725       0     1.6    -1.8     0.1       0        
GRID         726       0     1.6    -1.6     0.1       0        
GRID         727       0     1.6    -1.4     0.1       0        
GRID         728       0     1.6    -1.2     0.1       0        
GRID         729       0     1.6    -1.8    -0.1       0        
GRID         730       0     1.6    -1.6    -0.1       0        
GRID         731       0     1.6    -1.4    -0.1       0        
GRID         732       0     1.6    -1.2    -0.1       0        
GRID         733       0     1.6    -1.8    -0.3       0        
GRID         734       0     1.6    -1.6    -0.3       0        
GRID         735       0     1.6    -1.4    -0.3       0        
GRID         736       0     1.6    -1.2    -0.3       0        
GRID         737       0     1.8    -1.8     0.3       0        
GRID         738       0     1.8    -1.6     0.3       0        
GRID         739       0     1.8    -1.4     0.3       0        
GRID         740       0     1.8    -1.2     0.3       0        
GRID         741       0     1.8    -1.8     0.1       0        
GRID         742       0     1.8    -1.6     0.1       0        
GRID         743       0     1.8    -1.4     0.1       0        
GRID         744       0     1.8    -1.2     0.1       0        
GRID         745       0     1.8    -1.8    -0.1       0        
GRID         746       0     1.8    -1.6    -0.1       0        
GRID         747       0     1.8    -1.4    -0.1       0        
GRID         748       0     1.8    -1.2    -0.1       0        
GRID         749       0     1.8    -1.8    -0.3       0        
GRID         750       0     1.8    -1.6    -0.3       0        
GRID         751       0     1.8    -1.4    -0.3       0        
GRID         752       0     1.8    -1.2    -0.3       0        
GRID         753       0      2.    -1.8     0.3       0        
GRID         754       0      2.    -1.6     0.3       0        
GRID         755       0      2.    -1.4     0.3       0        
GRID         756       0      2.    -1.2     0.3       0        
GRID         757       0      2.    -1.8     0.1       0        
GRID         758       0      2.    -1.6     0.1       0        
GRID         759       0      2.    -1.4     0.1       0        
GRID         760       0      2.    -1.2     0.1       0        
GRID         761       0      2.    -1.8    -0.1       0        
GRID         762       0      2.    -1.6    -0.1       0        
GRID         763       0      2.    -1.4    -0.1       0        
GRID         764       0      2.    -1.2    -0.1       0        
GRID         765       0      2.    -1.8    -0.3       0        
GRID         766       0      2.    -1.6    -0.3       0        
GRID         767       0      2.    -1.4    -0.3       0        
GRID         768       0      2.    -1.2    -0.3       0        
GRID         769       0     2.2    -1.8     0.3       0        
GRID         770       0     2.2    -1.6     0.3       0        
GRID         771       0     2.2    -1.4     0.3       0        
GRID         772       0     2.2    -1.2     0.3       0        
GRID         773       0     2.2    -1.8     0.1       0        
GRID         774       0     2.2    -1.6     0.1       0        
GRID         775       0     2.2    -1.4     0.1       0        
GRID         776       0     2.2    -1.2     0.1       0        
GRID         777       0     2.2    -1.8    -0.1       0        
GRID         778       0     2.2    -1.6    -0.1       0        
GRID         779       0     2.2    -1.4    -0.1       0        
GRID         780       0     2.2    -1.2    -0.1       0        
GRID         781       0     2.2    -1.8    -0.3       0        
GRID         782       0     2.2    -1.6    -0.3       0        
GRID         783       0     2.2    -1.4    -0.3       0        
GRID         784       0     2.2    -1.2    -0.3       0        
GRID         785       0     0.2     0.2    -0.3       0        
GRID         786       0     0.4     0.2    -0.3       0        
GRID         787       0     0.6     0.2    -0.3       0        
GRID         788       0     0.8     0.2    -0.3       0        
GRID         789       0      1.     0.2    -0.3       0        
GRID         790       0     0.2     0.2    -0.1       0        
GRID         791       0     0.4     0.2    -0.1       0        
GRID         792       0     0.6     0.2    -0.1       0        
GRID         793       0     0.8     0.2    -0.1       0        
GRID         794       0      1.     0.2    -0.1       0        
GRID         795       0     0.2     0.2     0.1       0        
GRID         796       0     0.4     0.2     0.1       0        
GRID         797       0     0.6     0.2     0.1       0        
GRID         798       0     0.8     0.2     0.1       0        
GRID         799       0      1.     0.2     0.1       0        
GRID         800       0     0.2     0.2     0.3       0        
GRID         801       0     0.4     0.2     0.3       0        
GRID         802       0     0.6     0.2     0.3       0        
GRID         803       0     0.8     0.2     0.3       0        
GRID         804       0      1.     0.2     0.3       0        
GRID         805       0     0.2     0.4    -0.3       0        
GRID         806       0     0.4     0.4    -0.3       0        
GRID         807       0     0.6     0.4    -0.3       0        
GRID         808       0     0.8     0.4    -0.3       0        
GRID         809       0      1.     0.4    -0.3       0        
GRID         810       0     0.2     0.4    -0.1       0        
GRID         811       0     0.4     0.4    -0.1       0        
GRID         812       0     0.6     0.4    -0.1       0        
GRID         813       0     0.8     0.4    -0.1       0        
GRID         814       0      1.     0.4    -0.1       0        
GRID         815       0     0.2     0.4     0.1       0        
GRID         816       0     0.4     0.4     0.1       0        
GRID         817       0     0.6     0.4     0.1       0        
GRID         818       0     0.8     0.4     0.1       0        
GRID         819       0      1.     0.4     0.1       0        
GRID         820       0     0.2     0.4     0.3       0        
GRID         821       0     0.4     0.4     0.3       0        
GRID         822       0     0.6     0.4     0.3       0        
GRID         823       0     0.8     0.4     0.3       0        
GRID         824       0      1.     0.4     0.3       0        
GRID         825       0     0.2     0.6    -0.3       0        
GRID         826       0     0.4     0.6    -0.3       0        
GRID         827       0     0.6     0.6    -0.3       0        
GRID         828       0     0.8     0.6    -0.3       0        
GRID         829       0      1.     0.6    -0.3       0        
GRID         830       0     0.2     0.6    -0.1       0        
GRID         831       0     0.4     0.6    -0.1       0        
GRID         832       0     0.6     0.6    -0.1       0        
GRID         833       0     0.8     0.6    -0.1       0        
GRID         834       0      1.     0.6    -0.1       0        
GRID         835       0     0.2     0.6     0.1       0        
GRID         836       0     0.4     0.6     0.1       0        
GRID         837       0     0.6     0.6     0.1       0        
GRID         838       0     0.8     0.6     0.1       0        
GRID         839       0      1.     0.6     0.1       0        
GRID         840       0     0.2     0.6     0.3       0        
GRID         841       0     0.4     0.6     0.3       0        
GRID         842       0     0.6     0.6     0.3       0        
GRID         843       0     0.8     0.6     0.3       0        
GRID         844       0      1.     0.6     0.3       0        
GRID         845       0     0.2     0.8    -0.3       0        
GRID         846       0     0.4     0.8    -0.3       0        
GRID         847       0     0.6     0.8    -0.3       0        
GRID         848       0     0.8     0.8    -0.3       0        
GRID         849       0      1.     0.8    -0.3       0        
GRID         850       0     0.2     0.8    -0.1       0        
GRID         851       0     0.4     0.8    -0.1       0        
GRID         852       0     0.6     0.8    -0.1       0        
GRID         853       0     0.8     0.8    -0.1       0        
GRID         854       0      1.     0.8    -0.1       0        
GRID         855       0     0.2     0.8     0.1       0        
GRID         856       0     0.4     0.8     0.1       0        
GRID         857       0     0.6     0.8     0.1       0        
GRID         858       0     0.8     0.8     0.1       0        
GRID         859       0      1.     0.8     0.1       0        
GRID         860       0     0.2     0.8     0.3       0        
GRID         861       0     0.4     0.8     0.3       0        
GRID         862       0     0.6     0.8     0.3       0        
GRID         863       0     0.8     0.8     0.3       0        
GRID         864       0      1.     0.8     0.3       0        
GRID         865       0     1.4     0.2    -0.3       0        
GRID         866       0     1.6     0.2    -0.3       0        
GRID         867       0     1.8     0.2    -0.3       0        
GRID         868       0      2.     0.2    -0.3       0        
GRID         869       0     2.2     0.2    -0.3       0        
GRID         870       0     1.4     0.2    -0.1       0        
GRID         871       0     1.6     0.2    -0.1       0        
GRID         872       0     1.8     0.2    -0.1       0        
GRID         873       0      2.     0.2    -0.1       0        
GRID         874       0     2.2     0.2    -0.1       0        
GRID         875       0     1.4     0.2     0.1       0        
GRID         876       0     1.6     0.2     0.1       0        
GRID         877       0     1.8     0.2     0.1       0        
GRID         878       0      2.     0.2     0.1       0        
GRID         879       0     2.2     0.2     0.1       0        
GRID         880       0     1.4     0.2     0.3       0        
GRID         881       0     1.6     0.2     0.3       0        
GRID         882       0     1.8     0.2     0.3       0        
GRID         883       0      2.     0.2     0.3       0        
GRID         884       0     2.2     0.2     0.3       0        
GRID         885       0     1.4     0.4    -0.3       0        
GRID         886       0     1.6     0.4    -0.3       0        
GRID         887       0     1.8     0.4    -0.3       0        
GRID         888       0      2.     0.4    -0.3       0        
GRID         889       0     2.2     0.4    -0.3       0        
GRID         890       0     1.4     0.4    -0.1       0        
GRID         891       0     1.6     0.4    -0.1       0        
GRID         892       0     1.8     0.4    -0.1       0        
GRID         893       0      2.     0.4    -0.1       0        
GRID         894       0     2.2     0.4    -0.1       0        
GRID         895       0     1.4     0.4     0.1       0        
GRID         896       0     1.6     0.4     0.1       0        
GRID         897       0     1.8     0.4     0.1       0        
GRID         898       0      2.     0.4     0.1       0        
GRID         899       0     2.2     0.4     0.1       0        
GRID         900       0     1.4     0.4     0.3       0        
GRID         901       0     1.6     0.4     0.3       0        
GRID         902       0     1.8     0.4     0.3       0        
GRID         903       0      2.     0.4     0.3       0        
GRID         904       0     2.2     0.4     0.3       0        
GRID         905       0     1.4     0.6    -0.3       0        
GRID         906       0     1.6     0.6    -0.3       0        
GRID         907       0     1.8     0.6    -0.3       0        
GRID         908       0      2.     0.6    -0.3       0        
GRID         909       0     2.2     0.6    -0.3       0        
GRID         910       0     1.4     0.6    -0.1       0        
GRID         911       0     1.6     0.6    -0.1       0        
GRID         912       0     1.8     0.6    -0.1       0        
GRID         913       0      2.     0.6    -0.1       0        
GRID         914       0     2.2     0.6    -0.1       0        
GRID         915       0     1.4     0.6     0.1       0        
GRID         916       0     1.6     0.6     0.1       0        
GRID         917       0     1.8     0.6     0.1       0        
GRID         918       0      2.     0.6     0.1       0        
GRID         919       0     2.2     0.6     0.1       0        
GRID         920       0     1.4     0.6     0.3       0        
GRID         921       0     1.6     0.6     0.3       0        
GRID         922       0     1.8     0.6     0.3       0        
GRID         923       0      2.     0.6     0.3       0        
GRID         924       0     2.2     0.6     0.3       0        
GRID         925       0     1.4     0.8    -0.3       0        
GRID         926       0     1.6     0.8    -0.3       0        
GRID         927       0     1.8     0.8    -0.3       0        
GRID         928       0      2.     0.8    -0.3       0        
GRID         929       0     2.2     0.8    -0.3       0        
GRID         930       0     1.4     0.8    -0.1       0        
GRID         931       0     1.6     0.8    -0.1       0        
GRID         932       0     1.8     0.8    -0.1       0        
GRID         933       0      2.     0.8    -0.1       0        
GRID         934       0     2.2     0.8    -0.1       0        
GRID         935       0     1.4     0.8     0.1       0        
GRID         936       0     1.6     0.8     0.1       0        
GRID         937       0     1.8     0.8     0.1       0        
GRID         938       0      2.     0.8     0.1       0        
GRID         939       0     2.2     0.8     0.1       0        
GRID         940       0     1.4     0.8     0.3       0        
GRID         941       0     1.6     0.8     0.3       0        
GRID         942       0     1.8     0.8     0.3       0        
GRID         943       0      2.     0.8     0.3       0        
GRID         944       0     2.2     0.8     0.3       0        
GRID         945       0     2.6    -1.8     0.3       0        
GRID         946       0     2.6    -1.8     0.1       0        
GRID         947       0     2.6    -1.8    -0.1       0        
GRID         948       0     2.6    -1.8    -0.3       0        
GRID         949       0     2.8    -1.8     0.3       0        
GRID         950       0     2.8    -1.8     0.1       0        
GRID         951       0     2.8    -1.8    -0.1       0        
GRID         952       0     2.8    -1.8    -0.3       0        
GRID         953       0      3.    -1.8     0.3       0        
GRID         954       0      3.    -1.8     0.1       0        
GRID         955       0      3.    -1.8    -0.1       0        
GRID         956       0      3.    -1.8    -0.3       0        
GRID         957       0     3.2    -1.8     0.3       0        
GRID         958       0     3.2    -1.8     0.1       0        
GRID         959       0     3.2    -1.8    -0.1       0        
GRID         960       0     3.2    -1.8    -0.3       0        
GRID         961       0     3.4    -1.8     0.3       0        
GRID         962       0     3.4    -1.8     0.1       0        
GRID         963       0     3.4    -1.8    -0.1       0        
GRID         964       0     3.4    -1.8    -0.3       0        
GRID         965       0     2.6    -1.6     0.3       0        
GRID         966       0     2.6    -1.6     0.1       0        
GRID         967       0     2.6    -1.6    -0.1       0        
GRID         968       0     2.6    -1.6    -0.3       0        
GRID         969       0     2.8    -1.6     0.3       0        
GRID         970       0     2.8    -1.6     0.1       0        
GRID         971       0     2.8    -1.6    -0.1       0        
GRID         972       0     2.8    -1.6    -0.3       0        
GRID         973       0      3.    -1.6     0.3       0        
GRID         974       0      3.    -1.6     0.1       0        
GRID         975       0      3.    -1.6    -0.1       0        
GRID         976       0      3.    -1.6    -0.3       0        
GRID         977       0     3.2    -1.6     0.3       0        
GRID         978       0     3.2    -1.6     0.1       0        
GRID         979       0     3.2    -1.6    -0.1       0        
GRID         980       0     3.2    -1.6    -0.3       0        
GRID         981       0     3.4    -1.6     0.3       0        
GRID         982       0     3.4    -1.6     0.1       0        
GRID         983       0     3.4    -1.6    -0.1       0        
GRID         984       0     3.4    -1.6    -0.3       0        
GRID         985       0     2.6    -1.4     0.3       0        
GRID         986       0     2.6    -1.4     0.1       0        
GRID         987       0     2.6    -1.4    -0.1       0        
GRID         988       0     2.6    -1.4    -0.3       0        
GRID         989       0     2.8    -1.4     0.3       0        
GRID         990       0     2.8    -1.4     0.1       0        
GRID         991       0     2.8    -1.4    -0.1       0        
GRID         992       0     2.8    -1.4    -0.3       0        
GRID         993       0      3.    -1.4     0.3       0        
GRID         994       0      3.    -1.4     0.1       0        
GRID         995       0      3.    -1.4    -0.1       0        
GRID         996       0      3.    -1.4    -0.3       0        
GRID         997       0     3.2    -1.4     0.3       0        
GRID         998       0     3.2    -1.4     0.1       0        
GRID         999       0     3.2    -1.4    -0.1       0        
GRID        1000       0     3.2    -1.4    -0.3       0        
GRID        1001       0     3.4    -1.4     0.3       0        
GRID        1002       0     3.4    -1.4     0.1       0        
GRID        1003       0     3.4    -1.4    -0.1       0        
GRID        1004       0     3.4    -1.4    -0.3       0        
GRID        1005       0     2.6    -1.2     0.3       0        
GRID        1006       0     2.6    -1.2     0.1       0        
GRID        1007       0     2.6    -1.2    -0.1       0        
GRID        1008       0     2.6    -1.2    -0.3       0        
GRID        1009       0     2.8    -1.2     0.3       0        
GRID        1010       0     2.8    -1.2     0.1       0        
GRID        1011       0     2.8    -1.2    -0.1       0        
GRID        1012       0     2.8    -1.2    -0.3       0        
GRID        1013       0      3.    -1.2     0.3       0        
GRID        1014       0      3.    -1.2     0.1       0        
GRID        1015       0      3.    -1.2    -0.1       0        
GRID        1016       0      3.    -1.2    -0.3       0        
GRID        1017       0     3.2    -1.2     0.3       0        
GRID        1018       0     3.2    -1.2     0.1       0        
GRID        1019       0     3.2    -1.2    -0.1       0        
GRID        1020       0     3.2    -1.2    -0.3       0        
GRID        1021       0     3.4    -1.2     0.3       0        
GRID        1022       0     3.4    -1.2     0.1       0        
GRID        1023       0     3.4    -1.2    -0.1       0        
GRID        1024       0     3.4    -1.2    -0.3       0        
GRID        1025       0     2.2    -0.8     0.3       0        
GRID        1026       0      2.    -0.8     0.3       0        
GRID        1027       0     1.8    -0.8     0.3       0        
GRID        1028       0     1.6    -0.8     0.3       0        
GRID        1029       0     1.4    -0.8     0.3       0        
GRID        1030       0     2.2    -0.8     0.1       0        
GRID        1031       0      2.    -0.8     0.1       0        
GRID        1032       0     1.8    -0.8     0.1       0        
GRID        1033       0     1.6    -0.8     0.1       0        
GRID        1034       0     1.4    -0.8     0.1       0        
GRID        1035       0     2.2    -0.8    -0.1       0        
GRID        1036       0      2.    -0.8    -0.1       0        
GRID        1037       0     1.8    -0.8    -0.1       0        
GRID        1038       0     1.6    -0.8    -0.1       0        
GRID        1039       0     1.4    -0.8    -0.1       0        
GRID        1040       0     2.2    -0.8    -0.3       0        
GRID        1041       0      2.    -0.8    -0.3       0        
GRID        1042       0     1.8    -0.8    -0.3       0        
GRID        1043       0     1.6    -0.8    -0.3       0        
GRID        1044       0     1.4    -0.8    -0.3       0        
GRID        1045       0     2.2    -0.6     0.3       0        
GRID        1046       0      2.    -0.6     0.3       0        
GRID        1047       0     1.8    -0.6     0.3       0        
GRID        1048       0     1.6    -0.6     0.3       0        
GRID        1049       0     1.4    -0.6     0.3       0        
GRID        1050       0     2.2    -0.6     0.1       0        
GRID        1051       0      2.    -0.6     0.1       0        
GRID        1052       0     1.8    -0.6     0.1       0        
GRID        1053       0     1.6    -0.6     0.1       0        
GRID        1054       0     1.4    -0.6     0.1       0        
GRID        1055       0     2.2    -0.6    -0.1       0        
GRID        1056       0      2.    -0.6    -0.1       0        
GRID        1057       0     1.8    -0.6    -0.1       0        
GRID        1058       0     1.6    -0.6    -0.1       0        
GRID        1059       0     1.4    -0.6    -0.1       0        
GRID        1060       0     2.2    -0.6    -0.3       0        
GRID        1061       0      2.    -0.6    -0.3       0        
GRID        1062       0     1.8    -0.6    -0.3       0        
GRID        1063       0     1.6    -0.6    -0.3       0        
GRID        1064       0     1.4    -0.6    -0.3       0        
GRID        1065       0     2.2    -0.4     0.3       0        
GRID        1066       0      2.    -0.4     0.3       0        
GRID        1067       0     1.8    -0.4     0.3       0        
GRID        1068       0     1.6    -0.4     0.3       0        
GRID        1069       0     1.4    -0.4     0.3       0        
GRID        1070       0     2.2    -0.4     0.1       0        
GRID        1071       0      2.    -0.4     0.1       0        
GRID        1072       0     1.8    -0.4     0.1       0        
GRID        1073       0     1.6    -0.4     0.1       0        
GRID        1074       0     1.4    -0.4     0.1       0        
GRID        1075       0     2.2    -0.4    -0.1       0        
GRID        1076       0      2.    -0.4    -0.1       0        
GRID        1077       0     1.8    -0.4    -0.1       0        
GRID        1078       0     1.6    -0.4    -0.1       0        
GRID        1079       0     1.4    -0.4    -0.1       0        
GRID        1080       0     2.2    -0.4    -0.3       0        
GRID        1081       0      2.    -0.4    -0.3       0        
GRID        1082       0     1.8    -0.4    -0.3       0        
GRID        1083       0     1.6    -0.4    -0.3       0        
GRID        1084       0     1.4    -0.4    -0.3       0        
GRID        1085       0     2.2    -0.2     0.3       0        
GRID        1086       0      2.    -0.2     0.3       0        
GRID        1087       0     1.8    -0.2     0.3       0        
GRID        1088       0     1.6    -0.2     0.3       0        
GRID        1089       0     1.4    -0.2     0.3       0        
GRID        1090       0     2.2    -0.2     0.1       0        
GRID        1091       0      2.    -0.2     0.1       0        
GRID        1092       0     1.8    -0.2     0.1       0        
GRID        1093       0     1.6    -0.2     0.1       0        
GRID        1094       0     1.4    -0.2     0.1       0        
GRID        1095       0     2.2    -0.2    -0.1       0        
GRID        1096       0      2.    -0.2    -0.1       0        
GRID        1097       0     1.8    -0.2    -0.1       0        
GRID        1098       0     1.6    -0.2    -0.1       0        
GRID        1099       0     1.4    -0.2    -0.1       0        
GRID        1100       0     2.2    -0.2    -0.3       0        
GRID        1101       0      2.    -0.2    -0.3       0        
GRID        1102       0     1.8    -0.2    -0.3       0        
GRID        1103       0     1.6    -0.2    -0.3       0        
GRID        1104       0     1.4    -0.2    -0.3       0        
CHEXA          1       2       1      28       9       2      81      45+EL    1
+EL    1     705     102                                                        
CHEXA          2       2       2       9      10      80     102     705+EL    2
+EL    2     706     103                                                        
CHEXA          3       2      80      10      11       3     103     706+EL    3
+EL    3     707     104                                                        
CHEXA          4       2       3      11      12       4     104     707+EL    4
+EL    4     708     105                                                        
CHEXA          5       2       4      12     123     612     105     708+EL    5
+EL    5     124     120                                                        
CHEXA          6       2      28      29      13       9      45      40+EL    6
+EL    6     709     705                                                        
CHEXA          7       2       9      13      14      10     705     709+EL    7
+EL    7     710     706                                                        
CHEXA          8       2      10      14      15      11     706     710+EL    8
+EL    8     711     707                                                        
CHEXA          9       2      11      15      16      12     707     711+EL    9
+EL    9     712     708                                                        
CHEXA         10       2      12      16     122     123     708     712+EL    A
+EL    A     126     124                                                        
CHEXA         11       2      29       8      17      13      40      35+EL    B
+EL    B     713     709                                                        
CHEXA         12       2      13      17      18      14     709     713+EL    C
+EL    C     714     710                                                        
CHEXA         13       2      14      18      19      15     710     714+EL    D
+EL    D     715     711                                                        
CHEXA         14       2      15      19      20      16     711     715+EL    E
+EL    E     716     712                                                        
CHEXA         15       2      16      20     684     122     712     716+EL    F
+EL    F     622     126                                                        
CHEXA         16       2       8       7      21      17      35      30+EL    G
+EL    G     717     713                                                        
CHEXA         17       2      17      21      22      18     713     717+EL    H
+EL    H     718     714                                                        
CHEXA         18       2      18      22      23      19     714     718+EL    I
+EL    I     719     715                                                        
CHEXA         19       2      19      23      24      20     715     719+EL    J
+EL    J     720     716                                                        
CHEXA         20       2      20      24       5     684     716     720+EL    K
+EL    K     625     622                                                        
CHEXA         21       2       7      56      57      21      30      25+EL    L
+EL    L      75     717                                                        
CHEXA         22       2      21      57       6      22     717      75+EL    M
+EL    M      70     718                                                        
CHEXA         23       2      22       6      58      23     718      70+EL    N
+EL    N      65     719                                                        
CHEXA         24       2      23      58      59      24     719      65+EL    O
+EL    O      60     720                                                        
CHEXA         25       2      24      59     613       5     720      60+EL    P
+EL    P     631     625                                                        
CHEXA         26       2      81      45     705     102      82      46+EL    Q
+EL    Q     721      98                                                        
CHEXA         27       2     102     705     706     103      98     721+EL    R
+EL    R     722      99                                                        
CHEXA         28       2     103     706     707     104      99     722+EL    S
+EL    S     723     100                                                        
CHEXA         29       2     104     707     708     105     100     723+EL    T
+EL    T     724     101                                                        
CHEXA         30       2     105     708     124     120     101     724+EL    U
+EL    U     617     656                                                        
CHEXA         31       2      45      40     709     705      46      41+EL    V
+EL    V     725     721                                                        
CHEXA         32       2     705     709     710     706     721     725+EL    W
+EL    W     726     722                                                        
CHEXA         33       2     706     710     711     707     722     726+EL    X
+EL    X     727     723                                                        
CHEXA         34       2     707     711     712     708     723     727+EL    Y
+EL    Y     728     724                                                        
CHEXA         35       2     708     712     126     124     724     728+EL    Z
+EL    Z     127     617                                                        
CHEXA         36       2      40      35     713     709      41      36+EL   10
+EL   10     729     725                                                        
CHEXA         37       2     709     713     714     710     725     729+EL   11
+EL   11     730     726                                                        
CHEXA         38       2     710     714     715     711     726     730+EL   12
+EL   12     731     727                                                        
CHEXA         39       2     711     715     716     712     727     731+EL   13
+EL   13     732     728                                                        
CHEXA         40       2     712     716     622     126     728     732+EL   14
+EL   14     621     127                                                        
CHEXA         41       2      35      30     717     713      36      31+EL   15
+EL   15     733     729                                                        
CHEXA         42       2     713     717     718     714     729     733+EL   16
+EL   16     734     730                                                        
CHEXA         43       2     714     718     719     715     730     734+EL   17
+EL   17     735     731                                                        
CHEXA         44       2     715     719     720     716     731     735+EL   18
+EL   18     736     732                                                        
CHEXA         45       2     716     720     625     622     732     736+EL   19
+EL   19     131     621                                                        
CHEXA         46       2      30      25      75     717      31      55+EL   1A
+EL   1A      76     733                                                        
CHEXA         47       2     717      75      70     718     733      76+EL   1B
+EL   1B      71     734                                                        
CHEXA         48       2     718      70      65     719     734      71+EL   1C
+EL   1C      66     735                                                        
CHEXA         49       2     719      65      60     720     735      66+EL   1D
+EL   1D      61     736                                                        
CHEXA         50       2     720      60     631     625     736      61+EL   1E
+EL   1E      50     131                                                        
CHEXA         51       2      82      46     721      98      83      47+EL   1F
+EL   1F     737      94                                                        
CHEXA         52       2      98     721     722      99      94     737+EL   1G
+EL   1G     738      95                                                        
CHEXA         53       2      99     722     723     100      95     738+EL   1H
+EL   1H     739      96                                                        
CHEXA         54       2     100     723     724     101      96     739+EL   1I
+EL   1I     740      97                                                        
CHEXA         55       2     101     724     617     656      97     740+EL   1J
+EL   1J     616     657                                                        
CHEXA         56       2      46      41     725     721      47      42+EL   1K
+EL   1K     741     737                                                        
CHEXA         57       2     721     725     726     722     737     741+EL   1L
+EL   1L     742     738                                                        
CHEXA         58       2     722     726     727     723     738     742+EL   1M
+EL   1M     743     739                                                        
CHEXA         59       2     723     727     728     724     739     743+EL   1N
+EL   1N     744     740                                                        
CHEXA         60       2     724     728     127     617     740     744+EL   1O
+EL   1O     618     616                                                        
CHEXA         61       2      41      36     729     725      42      37+EL   1P
+EL   1P     745     741                                                        
CHEXA         62       2     725     729     730     726     741     745+EL   1Q
+EL   1Q     746     742                                                        
CHEXA         63       2     726     730     731     727     742     746+EL   1R
+EL   1R     747     743                                                        
CHEXA         64       2     727     731     732     728     743     747+EL   1S
+EL   1S     748     744                                                        
CHEXA         65       2     728     732     621     127     744     748+EL   1T
+EL   1T     130     618                                                        
CHEXA         66       2      36      31     733     729      37      32+EL   1U
+EL   1U     749     745                                                        
CHEXA         67       2     729     733     734     730     745     749+EL   1V
+EL   1V     750     746                                                        
CHEXA         68       2     730     734     735     731     746     750+EL   1W
+EL   1W     751     747                                                        
CHEXA         69       2     731     735     736     732     747     751+EL   1X
+EL   1X     752     748                                                        
CHEXA         70       2     732     736     131     621     748     752+EL   1Y
+EL   1Y     132     130                                                        
CHEXA         71       2      31      55      76     733      32      54+EL   1Z
+EL   1Z      77     749                                                        
CHEXA         72       2     733      76      71     734     749      77+EL   20
+EL   20      72     750                                                        
CHEXA         73       2     734      71      66     735     750      72+EL   21
+EL   21      67     751                                                        
CHEXA         74       2     735      66      61     736     751      67+EL   22
+EL   22      62     752                                                        
CHEXA         75       2     736      61      50     131     752      62+EL   23
+EL   23     614     132                                                        
CHEXA         76       2      83      47     737      94      84      48+EL   24
+EL   24     753      90                                                        
CHEXA         77       2      94     737     738      95      90     753+EL   25
+EL   25     754      91                                                        
CHEXA         78       2      95     738     739      96      91     754+EL   26
+EL   26     755      92                                                        
CHEXA         79       2      96     739     740      97      92     755+EL   27
+EL   27     756      93                                                        
CHEXA         80       2      97     740     616     657      93     756+EL   28
+EL   28     615     611                                                        
CHEXA         81       2      47      42     741     737      48      43+EL   29
+EL   29     757     753                                                        
CHEXA         82       2     737     741     742     738     753     757+EL   2A
+EL   2A     758     754                                                        
CHEXA         83       2     738     742     743     739     754     758+EL   2B
+EL   2B     759     755                                                        
CHEXA         84       2     739     743     744     740     755     759+EL   2C
+EL   2C     760     756                                                        
CHEXA         85       2     740     744     618     616     756     760+EL   2D
+EL   2D     128     615                                                        
CHEXA         86       2      42      37     745     741      43      38+EL   2E
+EL   2E     761     757                                                        
CHEXA         87       2     741     745     746     742     757     761+EL   2F
+EL   2F     762     758                                                        
CHEXA         88       2     742     746     747     743     758     762+EL   2G
+EL   2G     763     759                                                        
CHEXA         89       2     743     747     748     744     759     763+EL   2H
+EL   2H     764     760                                                        
CHEXA         90       2     744     748     130     618     760     764+EL   2I
+EL   2I     620     128                                                        
CHEXA         91       2      37      32     749     745      38      33+EL   2J
+EL   2J     765     761                                                        
CHEXA         92       2     745     749     750     746     761     765+EL   2K
+EL   2K     766     762                                                        
CHEXA         93       2     746     750     751     747     762     766+EL   2L
+EL   2L     767     763                                                        
CHEXA         94       2     747     751     752     748     763     767+EL   2M
+EL   2M     768     764                                                        
CHEXA         95       2     748     752     132     130     764     768+EL   2N
+EL   2N     624     620                                                        
CHEXA         96       2      32      54      77     749      33      53+EL   2O
+EL   2O      78     765                                                        
CHEXA         97       2     749      77      72     750     765      78+EL   2P
+EL   2P      73     766                                                        
CHEXA         98       2     750      72      67     751     766      73+EL   2Q
+EL   2Q      68     767                                                        
CHEXA         99       2     751      67      62     752     767      68+EL   2R
+EL   2R      63     768                                                        
CHEXA        100       2     752      62     614     132     768      63+EL   2S
+EL   2S     630     624                                                        
CHEXA        101       2      84      48     753      90      85      49+EL   2T
+EL   2T     769      86                                                        
CHEXA        102       2      90     753     754      91      86     769+EL   2U
+EL   2U     770      87                                                        
CHEXA        103       2      91     754     755      92      87     770+EL   2V
+EL   2V     771      88                                                        
CHEXA        104       2      92     755     756      93      88     771+EL   2W
+EL   2W     772      89                                                        
CHEXA        105       2      93     756     615     611      89     772+EL   2X
+EL   2X     125     658                                                        
CHEXA        106       2      48      43     757     753      49      44+EL   2Y
+EL   2Y     773     769                                                        
CHEXA        107       2     753     757     758     754     769     773+EL   2Z
+EL   2Z     774     770                                                        
CHEXA        108       2     754     758     759     755     770     774+EL   30
+EL   30     775     771                                                        
CHEXA        109       2     755     759     760     756     771     775+EL   31
+EL   31     776     772                                                        
CHEXA        110       2     756     760     128     615     772     776+EL   32
+EL   32     129     125                                                        
CHEXA        111       2      43      38     761     757      44      39+EL   33
+EL   33     777     773                                                        
CHEXA        112       2     757     761     762     758     773     777+EL   34
+EL   34     778     774                                                        
CHEXA        113       2     758     762     763     759     774     778+EL   35
+EL   35     779     775                                                        
CHEXA        114       2     759     763     764     760     775     779+EL   36
+EL   36     780     776                                                        
CHEXA        115       2     760     764     620     128     776     780+EL   37
+EL   37     619     129                                                        
CHEXA        116       2      38      33     765     761      39      34+EL   38
+EL   38     781     777                                                        
CHEXA        117       2     761     765     766     762     777     781+EL   39
+EL   39     782     778                                                        
CHEXA        118       2     762     766     767     763     778     782+EL   3A
+EL   3A     783     779                                                        
CHEXA        119       2     763     767     768     764     779     783+EL   3B
+EL   3B     784     780                                                        
CHEXA        120       2     764     768     624     620     780     784+EL   3C
+EL   3C     623     619                                                        
CHEXA        121       2      33      53      78     765      34      26+EL   3D
+EL   3D      79     781                                                        
CHEXA        122       2     765      78      73     766     781      79+EL   3E
+EL   3E      74     782                                                        
CHEXA        123       2     766      73      68     767     782      74+EL   3F
+EL   3F      69     783                                                        
CHEXA        124       2     767      68      63     768     783      69+EL   3G
+EL   3G      64     784                                                        
CHEXA        125       2     768      63     630     624     784      64+EL   3H
+EL   3H      51     623                                                        
CHEXA        126       2      85      49     769      86     469     109+EL   3I
+EL   3I     558     108                                                        
CHEXA        127       2      86     769     770      87     108     558+EL   3J
+EL   3J     113     107                                                        
CHEXA        128       2      87     770     771      88     107     113+EL   3K
+EL   3K     559     106                                                        
CHEXA        129       2      88     771     772      89     106     559+EL   3L
+EL   3L     560     528                                                        
CHEXA        130       2      89     772     125     658     528     560+EL   3M
+EL   3M     121     554                                                        
CHEXA        131       2      49      44     773     769     109     557+EL   3N
+EL   3N     116     558                                                        
CHEXA        132       2     769     773     774     770     558     116+EL   3O
+EL   3O     561     113                                                        
CHEXA        133       2     770     774     775     771     113     561+EL   3P
+EL   3P     115     559                                                        
CHEXA        134       2     771     775     776     772     559     115+EL   3Q
+EL   3Q     114     560                                                        
CHEXA        135       2     772     776     129     125     560     114+EL   3R
+EL   3R     571     121                                                        
CHEXA        136       2      44      39     777     773     557      27+EL   3S
+EL   3S     562     116                                                        
CHEXA        137       2     773     777     778     774     116     562+EL   3T
+EL   3T     117     561                                                        
CHEXA        138       2     774     778     779     775     561     117+EL   3U
+EL   3U     563     115                                                        
CHEXA        139       2     775     779     780     776     115     563+EL   3V
+EL   3V     564     114                                                        
CHEXA        140       2     776     780     619     129     114     564+EL   3W
+EL   3W     112     571                                                        
CHEXA        141       2      39      34     781     777      27     470+EL   3X
+EL   3X     565     562                                                        
CHEXA        142       2     777     781     782     778     562     565+EL   3Y
+EL   3Y     119     117                                                        
CHEXA        143       2     778     782     783     779     117     119+EL   3Z
+EL   3Z     566     563                                                        
CHEXA        144       2     779     783     784     780     563     566+EL   40
+EL   40     118     564                                                        
CHEXA        145       2     780     784     623     619     564     118+EL   41
+EL   41     555     112                                                        
CHEXA        146       2      34      26      79     781     470     110+EL   42
+EL   42      52     565                                                        
CHEXA        147       2     781      79      74     782     565      52+EL   43
+EL   43     502     119                                                        
CHEXA        148       2     782      74      69     783     119     502+EL   44
+EL   44     556     566                                                        
CHEXA        149       2     783      69      64     784     566     556+EL   45
+EL   45     111     118                                                        
CHEXA        150       2     784      64      51     623     118     111+EL   46
+EL   46     503     555                                                        
CHEXA        151       2     193     194     142     133     251     204+EL   47
+EL   47     785     274                                                        
CHEXA        152       2     133     142     143     134     274     785+EL   48
+EL   48     786     270                                                        
CHEXA        153       2     134     143     144     135     270     786+EL   49
+EL   49     787     266                                                        
CHEXA        154       2     135     144     145     136     266     787+EL   4A
+EL   4A     788     262                                                        
CHEXA        155       2     136     145     146     137     262     788+EL   4B
+EL   4B     789     258                                                        
CHEXA        156       2     137     146     590     589     258     789+EL   4C
+EL   4C     283     277                                                        
CHEXA        157       2     194     141     147     142     204     208+EL   4D
+EL   4D     790     785                                                        
CHEXA        158       2     142     147     148     143     785     790+EL   4E
+EL   4E     791     786                                                        
CHEXA        159       2     143     148     149     144     786     791+EL   4F
+EL   4F     792     787                                                        
CHEXA        160       2     144     149     150     145     787     792+EL   4G
+EL   4G     793     788                                                        
CHEXA        161       2     145     150     151     146     788     793+EL   4H
+EL   4H     794     789                                                        
CHEXA        162       2     146     151     276     590     789     794+EL   4I
+EL   4I     282     283                                                        
CHEXA        163       2     141     195     152     147     208     212+EL   4J
+EL   4J     795     790                                                        
CHEXA        164       2     147     152     153     148     790     795+EL   4K
+EL   4K     796     791                                                        
CHEXA        165       2     148     153     154     149     791     796+EL   4L
+EL   4L     797     792                                                        
CHEXA        166       2     149     154     155     150     792     797+EL   4M
+EL   4M     798     793                                                        
CHEXA        167       2     150     155     156     151     793     798+EL   4N
+EL   4N     799     794                                                        
CHEXA        168       2     151     156     591     276     794     799+EL   4O
+EL   4O     281     282                                                        
CHEXA        169       2     195     196     157     152     212     216+EL   4P
+EL   4P     800     795                                                        
CHEXA        170       2     152     157     158     153     795     800+EL   4Q
+EL   4Q     801     796                                                        
CHEXA        171       2     153     158     159     154     796     801+EL   4R
+EL   4R     802     797                                                        
CHEXA        172       2     154     159     160     155     797     802+EL   4S
+EL   4S     803     798                                                        
CHEXA        173       2     155     160     161     156     798     803+EL   4T
+EL   4T     804     799                                                        
CHEXA        174       2     156     161     592     591     799     804+EL   4U
+EL   4U     280     281                                                        
CHEXA        175       2     196     217     218     157     216     226+EL   4V
+EL   4V     227     800                                                        
CHEXA        176       2     157     218     140     158     800     227+EL   4W
+EL   4W     228     801                                                        
CHEXA        177       2     158     140     139     159     801     228+EL   4X
+EL   4X     229     802                                                        
CHEXA        178       2     159     139     138     160     802     229+EL   4Y
+EL   4Y     230     803                                                        
CHEXA        179       2     160     138     219     161     803     230+EL   4Z
+EL   4Z     231     804                                                        
CHEXA        180       2     161     219     275     592     804     231+EL   50
+EL   50     345     280                                                        
CHEXA        181       2     251     204     785     274     192     203+EL   51
+EL   51     805     273                                                        
CHEXA        182       2     274     785     786     270     273     805+EL   52
+EL   52     806     269                                                        
CHEXA        183       2     270     786     787     266     269     806+EL   53
+EL   53     807     265                                                        
CHEXA        184       2     266     787     788     262     265     807+EL   54
+EL   54     808     261                                                        
CHEXA        185       2     262     788     789     258     261     808+EL   55
+EL   55     809     257                                                        
CHEXA        186       2     258     789     283     277     257     809+EL   56
+EL   56     285     250                                                        
CHEXA        187       2     204     208     790     785     203     207+EL   57
+EL   57     810     805                                                        
CHEXA        188       2     785     790     791     786     805     810+EL   58
+EL   58     811     806                                                        
CHEXA        189       2     786     791     792     787     806     811+EL   59
+EL   59     812     807                                                        
CHEXA        190       2     787     792     793     788     807     812+EL   5A
+EL   5A     813     808                                                        
CHEXA        191       2     788     793     794     789     808     813+EL   5B
+EL   5B     814     809                                                        
CHEXA        192       2     789     794     282     283     809     814+EL   5C
+EL   5C     324     285                                                        
CHEXA        193       2     208     212     795     790     207     211+EL   5D
+EL   5D     815     810                                                        
CHEXA        194       2     790     795     796     791     810     815+EL   5E
+EL   5E     816     811                                                        
CHEXA        195       2     791     796     797     792     811     816+EL   5F
+EL   5F     817     812                                                        
CHEXA        196       2     792     797     798     793     812     817+EL   5G
+EL   5G     818     813                                                        
CHEXA        197       2     793     798     799     794     813     818+EL   5H
+EL   5H     819     814                                                        
CHEXA        198       2     794     799     281     282     814     819+EL   5I
+EL   5I     325     324                                                        
CHEXA        199       2     212     216     800     795     211     215+EL   5J
+EL   5J     820     815                                                        
CHEXA        200       2     795     800     801     796     815     820+EL   5K
+EL   5K     821     816                                                        
CHEXA        201       2     796     801     802     797     816     821+EL   5L
+EL   5L     822     817                                                        
CHEXA        202       2     797     802     803     798     817     822+EL   5M
+EL   5M     823     818                                                        
CHEXA        203       2     798     803     804     799     818     823+EL   5N
+EL   5N     824     819                                                        
CHEXA        204       2     799     804     280     281     819     824+EL   5O
+EL   5O     284     325                                                        
CHEXA        205       2     216     226     227     800     215     197+EL   5P
+EL   5P     232     820                                                        
CHEXA        206       2     800     227     228     801     820     232+EL   5Q
+EL   5Q     233     821                                                        
CHEXA        207       2     801     228     229     802     821     233+EL   5R
+EL   5R     234     822                                                        
CHEXA        208       2     802     229     230     803     822     234+EL   5S
+EL   5S     235     823                                                        
CHEXA        209       2     803     230     231     804     823     235+EL   5T
+EL   5T     236     824                                                        
CHEXA        210       2     804     231     345     280     824     236+EL   5U
+EL   5U     279     284                                                        
CHEXA        211       2     192     203     805     273     252     202+EL   5V
+EL   5V     825     272                                                        
CHEXA        212       2     273     805     806     269     272     825+EL   5W
+EL   5W     826     268                                                        
CHEXA        213       2     269     806     807     265     268     826+EL   5X
+EL   5X     827     264                                                        
CHEXA        214       2     265     807     808     261     264     827+EL   5Y
+EL   5Y     828     260                                                        
CHEXA        215       2     261     808     809     257     260     828+EL   5Z
+EL   5Z     829     256                                                        
CHEXA        216       2     257     809     285     250     256     829+EL   60
+EL   60     287     249                                                        
CHEXA        217       2     203     207     810     805     202     206+EL   61
+EL   61     830     825                                                        
CHEXA        218       2     805     810     811     806     825     830+EL   62
+EL   62     831     826                                                        
CHEXA        219       2     806     811     812     807     826     831+EL   63
+EL   63     832     827                                                        
CHEXA        220       2     807     812     813     808     827     832+EL   64
+EL   64     833     828                                                        
CHEXA        221       2     808     813     814     809     828     833+EL   65
+EL   65     834     829                                                        
CHEXA        222       2     809     814     324     285     829     834+EL   66
+EL   66     286     287                                                        
CHEXA        223       2     207     211     815     810     206     210+EL   67
+EL   67     835     830                                                        
CHEXA        224       2     810     815     816     811     830     835+EL   68
+EL   68     836     831                                                        
CHEXA        225       2     811     816     817     812     831     836+EL   69
+EL   69     837     832                                                        
CHEXA        226       2     812     817     818     813     832     837+EL   6A
+EL   6A     838     833                                                        
CHEXA        227       2     813     818     819     814     833     838+EL   6B
+EL   6B     839     834                                                        
CHEXA        228       2     814     819     325     324     834     839+EL   6C
+EL   6C     326     286                                                        
CHEXA        229       2     211     215     820     815     210     214+EL   6D
+EL   6D     840     835                                                        
CHEXA        230       2     815     820     821     816     835     840+EL   6E
+EL   6E     841     836                                                        
CHEXA        231       2     816     821     822     817     836     841+EL   6F
+EL   6F     842     837                                                        
CHEXA        232       2     817     822     823     818     837     842+EL   6G
+EL   6G     843     838                                                        
CHEXA        233       2     818     823     824     819     838     843+EL   6H
+EL   6H     844     839                                                        
CHEXA        234       2     819     824     284     325     839     844+EL   6I
+EL   6I     327     326                                                        
CHEXA        235       2     215     197     232     820     214     225+EL   6J
+EL   6J     237     840                                                        
CHEXA        236       2     820     232     233     821     840     237+EL   6K
+EL   6K     238     841                                                        
CHEXA        237       2     821     233     234     822     841     238+EL   6L
+EL   6L     239     842                                                        
CHEXA        238       2     822     234     235     823     842     239+EL   6M
+EL   6M     240     843                                                        
CHEXA        239       2     823     235     236     824     843     240+EL   6N
+EL   6N     241     844                                                        
CHEXA        240       2     824     236     279     284     844     241+EL   6O
+EL   6O     220     327                                                        
CHEXA        241       2     252     202     825     272     191     201+EL   6P
+EL   6P     845     271                                                        
CHEXA        242       2     272     825     826     268     271     845+EL   6Q
+EL   6Q     846     267                                                        
CHEXA        243       2     268     826     827     264     267     846+EL   6R
+EL   6R     847     263                                                        
CHEXA        244       2     264     827     828     260     263     847+EL   6S
+EL   6S     848     259                                                        
CHEXA        245       2     260     828     829     256     259     848+EL   6T
+EL   6T     849     255                                                        
CHEXA        246       2     256     829     287     249     255     849+EL   6U
+EL   6U     289     248                                                        
CHEXA        247       2     202     206     830     825     201     205+EL   6V
+EL   6V     850     845                                                        
CHEXA        248       2     825     830     831     826     845     850+EL   6W
+EL   6W     851     846                                                        
CHEXA        249       2     826     831     832     827     846     851+EL   6X
+EL   6X     852     847                                                        
CHEXA        250       2     827     832     833     828     847     852+EL   6Y
+EL   6Y     853     848                                                        
CHEXA        251       2     828     833     834     829     848     853+EL   6Z
+EL   6Z     854     849                                                        
CHEXA        252       2     829     834     286     287     849     854+EL   70
+EL   70     328     289                                                        
CHEXA        253       2     206     210     835     830     205     209+EL   71
+EL   71     855     850                                                        
CHEXA        254       2     830     835     836     831     850     855+EL   72
+EL   72     856     851                                                        
CHEXA        255       2     831     836     837     832     851     856+EL   73
+EL   73     857     852                                                        
CHEXA        256       2     832     837     838     833     852     857+EL   74
+EL   74     858     853                                                        
CHEXA        257       2     833     838     839     834     853     858+EL   75
+EL   75     859     854                                                        
CHEXA        258       2     834     839     326     286     854     859+EL   76
+EL   76     329     328                                                        
CHEXA        259       2     210     214     840     835     209     213+EL   77
+EL   77     860     855                                                        
CHEXA        260       2     835     840     841     836     855     860+EL   78
+EL   78     861     856                                                        
CHEXA        261       2     836     841     842     837     856     861+EL   79
+EL   79     862     857                                                        
CHEXA        262       2     837     842     843     838     857     862+EL   7A
+EL   7A     863     858                                                        
CHEXA        263       2     838     843     844     839     858     863+EL   7B
+EL   7B     864     859                                                        
CHEXA        264       2     839     844     327     326     859     864+EL   7C
+EL   7C     288     329                                                        
CHEXA        265       2     214     225     237     840     213     198+EL   7D
+EL   7D     242     860                                                        
CHEXA        266       2     840     237     238     841     860     242+EL   7E
+EL   7E     243     861                                                        
CHEXA        267       2     841     238     239     842     861     243+EL   7F
+EL   7F     244     862                                                        
CHEXA        268       2     842     239     240     843     862     244+EL   7G
+EL   7G     245     863                                                        
CHEXA        269       2     843     240     241     844     863     245+EL   7H
+EL   7H     246     864                                                        
CHEXA        270       2     844     241     220     327     864     246+EL   7I
+EL   7I     221     288                                                        
CHEXA        271       2     191     201     845     271     253     200+EL   7J
+EL   7J     186     167                                                        
CHEXA        272       2     271     845     846     267     167     186+EL   7K
+EL   7K     187     166                                                        
CHEXA        273       2     267     846     847     263     166     187+EL   7L
+EL   7L     188     165                                                        
CHEXA        274       2     263     847     848     259     165     188+EL   7M
+EL   7M     189     164                                                        
CHEXA        275       2     259     848     849     255     164     189+EL   7N
+EL   7N     190     254                                                        
CHEXA        276       2     255     849     289     248     254     190+EL   7O
+EL   7O     278     247                                                        
CHEXA        277       2     201     205     850     845     200     168+EL   7P
+EL   7P     181     186                                                        
CHEXA        278       2     845     850     851     846     186     181+EL   7Q
+EL   7Q     182     187                                                        
CHEXA        279       2     846     851     852     847     187     182+EL   7R
+EL   7R     183     188                                                        
CHEXA        280       2     847     852     853     848     188     183+EL   7S
+EL   7S     184     189                                                        
CHEXA        281       2     848     853     854     849     189     184+EL   7T
+EL   7T     185     190                                                        
CHEXA        282       2     849     854     328     289     190     185+EL   7U
+EL   7U     290     278                                                        
CHEXA        283       2     205     209     855     850     168     169+EL   7V
+EL   7V     176     181                                                        
CHEXA        284       2     850     855     856     851     181     176+EL   7W
+EL   7W     177     182                                                        
CHEXA        285       2     851     856     857     852     182     177+EL   7X
+EL   7X     178     183                                                        
CHEXA        286       2     852     857     858     853     183     178+EL   7Y
+EL   7Y     179     184                                                        
CHEXA        287       2     853     858     859     854     184     179+EL   7Z
+EL   7Z     180     185                                                        
CHEXA        288       2     854     859     329     328     185     180+EL   80
+EL   80     291     290                                                        
CHEXA        289       2     209     213     860     855     169     170+EL   81
+EL   81     171     176                                                        
CHEXA        290       2     855     860     861     856     176     171+EL   82
+EL   82     172     177                                                        
CHEXA        291       2     856     861     862     857     177     172+EL   83
+EL   83     173     178                                                        
CHEXA        292       2     857     862     863     858     178     173+EL   84
+EL   84     174     179                                                        
CHEXA        293       2     858     863     864     859     179     174+EL   85
+EL   85     175     180                                                        
CHEXA        294       2     859     864     288     329     180     175+EL   86
+EL   86     292     291                                                        
CHEXA        295       2     213     198     242     860     170     199+EL   87
+EL   87     224     171                                                        
CHEXA        296       2     860     242     243     861     171     224+EL   88
+EL   88     223     172                                                        
CHEXA        297       2     861     243     244     862     172     223+EL   89
+EL   89     222     173                                                        
CHEXA        298       2     862     244     245     863     173     222+EL   8A
+EL   8A     162     174                                                        
CHEXA        299       2     863     245     246     864     174     162+EL   8B
+EL   8B     163     175                                                        
CHEXA        300       2     864     246     221     288     175     163+EL   8C
+EL   8C     293     292                                                        
CHEXA        301       2     589     590     602     330     277     283+EL   8D
+EL   8D     865     391                                                        
CHEXA        302       2     330     602     601     331     391     865+EL   8E
+EL   8E     866     387                                                        
CHEXA        303       2     331     601     600     368     387     866+EL   8F
+EL   8F     867     383                                                        
CHEXA        304       2     368     600     599     635     383     867+EL   8G
+EL   8G     868     379                                                        
CHEXA        305       2     635     599     598     588     379     868+EL   8H
+EL   8H     869     375                                                        
CHEXA        306       2     588     598     332     567     375     869+EL   8I
+EL   8I     400     392                                                        
CHEXA        307       2     590     276     334     602     283     282+EL   8J
+EL   8J     870     865                                                        
CHEXA        308       2     602     334     605     601     865     870+EL   8K
+EL   8K     871     866                                                        
CHEXA        309       2     601     605     604     600     866     871+EL   8L
+EL   8L     872     867                                                        
CHEXA        310       2     600     604     603     599     867     872+EL   8M
+EL   8M     873     868                                                        
CHEXA        311       2     599     603     335     598     868     873+EL   8N
+EL   8N     874     869                                                        
CHEXA        312       2     598     335     597     332     869     874+EL   8O
+EL   8O     399     400                                                        
CHEXA        313       2     276     591     608     334     282     281+EL   8P
+EL   8P     875     870                                                        
CHEXA        314       2     334     608     336     605     870     875+EL   8Q
+EL   8Q     876     871                                                        
CHEXA        315       2     605     336     607     604     871     876+EL   8R
+EL   8R     877     872                                                        
CHEXA        316       2     604     607     606     603     872     877+EL   8S
+EL   8S     878     873                                                        
CHEXA        317       2     603     606     337     335     873     878+EL   8T
+EL   8T     879     874                                                        
CHEXA        318       2     335     337     596     597     874     879+EL   8U
+EL   8U     398     399                                                        
CHEXA        319       2     591     592     338     608     281     280+EL   8V
+EL   8V     880     875                                                        
CHEXA        320       2     608     338     610     336     875     880+EL   8W
+EL   8W     881     876                                                        
CHEXA        321       2     336     610     339     607     876     881+EL   8X
+EL   8X     882     877                                                        
CHEXA        322       2     607     339     340     606     877     882+EL   8Y
+EL   8Y     883     878                                                        
CHEXA        323       2     606     340     609     337     878     883+EL   8Z
+EL   8Z     884     879                                                        
CHEXA        324       2     337     609     568     596     879     884+EL   90
+EL   90     397     398                                                        
CHEXA        325       2     592     275     333     338     280     345+EL   91
+EL   91     346     880                                                        
CHEXA        326       2     338     333     593     610     880     346+EL   92
+EL   92     347     881                                                        
CHEXA        327       2     610     593     594     339     881     347+EL   93
+EL   93     348     882                                                        
CHEXA        328       2     339     594     662     340     882     348+EL   94
+EL   94     349     883                                                        
CHEXA        329       2     340     662     595     609     883     349+EL   95
+EL   95     350     884                                                        
CHEXA        330       2     609     595     661     568     884     350+EL   96
+EL   96     396     397                                                        
CHEXA        331       2     277     283     865     391     250     285+EL   97
+EL   97     885     390                                                        
CHEXA        332       2     391     865     866     387     390     885+EL   98
+EL   98     886     386                                                        
CHEXA        333       2     387     866     867     383     386     886+EL   99
+EL   99     887     382                                                        
CHEXA        334       2     383     867     868     379     382     887+EL   9A
+EL   9A     888     378                                                        
CHEXA        335       2     379     868     869     375     378     888+EL   9B
+EL   9B     889     374                                                        
CHEXA        336       2     375     869     400     392     374     889+EL   9C
+EL   9C     404     367                                                        
CHEXA        337       2     283     282     870     865     285     324+EL   9D
+EL   9D     890     885                                                        
CHEXA        338       2     865     870     871     866     885     890+EL   9E
+EL   9E     891     886                                                        
CHEXA        339       2     866     871     872     867     886     891+EL   9F
+EL   9F     892     887                                                        
CHEXA        340       2     867     872     873     868     887     892+EL   9G
+EL   9G     893     888                                                        
CHEXA        341       2     868     873     874     869     888     893+EL   9H
+EL   9H     894     889                                                        
CHEXA        342       2     869     874     399     400     889     894+EL   9I
+EL   9I     403     404                                                        
CHEXA        343       2     282     281     875     870     324     325+EL   9J
+EL   9J     895     890                                                        
CHEXA        344       2     870     875     876     871     890     895+EL   9K
+EL   9K     896     891                                                        
CHEXA        345       2     871     876     877     872     891     896+EL   9L
+EL   9L     897     892                                                        
CHEXA        346       2     872     877     878     873     892     897+EL   9M
+EL   9M     898     893                                                        
CHEXA        347       2     873     878     879     874     893     898+EL   9N
+EL   9N     899     894                                                        
CHEXA        348       2     874     879     398     399     894     899+EL   9O
+EL   9O     402     403                                                        
CHEXA        349       2     281     280     880     875     325     284+EL   9P
+EL   9P     900     895                                                        
CHEXA        350       2     875     880     881     876     895     900+EL   9Q
+EL   9Q     901     896                                                        
CHEXA        351       2     876     881     882     877     896     901+EL   9R
+EL   9R     902     897                                                        
CHEXA        352       2     877     882     883     878     897     902+EL   9S
+EL   9S     903     898                                                        
CHEXA        353       2     878     883     884     879     898     903+EL   9T
+EL   9T     904     899                                                        
CHEXA        354       2     879     884     397     398     899     904+EL   9U
+EL   9U     401     402                                                        
CHEXA        355       2     280     345     346     880     284     279+EL   9V
+EL   9V     351     900                                                        
CHEXA        356       2     880     346     347     881     900     351+EL   9W
+EL   9W     352     901                                                        
CHEXA        357       2     881     347     348     882     901     352+EL   9X
+EL   9X     353     902                                                        
CHEXA        358       2     882     348     349     883     902     353+EL   9Y
+EL   9Y     354     903                                                        
CHEXA        359       2     883     349     350     884     903     354+EL   9Z
+EL   9Z     355     904                                                        
CHEXA        360       2     884     350     396     397     904     355+EL   A0
+EL   A0     341     401                                                        
CHEXA        361       2     250     285     885     390     249     287+EL   A1
+EL   A1     905     389                                                        
CHEXA        362       2     390     885     886     386     389     905+EL   A2
+EL   A2     906     385                                                        
CHEXA        363       2     386     886     887     382     385     906+EL   A3
+EL   A3     907     381                                                        
CHEXA        364       2     382     887     888     378     381     907+EL   A4
+EL   A4     908     377                                                        
CHEXA        365       2     378     888     889     374     377     908+EL   A5
+EL   A5     909     373                                                        
CHEXA        366       2     374     889     404     367     373     909+EL   A6
+EL   A6     408     393                                                        
CHEXA        367       2     285     324     890     885     287     286+EL   A7
+EL   A7     910     905                                                        
CHEXA        368       2     885     890     891     886     905     910+EL   A8
+EL   A8     911     906                                                        
CHEXA        369       2     886     891     892     887     906     911+EL   A9
+EL   A9     912     907                                                        
CHEXA        370       2     887     892     893     888     907     912+EL   AA
+EL   AA     913     908                                                        
CHEXA        371       2     888     893     894     889     908     913+EL   AB
+EL   AB     914     909                                                        
CHEXA        372       2     889     894     403     404     909     914+EL   AC
+EL   AC     407     408                                                        
CHEXA        373       2     324     325     895     890     286     326+EL   AD
+EL   AD     915     910                                                        
CHEXA        374       2     890     895     896     891     910     915+EL   AE
+EL   AE     916     911                                                        
CHEXA        375       2     891     896     897     892     911     916+EL   AF
+EL   AF     917     912                                                        
CHEXA        376       2     892     897     898     893     912     917+EL   AG
+EL   AG     918     913                                                        
CHEXA        377       2     893     898     899     894     913     918+EL   AH
+EL   AH     919     914                                                        
CHEXA        378       2     894     899     402     403     914     919+EL   AI
+EL   AI     406     407                                                        
CHEXA        379       2     325     284     900     895     326     327+EL   AJ
+EL   AJ     920     915                                                        
CHEXA        380       2     895     900     901     896     915     920+EL   AK
+EL   AK     921     916                                                        
CHEXA        381       2     896     901     902     897     916     921+EL   AL
+EL   AL     922     917                                                        
CHEXA        382       2     897     902     903     898     917     922+EL   AM
+EL   AM     923     918                                                        
CHEXA        383       2     898     903     904     899     918     923+EL   AN
+EL   AN     924     919                                                        
CHEXA        384       2     899     904     401     402     919     924+EL   AO
+EL   AO     405     406                                                        
CHEXA        385       2     284     279     351     900     327     220+EL   AP
+EL   AP     356     920                                                        
CHEXA        386       2     900     351     352     901     920     356+EL   AQ
+EL   AQ     357     921                                                        
CHEXA        387       2     901     352     353     902     921     357+EL   AR
+EL   AR     358     922                                                        
CHEXA        388       2     902     353     354     903     922     358+EL   AS
+EL   AS     359     923                                                        
CHEXA        389       2     903     354     355     904     923     359+EL   AT
+EL   AT     360     924                                                        
CHEXA        390       2     904     355     341     401     924     360+EL   AU
+EL   AU     342     405                                                        
CHEXA        391       2     249     287     905     389     248     289+EL   AV
+EL   AV     925     388                                                        
CHEXA        392       2     389     905     906     385     388     925+EL   AW
+EL   AW     926     384                                                        
CHEXA        393       2     385     906     907     381     384     926+EL   AX
+EL   AX     927     380                                                        
CHEXA        394       2     381     907     908     377     380     927+EL   AY
+EL   AY     928     376                                                        
CHEXA        395       2     377     908     909     373     376     928+EL   AZ
+EL   AZ     929     372                                                        
CHEXA        396       2     373     909     408     393     372     929+EL   B0
+EL   B0     412     394                                                        
CHEXA        397       2     287     286     910     905     289     328+EL   B1
+EL   B1     930     925                                                        
CHEXA        398       2     905     910     911     906     925     930+EL   B2
+EL   B2     931     926                                                        
CHEXA        399       2     906     911     912     907     926     931+EL   B3
+EL   B3     932     927                                                        
CHEXA        400       2     907     912     913     908     927     932+EL   B4
+EL   B4     933     928                                                        
CHEXA        401       2     908     913     914     909     928     933+EL   B5
+EL   B5     934     929                                                        
CHEXA        402       2     909     914     407     408     929     934+EL   B6
+EL   B6     411     412                                                        
CHEXA        403       2     286     326     915     910     328     329+EL   B7
+EL   B7     935     930                                                        
CHEXA        404       2     910     915     916     911     930     935+EL   B8
+EL   B8     936     931                                                        
CHEXA        405       2     911     916     917     912     931     936+EL   B9
+EL   B9     937     932                                                        
CHEXA        406       2     912     917     918     913     932     937+EL   BA
+EL   BA     938     933                                                        
CHEXA        407       2     913     918     919     914     933     938+EL   BB
+EL   BB     939     934                                                        
CHEXA        408       2     914     919     406     407     934     939+EL   BC
+EL   BC     410     411                                                        
CHEXA        409       2     326     327     920     915     329     288+EL   BD
+EL   BD     940     935                                                        
CHEXA        410       2     915     920     921     916     935     940+EL   BE
+EL   BE     941     936                                                        
CHEXA        411       2     916     921     922     917     936     941+EL   BF
+EL   BF     942     937                                                        
CHEXA        412       2     917     922     923     918     937     942+EL   BG
+EL   BG     943     938                                                        
CHEXA        413       2     918     923     924     919     938     943+EL   BH
+EL   BH     944     939                                                        
CHEXA        414       2     919     924     405     406     939     944+EL   BI
+EL   BI     409     410                                                        
CHEXA        415       2     327     220     356     920     288     221+EL   BJ
+EL   BJ     361     940                                                        
CHEXA        416       2     920     356     357     921     940     361+EL   BK
+EL   BK     362     941                                                        
CHEXA        417       2     921     357     358     922     941     362+EL   BL
+EL   BL     363     942                                                        
CHEXA        418       2     922     358     359     923     942     363+EL   BM
+EL   BM     364     943                                                        
CHEXA        419       2     923     359     360     924     943     364+EL   BN
+EL   BN     365     944                                                        
CHEXA        420       2     924     360     342     405     944     365+EL   BO
+EL   BO     343     409                                                        
CHEXA        421       2     248     289     925     388     247     278+EL   BP
+EL   BP     304     369                                                        
CHEXA        422       2     388     925     926     384     369     304+EL   BQ
+EL   BQ     308     370                                                        
CHEXA        423       2     384     926     927     380     370     308+EL   BR
+EL   BR     312     371                                                        
CHEXA        424       2     380     927     928     376     371     312+EL   BS
+EL   BS     316     303                                                        
CHEXA        425       2     376     928     929     372     303     316+EL   BT
+EL   BT     320     302                                                        
CHEXA        426       2     372     929     412     394     302     320+EL   BU
+EL   BU     301     366                                                        
CHEXA        427       2     289     328     930     925     278     290+EL   BV
+EL   BV     305     304                                                        
CHEXA        428       2     925     930     931     926     304     305+EL   BW
+EL   BW     309     308                                                        
CHEXA        429       2     926     931     932     927     308     309+EL   BX
+EL   BX     313     312                                                        
CHEXA        430       2     927     932     933     928     312     313+EL   BY
+EL   BY     317     316                                                        
CHEXA        431       2     928     933     934     929     316     317+EL   BZ
+EL   BZ     321     320                                                        
CHEXA        432       2     929     934     411     412     320     321+EL   C0
+EL   C0     300     301                                                        
CHEXA        433       2     328     329     935     930     290     291+EL   C1
+EL   C1     306     305                                                        
CHEXA        434       2     930     935     936     931     305     306+EL   C2
+EL   C2     310     309                                                        
CHEXA        435       2     931     936     937     932     309     310+EL   C3
+EL   C3     314     313                                                        
CHEXA        436       2     932     937     938     933     313     314+EL   C4
+EL   C4     318     317                                                        
CHEXA        437       2     933     938     939     934     317     318+EL   C5
+EL   C5     322     321                                                        
CHEXA        438       2     934     939     410     411     321     322+EL   C6
+EL   C6     299     300                                                        
CHEXA        439       2     329     288     940     935     291     292+EL   C7
+EL   C7     307     306                                                        
CHEXA        440       2     935     940     941     936     306     307+EL   C8
+EL   C8     311     310                                                        
CHEXA        441       2     936     941     942     937     310     311+EL   C9
+EL   C9     315     314                                                        
CHEXA        442       2     937     942     943     938     314     315+EL   CA
+EL   CA     319     318                                                        
CHEXA        443       2     938     943     944     939     318     319+EL   CB
+EL   CB     323     322                                                        
CHEXA        444       2     939     944     409     410     322     323+EL   CC
+EL   CC     395     299                                                        
CHEXA        445       2     288     221     361     940     292     293+EL   CD
+EL   CD     344     307                                                        
CHEXA        446       2     940     361     362     941     307     344+EL   CE
+EL   CE     294     311                                                        
CHEXA        447       2     941     362     363     942     311     294+EL   CF
+EL   CF     295     315                                                        
CHEXA        448       2     942     363     364     943     315     295+EL   CG
+EL   CG     296     319                                                        
CHEXA        449       2     943     364     365     944     319     296+EL   CH
+EL   CH     297     323                                                        
CHEXA        450       2     944     365     343     409     323     297+EL   CI
+EL   CI     298     395                                                        
CHEXA        451       2     469     529     479     109     108     537+EL   CJ
+EL   CJ     945     558                                                        
CHEXA        452       2     109     479     480     557     558     945+EL   CK
+EL   CK     946     116                                                        
CHEXA        453       2     557     480     481      27     116     946+EL   CL
+EL   CL     947     562                                                        
CHEXA        454       2      27     481     482     470     562     947+EL   CM
+EL   CM     948     565                                                        
CHEXA        455       2     470     482     471     110     565     948+EL   CN
+EL   CN     512      52                                                        
CHEXA        456       2     529     478     483     479     537     541+EL   CO
+EL   CO     949     945                                                        
CHEXA        457       2     479     483     484     480     945     949+EL   CP
+EL   CP     950     946                                                        
CHEXA        458       2     480     484     485     481     946     950+EL   CQ
+EL   CQ     951     947                                                        
CHEXA        459       2     481     485     486     482     947     951+EL   CR
+EL   CR     952     948                                                        
CHEXA        460       2     482     486     472     471     948     952+EL   CS
+EL   CS     511     512                                                        
CHEXA        461       2     478     530     487     483     541     545+EL   CT
+EL   CT     953     949                                                        
CHEXA        462       2     483     487     488     484     949     953+EL   CU
+EL   CU     954     950                                                        
CHEXA        463       2     484     488     489     485     950     954+EL   CV
+EL   CV     955     951                                                        
CHEXA        464       2     485     489     490     486     951     955+EL   CW
+EL   CW     956     952                                                        
CHEXA        465       2     486     490     501     472     952     956+EL   CX
+EL   CX     510     511                                                        
CHEXA        466       2     530     477     491     487     545     549+EL   CY
+EL   CY     957     953                                                        
CHEXA        467       2     487     491     492     488     953     957+EL   CZ
+EL   CZ     958     954                                                        
CHEXA        468       2     488     492     493     489     954     958+EL   D0
+EL   D0     959     955                                                        
CHEXA        469       2     489     493     494     490     955     959+EL   D1
+EL   D1     960     956                                                        
CHEXA        470       2     490     494     500     501     956     960+EL   D2
+EL   D2     509     510                                                        
CHEXA        471       2     477     476     495     491     549     553+EL   D3
+EL   D3     961     957                                                        
CHEXA        472       2     491     495     496     492     957     961+EL   D4
+EL   D4     962     958                                                        
CHEXA        473       2     492     496     497     493     958     962+EL   D5
+EL   D5     963     959                                                        
CHEXA        474       2     493     497     498     494     959     963+EL   D6
+EL   D6     964     960                                                        
CHEXA        475       2     494     498     499     500     960     964+EL   D7
+EL   D7     508     509                                                        
CHEXA        476       2     476     531     475     495     553     422+EL   D8
+EL   D8     436     961                                                        
CHEXA        477       2     495     475     423     496     961     436+EL   D9
+EL   D9     432     962                                                        
CHEXA        478       2     496     423     474     497     962     432+EL   DA
+EL   DA     428     963                                                        
CHEXA        479       2     497     474     473     498     963     428+EL   DB
+EL   DB     424     964                                                        
CHEXA        480       2     498     473     413     499     964     424+EL   DC
+EL   DC     507     508                                                        
CHEXA        481       2     108     537     945     558     107     536+EL   DD
+EL   DD     965     113                                                        
CHEXA        482       2     558     945     946     116     113     965+EL   DE
+EL   DE     966     561                                                        
CHEXA        483       2     116     946     947     562     561     966+EL   DF
+EL   DF     967     117                                                        
CHEXA        484       2     562     947     948     565     117     967+EL   DG
+EL   DG     968     119                                                        
CHEXA        485       2     565     948     512      52     119     968+EL   DH
+EL   DH     517     502                                                        
CHEXA        486       2     537     541     949     945     536     540+EL   DI
+EL   DI     969     965                                                        
CHEXA        487       2     945     949     950     946     965     969+EL   DJ
+EL   DJ     970     966                                                        
CHEXA        488       2     946     950     951     947     966     970+EL   DK
+EL   DK     971     967                                                        
CHEXA        489       2     947     951     952     948     967     971+EL   DL
+EL   DL     972     968                                                        
CHEXA        490       2     948     952     511     512     968     972+EL   DM
+EL   DM     516     517                                                        
CHEXA        491       2     541     545     953     949     540     544+EL   DN
+EL   DN     973     969                                                        
CHEXA        492       2     949     953     954     950     969     973+EL   DO
+EL   DO     974     970                                                        
CHEXA        493       2     950     954     955     951     970     974+EL   DP
+EL   DP     975     971                                                        
CHEXA        494       2     951     955     956     952     971     975+EL   DQ
+EL   DQ     976     972                                                        
CHEXA        495       2     952     956     510     511     972     976+EL   DR
+EL   DR     515     516                                                        
CHEXA        496       2     545     549     957     953     544     548+EL   DS
+EL   DS     977     973                                                        
CHEXA        497       2     953     957     958     954     973     977+EL   DT
+EL   DT     978     974                                                        
CHEXA        498       2     954     958     959     955     974     978+EL   DU
+EL   DU     979     975                                                        
CHEXA        499       2     955     959     960     956     975     979+EL   DV
+EL   DV     980     976                                                        
CHEXA        500       2     956     960     509     510     976     980+EL   DW
+EL   DW     514     515                                                        
CHEXA        501       2     549     553     961     957     548     552+EL   DX
+EL   DX     981     977                                                        
CHEXA        502       2     957     961     962     958     977     981+EL   DY
+EL   DY     982     978                                                        
CHEXA        503       2     958     962     963     959     978     982+EL   DZ
+EL   DZ     983     979                                                        
CHEXA        504       2     959     963     964     960     979     983+EL   E0
+EL   E0     984     980                                                        
CHEXA        505       2     960     964     508     509     980     984+EL   E1
+EL   E1     513     514                                                        
CHEXA        506       2     553     422     436     961     552     532+EL   E2
+EL   E2     437     981                                                        
CHEXA        507       2     961     436     432     962     981     437+EL   E3
+EL   E3     433     982                                                        
CHEXA        508       2     962     432     428     963     982     433+EL   E4
+EL   E4     429     983                                                        
CHEXA        509       2     963     428     424     964     983     429+EL   E5
+EL   E5     425     984                                                        
CHEXA        510       2     964     424     507     508     984     425+EL   E6
+EL   E6     414     513                                                        
CHEXA        511       2     107     536     965     113     106     535+EL   E7
+EL   E7     985     559                                                        
CHEXA        512       2     113     965     966     561     559     985+EL   E8
+EL   E8     986     115                                                        
CHEXA        513       2     561     966     967     117     115     986+EL   E9
+EL   E9     987     563                                                        
CHEXA        514       2     117     967     968     119     563     987+EL   EA
+EL   EA     988     566                                                        
CHEXA        515       2     119     968     517     502     566     988+EL   EB
+EL   EB     522     556                                                        
CHEXA        516       2     536     540     969     965     535     539+EL   EC
+EL   EC     989     985                                                        
CHEXA        517       2     965     969     970     966     985     989+EL   ED
+EL   ED     990     986                                                        
CHEXA        518       2     966     970     971     967     986     990+EL   EE
+EL   EE     991     987                                                        
CHEXA        519       2     967     971     972     968     987     991+EL   EF
+EL   EF     992     988                                                        
CHEXA        520       2     968     972     516     517     988     992+EL   EG
+EL   EG     521     522                                                        
CHEXA        521       2     540     544     973     969     539     543+EL   EH
+EL   EH     993     989                                                        
CHEXA        522       2     969     973     974     970     989     993+EL   EI
+EL   EI     994     990                                                        
CHEXA        523       2     970     974     975     971     990     994+EL   EJ
+EL   EJ     995     991                                                        
CHEXA        524       2     971     975     976     972     991     995+EL   EK
+EL   EK     996     992                                                        
CHEXA        525       2     972     976     515     516     992     996+EL   EL
+EL   EL     520     521                                                        
CHEXA        526       2     544     548     977     973     543     547+EL   EM
+EL   EM     997     993                                                        
CHEXA        527       2     973     977     978     974     993     997+EL   EN
+EL   EN     998     994                                                        
CHEXA        528       2     974     978     979     975     994     998+EL   EO
+EL   EO     999     995                                                        
CHEXA        529       2     975     979     980     976     995     999+EL   EP
+EL   EP    1000     996                                                        
CHEXA        530       2     976     980     514     515     996    1000+EL   EQ
+EL   EQ     519     520                                                        
CHEXA        531       2     548     552     981     977     547     551+EL   ER
+EL   ER    1001     997                                                        
CHEXA        532       2     977     981     982     978     997    1001+EL   ES
+EL   ES    1002     998                                                        
CHEXA        533       2     978     982     983     979     998    1002+EL   ET
+EL   ET    1003     999                                                        
CHEXA        534       2     979     983     984     980     999    1003+EL   EU
+EL   EU    1004    1000                                                        
CHEXA        535       2     980     984     513     514    1000    1004+EL   EV
+EL   EV     518     519                                                        
CHEXA        536       2     552     532     437     981     551     421+EL   EW
+EL   EW     438    1001                                                        
CHEXA        537       2     981     437     433     982    1001     438+EL   EX
+EL   EX     434    1002                                                        
CHEXA        538       2     982     433     429     983    1002     434+EL   EY
+EL   EY     430    1003                                                        
CHEXA        539       2     983     429     425     984    1003     430+EL   EZ
+EL   EZ     426    1004                                                        
CHEXA        540       2     984     425     414     513    1004     426+EL   F0
+EL   F0     415     518                                                        
CHEXA        541       2     106     535     985     559     528     534+EL   F1
+EL   F1    1005     560                                                        
CHEXA        542       2     559     985     986     115     560    1005+EL   F2
+EL   F2    1006     114                                                        
CHEXA        543       2     115     986     987     563     114    1006+EL   F3
+EL   F3    1007     564                                                        
CHEXA        544       2     563     987     988     566     564    1007+EL   F4
+EL   F4    1008     118                                                        
CHEXA        545       2     566     988     522     556     118    1008+EL   F5
+EL   F5     527     111                                                        
CHEXA        546       2     535     539     989     985     534     538+EL   F6
+EL   F6    1009    1005                                                        
CHEXA        547       2     985     989     990     986    1005    1009+EL   F7
+EL   F7    1010    1006                                                        
CHEXA        548       2     986     990     991     987    1006    1010+EL   F8
+EL   F8    1011    1007                                                        
CHEXA        549       2     987     991     992     988    1007    1011+EL   F9
+EL   F9    1012    1008                                                        
CHEXA        550       2     988     992     521     522    1008    1012+EL   FA
+EL   FA     526     527                                                        
CHEXA        551       2     539     543     993     989     538     542+EL   FB
+EL   FB    1013    1009                                                        
CHEXA        552       2     989     993     994     990    1009    1013+EL   FC
+EL   FC    1014    1010                                                        
CHEXA        553       2     990     994     995     991    1010    1014+EL   FD
+EL   FD    1015    1011                                                        
CHEXA        554       2     991     995     996     992    1011    1015+EL   FE
+EL   FE    1016    1012                                                        
CHEXA        555       2     992     996     520     521    1012    1016+EL   FF
+EL   FF     525     526                                                        
CHEXA        556       2     543     547     997     993     542     546+EL   FG
+EL   FG    1017    1013                                                        
CHEXA        557       2     993     997     998     994    1013    1017+EL   FH
+EL   FH    1018    1014                                                        
CHEXA        558       2     994     998     999     995    1014    1018+EL   FI
+EL   FI    1019    1015                                                        
CHEXA        559       2     995     999    1000     996    1015    1019+EL   FJ
+EL   FJ    1020    1016                                                        
CHEXA        560       2     996    1000     519     520    1016    1020+EL   FK
+EL   FK     524     525                                                        
CHEXA        561       2     547     551    1001     997     546     550+EL   FL
+EL   FL    1021    1017                                                        
CHEXA        562       2     997    1001    1002     998    1017    1021+EL   FM
+EL   FM    1022    1018                                                        
CHEXA        563       2     998    1002    1003     999    1018    1022+EL   FN
+EL   FN    1023    1019                                                        
CHEXA        564       2     999    1003    1004    1000    1019    1023+EL   FO
+EL   FO    1024    1020                                                        
CHEXA        565       2    1000    1004     518     519    1020    1024+EL   FP
+EL   FP     523     524                                                        
CHEXA        566       2     551     421     438    1001     550     420+EL   FQ
+EL   FQ     439    1021                                                        
CHEXA        567       2    1001     438     434    1002    1021     439+EL   FR
+EL   FR     435    1022                                                        
CHEXA        568       2    1002     434     430    1003    1022     435+EL   FS
+EL   FS     431    1023                                                        
CHEXA        569       2    1003     430     426    1004    1023     431+EL   FT
+EL   FT     427    1024                                                        
CHEXA        570       2    1004     426     415     518    1024     427+EL   FU
+EL   FU     416     523                                                        
CHEXA        571       2     528     534    1005     560     554     442+EL   FV
+EL   FV     468     121                                                        
CHEXA        572       2     560    1005    1006     114     121     468+EL   FW
+EL   FW     463     571                                                        
CHEXA        573       2     114    1006    1007     564     571     463+EL   FX
+EL   FX     458     112                                                        
CHEXA        574       2     564    1007    1008     118     112     458+EL   FY
+EL   FY     453     555                                                        
CHEXA        575       2     118    1008     527     111     555     453+EL   FZ
+EL   FZ     441     503                                                        
CHEXA        576       2     534     538    1009    1005     442     443+EL   G0
+EL   G0     467     468                                                        
CHEXA        577       2    1005    1009    1010    1006     468     467+EL   G1
+EL   G1     462     463                                                        
CHEXA        578       2    1006    1010    1011    1007     463     462+EL   G2
+EL   G2     457     458                                                        
CHEXA        579       2    1007    1011    1012    1008     458     457+EL   G3
+EL   G3     452     453                                                        
CHEXA        580       2    1008    1012     526     527     453     452+EL   G4
+EL   G4     504     441                                                        
CHEXA        581       2     538     542    1013    1009     443     444+EL   G5
+EL   G5     466     467                                                        
CHEXA        582       2    1009    1013    1014    1010     467     466+EL   G6
+EL   G6     461     462                                                        
CHEXA        583       2    1010    1014    1015    1011     462     461+EL   G7
+EL   G7     456     457                                                        
CHEXA        584       2    1011    1015    1016    1012     457     456+EL   G8
+EL   G8     451     452                                                        
CHEXA        585       2    1012    1016     525     526     452     451+EL   G9
+EL   G9     440     504                                                        
CHEXA        586       2     542     546    1017    1013     444     445+EL   GA
+EL   GA     465     466                                                        
CHEXA        587       2    1013    1017    1018    1014     466     465+EL   GB
+EL   GB     460     461                                                        
CHEXA        588       2    1014    1018    1019    1015     461     460+EL   GC
+EL   GC     455     456                                                        
CHEXA        589       2    1015    1019    1020    1016     456     455+EL   GD
+EL   GD     450     451                                                        
CHEXA        590       2    1016    1020     524     525     451     450+EL   GE
+EL   GE     505     440                                                        
CHEXA        591       2     546     550    1021    1017     445     533+EL   GF
+EL   GF     464     465                                                        
CHEXA        592       2    1017    1021    1022    1018     465     464+EL   GG
+EL   GG     459     460                                                        
CHEXA        593       2    1018    1022    1023    1019     460     459+EL   GH
+EL   GH     454     455                                                        
CHEXA        594       2    1019    1023    1024    1020     455     454+EL   GI
+EL   GI     449     450                                                        
CHEXA        595       2    1020    1024     523     524     450     449+EL   GJ
+EL   GJ     506     505                                                        
CHEXA        596       2     550     420     439    1021     533     419+EL   GK
+EL   GK     446     464                                                        
CHEXA        597       2    1021     439     435    1022     464     446+EL   GL
+EL   GL     447     459                                                        
CHEXA        598       2    1022     435     431    1023     459     447+EL   GM
+EL   GM     418     454                                                        
CHEXA        599       2    1023     431     427    1024     454     418+EL   GN
+EL   GN     448     449                                                        
CHEXA        600       2    1024     427     416     523     449     448+EL   GO
+EL   GO     417     506                                                        
CHEXA        601       2     554     121     125     658     570     584+EL   GP
+EL   GP    1025     668                                                        
CHEXA        602       2     658     125     615     611     668    1025+EL   GQ
+EL   GQ    1026     667                                                        
CHEXA        603       2     611     615     616     657     667    1026+EL   GR
+EL   GR    1027     666                                                        
CHEXA        604       2     657     616     617     656     666    1027+EL   GS
+EL   GS    1028     665                                                        
CHEXA        605       2     656     617     124     120     665    1028+EL   GT
+EL   GT    1029     664                                                        
CHEXA        606       2     120     124     123     612     664    1029+EL   GU
+EL   GU     692     685                                                        
CHEXA        607       2     121     571     129     125     584     580+EL   GV
+EL   GV    1030    1025                                                        
CHEXA        608       2     125     129     128     615    1025    1030+EL   GW
+EL   GW    1031    1026                                                        
CHEXA        609       2     615     128     618     616    1026    1031+EL   GX
+EL   GX    1032    1027                                                        
CHEXA        610       2     616     618     127     617    1027    1032+EL   GY
+EL   GY    1033    1028                                                        
CHEXA        611       2     617     127     126     124    1028    1033+EL   GZ
+EL   GZ    1034    1029                                                        
CHEXA        612       2     124     126     122     123    1029    1034+EL   H0
+EL   H0     691     692                                                        
CHEXA        613       2     571     112     619     129     580     576+EL   H1
+EL   H1    1035    1030                                                        
CHEXA        614       2     129     619     620     128    1030    1035+EL   H2
+EL   H2    1036    1031                                                        
CHEXA        615       2     128     620     130     618    1031    1036+EL   H3
+EL   H3    1037    1032                                                        
CHEXA        616       2     618     130     621     127    1032    1037+EL   H4
+EL   H4    1038    1033                                                        
CHEXA        617       2     127     621     622     126    1033    1038+EL   H5
+EL   H5    1039    1034                                                        
CHEXA        618       2     126     622     684     122    1034    1039+EL   H6
+EL   H6     690     691                                                        
CHEXA        619       2     112     555     623     619     576     572+EL   H7
+EL   H7    1040    1035                                                        
CHEXA        620       2     619     623     624     620    1035    1040+EL   H8
+EL   H8    1041    1036                                                        
CHEXA        621       2     620     624     132     130    1036    1041+EL   H9
+EL   H9    1042    1037                                                        
CHEXA        622       2     130     132     131     621    1037    1042+EL   HA
+EL   HA    1043    1038                                                        
CHEXA        623       2     621     131     625     622    1038    1043+EL   HB
+EL   HB    1044    1039                                                        
CHEXA        624       2     622     625       5     684    1039    1044+EL   HC
+EL   HC     689     690                                                        
CHEXA        625       2     555     503      51     623     572     629+EL   HD
+EL   HD     639    1040                                                        
CHEXA        626       2     623      51     630     624    1040     639+EL   HE
+EL   HE     643    1041                                                        
CHEXA        627       2     624     630     614     132    1041     643+EL   HF
+EL   HF     647    1042                                                        
CHEXA        628       2     132     614      50     131    1042     647+EL   HG
+EL   HG     651    1043                                                        
CHEXA        629       2     131      50     631     625    1043     651+EL   HH
+EL   HH     655    1044                                                        
CHEXA        630       2     625     631     613       5    1044     655+EL   HI
+EL   HI     632     689                                                        
CHEXA        631       2     570     584    1025     668     659     585+EL   HJ
+EL   HJ    1045     673                                                        
CHEXA        632       2     668    1025    1026     667     673    1045+EL   HK
+EL   HK    1046     672                                                        
CHEXA        633       2     667    1026    1027     666     672    1046+EL   HL
+EL   HL    1047     671                                                        
CHEXA        634       2     666    1027    1028     665     671    1047+EL   HM
+EL   HM    1048     670                                                        
CHEXA        635       2     665    1028    1029     664     670    1048+EL   HN
+EL   HN    1049     669                                                        
CHEXA        636       2     664    1029     692     685     669    1049+EL   HO
+EL   HO     696     663                                                        
CHEXA        637       2     584     580    1030    1025     585     581+EL   HP
+EL   HP    1050    1045                                                        
CHEXA        638       2    1025    1030    1031    1026    1045    1050+EL   HQ
+EL   HQ    1051    1046                                                        
CHEXA        639       2    1026    1031    1032    1027    1046    1051+EL   HR
+EL   HR    1052    1047                                                        
CHEXA        640       2    1027    1032    1033    1028    1047    1052+EL   HS
+EL   HS    1053    1048                                                        
CHEXA        641       2    1028    1033    1034    1029    1048    1053+EL   HT
+EL   HT    1054    1049                                                        
CHEXA        642       2    1029    1034     691     692    1049    1054+EL   HU
+EL   HU     695     696                                                        
CHEXA        643       2     580     576    1035    1030     581     577+EL   HV
+EL   HV    1055    1050                                                        
CHEXA        644       2    1030    1035    1036    1031    1050    1055+EL   HW
+EL   HW    1056    1051                                                        
CHEXA        645       2    1031    1036    1037    1032    1051    1056+EL   HX
+EL   HX    1057    1052                                                        
CHEXA        646       2    1032    1037    1038    1033    1052    1057+EL   HY
+EL   HY    1058    1053                                                        
CHEXA        647       2    1033    1038    1039    1034    1053    1058+EL   HZ
+EL   HZ    1059    1054                                                        
CHEXA        648       2    1034    1039     690     691    1054    1059+EL   I0
+EL   I0     694     695                                                        
CHEXA        649       2     576     572    1040    1035     577     573+EL   I1
+EL   I1    1060    1055                                                        
CHEXA        650       2    1035    1040    1041    1036    1055    1060+EL   I2
+EL   I2    1061    1056                                                        
CHEXA        651       2    1036    1041    1042    1037    1056    1061+EL   I3
+EL   I3    1062    1057                                                        
CHEXA        652       2    1037    1042    1043    1038    1057    1062+EL   I4
+EL   I4    1063    1058                                                        
CHEXA        653       2    1038    1043    1044    1039    1058    1063+EL   I5
+EL   I5    1064    1059                                                        
CHEXA        654       2    1039    1044     689     690    1059    1064+EL   I6
+EL   I6     693     694                                                        
CHEXA        655       2     572     629     639    1040     573     628+EL   I7
+EL   I7     638    1060                                                        
CHEXA        656       2    1040     639     643    1041    1060     638+EL   I8
+EL   I8     642    1061                                                        
CHEXA        657       2    1041     643     647    1042    1061     642+EL   I9
+EL   I9     646    1062                                                        
CHEXA        658       2    1042     647     651    1043    1062     646+EL   IA
+EL   IA     650    1063                                                        
CHEXA        659       2    1043     651     655    1044    1063     650+EL   IB
+EL   IB     654    1064                                                        
CHEXA        660       2    1044     655     632     689    1064     654+EL   IC
+EL   IC     633     693                                                        
CHEXA        661       2     659     585    1045     673     569     586+EL   ID
+EL   ID    1065     678                                                        
CHEXA        662       2     673    1045    1046     672     678    1065+EL   IE
+EL   IE    1066     677                                                        
CHEXA        663       2     672    1046    1047     671     677    1066+EL   IF
+EL   IF    1067     676                                                        
CHEXA        664       2     671    1047    1048     670     676    1067+EL   IG
+EL   IG    1068     675                                                        
CHEXA        665       2     670    1048    1049     669     675    1068+EL   IH
+EL   IH    1069     674                                                        
CHEXA        666       2     669    1049     696     663     674    1069+EL   II
+EL   II     700     686                                                        
CHEXA        667       2     585     581    1050    1045     586     582+EL   IJ
+EL   IJ    1070    1065                                                        
CHEXA        668       2    1045    1050    1051    1046    1065    1070+EL   IK
+EL   IK    1071    1066                                                        
CHEXA        669       2    1046    1051    1052    1047    1066    1071+EL   IL
+EL   IL    1072    1067                                                        
CHEXA        670       2    1047    1052    1053    1048    1067    1072+EL   IM
+EL   IM    1073    1068                                                        
CHEXA        671       2    1048    1053    1054    1049    1068    1073+EL   IN
+EL   IN    1074    1069                                                        
CHEXA        672       2    1049    1054     695     696    1069    1074+EL   IO
+EL   IO     699     700                                                        
CHEXA        673       2     581     577    1055    1050     582     578+EL   IP
+EL   IP    1075    1070                                                        
CHEXA        674       2    1050    1055    1056    1051    1070    1075+EL   IQ
+EL   IQ    1076    1071                                                        
CHEXA        675       2    1051    1056    1057    1052    1071    1076+EL   IR
+EL   IR    1077    1072                                                        
CHEXA        676       2    1052    1057    1058    1053    1072    1077+EL   IS
+EL   IS    1078    1073                                                        
CHEXA        677       2    1053    1058    1059    1054    1073    1078+EL   IT
+EL   IT    1079    1074                                                        
CHEXA        678       2    1054    1059     694     695    1074    1079+EL   IU
+EL   IU     698     699                                                        
CHEXA        679       2     577     573    1060    1055     578     574+EL   IV
+EL   IV    1080    1075                                                        
CHEXA        680       2    1055    1060    1061    1056    1075    1080+EL   IW
+EL   IW    1081    1076                                                        
CHEXA        681       2    1056    1061    1062    1057    1076    1081+EL   IX
+EL   IX    1082    1077                                                        
CHEXA        682       2    1057    1062    1063    1058    1077    1082+EL   IY
+EL   IY    1083    1078                                                        
CHEXA        683       2    1058    1063    1064    1059    1078    1083+EL   IZ
+EL   IZ    1084    1079                                                        
CHEXA        684       2    1059    1064     693     694    1079    1084+EL   J0
+EL   J0     697     698                                                        
CHEXA        685       2     573     628     638    1060     574     627+EL   J1
+EL   J1     637    1080                                                        
CHEXA        686       2    1060     638     642    1061    1080     637+EL   J2
+EL   J2     641    1081                                                        
CHEXA        687       2    1061     642     646    1062    1081     641+EL   J3
+EL   J3     645    1082                                                        
CHEXA        688       2    1062     646     650    1063    1082     645+EL   J4
+EL   J4     649    1083                                                        
CHEXA        689       2    1063     650     654    1064    1083     649+EL   J5
+EL   J5     653    1084                                                        
CHEXA        690       2    1064     654     633     693    1084     653+EL   J6
+EL   J6     634     697                                                        
CHEXA        691       2     569     586    1065     678     660     587+EL   J7
+EL   J7    1085     683                                                        
CHEXA        692       2     678    1065    1066     677     683    1085+EL   J8
+EL   J8    1086     682                                                        
CHEXA        693       2     677    1066    1067     676     682    1086+EL   J9
+EL   J9    1087     681                                                        
CHEXA        694       2     676    1067    1068     675     681    1087+EL   JA
+EL   JA    1088     680                                                        
CHEXA        695       2     675    1068    1069     674     680    1088+EL   JB
+EL   JB    1089     679                                                        
CHEXA        696       2     674    1069     700     686     679    1089+EL   JC
+EL   JC     704     687                                                        
CHEXA        697       2     586     582    1070    1065     587     583+EL   JD
+EL   JD    1090    1085                                                        
CHEXA        698       2    1065    1070    1071    1066    1085    1090+EL   JE
+EL   JE    1091    1086                                                        
CHEXA        699       2    1066    1071    1072    1067    1086    1091+EL   JF
+EL   JF    1092    1087                                                        
CHEXA        700       2    1067    1072    1073    1068    1087    1092+EL   JG
+EL   JG    1093    1088                                                        
CHEXA        701       2    1068    1073    1074    1069    1088    1093+EL   JH
+EL   JH    1094    1089                                                        
CHEXA        702       2    1069    1074     699     700    1089    1094+EL   JI
+EL   JI     703     704                                                        
CHEXA        703       2     582     578    1075    1070     583     579+EL   JJ
+EL   JJ    1095    1090                                                        
CHEXA        704       2    1070    1075    1076    1071    1090    1095+EL   JK
+EL   JK    1096    1091                                                        
CHEXA        705       2    1071    1076    1077    1072    1091    1096+EL   JL
+EL   JL    1097    1092                                                        
CHEXA        706       2    1072    1077    1078    1073    1092    1097+EL   JM
+EL   JM    1098    1093                                                        
CHEXA        707       2    1073    1078    1079    1074    1093    1098+EL   JN
+EL   JN    1099    1094                                                        
CHEXA        708       2    1074    1079     698     699    1094    1099+EL   JO
+EL   JO     702     703                                                        
CHEXA        709       2     578     574    1080    1075     579     575+EL   JP
+EL   JP    1100    1095                                                        
CHEXA        710       2    1075    1080    1081    1076    1095    1100+EL   JQ
+EL   JQ    1101    1096                                                        
CHEXA        711       2    1076    1081    1082    1077    1096    1101+EL   JR
+EL   JR    1102    1097                                                        
CHEXA        712       2    1077    1082    1083    1078    1097    1102+EL   JS
+EL   JS    1103    1098                                                        
CHEXA        713       2    1078    1083    1084    1079    1098    1103+EL   JT
+EL   JT    1104    1099                                                        
CHEXA        714       2    1079    1084     697     698    1099    1104+EL   JU
+EL   JU     701     702                                                        
CHEXA        715       2     574     627     637    1080     575     626+EL   JV
+EL   JV     636    1100                                                        
CHEXA        716       2    1080     637     641    1081    1100     636+EL   JW
+EL   JW     640    1101                                                        
CHEXA        717       2    1081     641     645    1082    1101     640+EL   JX
+EL   JX     644    1102                                                        
CHEXA        718       2    1082     645     649    1083    1102     644+EL   JY
+EL   JY     648    1103                                                        
CHEXA        719       2    1083     649     653    1084    1103     648+EL   JZ
+EL   JZ     652    1104                                                        
CHEXA        720       2    1084     653     634     697    1104     652+EL   K0
+EL   K0     688     701                                                        
CHEXA        721       2     660     587    1085     683     661     568+EL   K1
+EL   K1     609     595                                                        
CHEXA        722       2     683    1085    1086     682     595     609+EL   K2
+EL   K2     340     662                                                        
CHEXA        723       2     682    1086    1087     681     662     340+EL   K3
+EL   K3     339     594                                                        
CHEXA        724       2     681    1087    1088     680     594     339+EL   K4
+EL   K4     610     593                                                        
CHEXA        725       2     680    1088    1089     679     593     610+EL   K5
+EL   K5     338     333                                                        
CHEXA        726       2     679    1089     704     687     333     338+EL   K6
+EL   K6     592     275                                                        
CHEXA        727       2     587     583    1090    1085     568     596+EL   K7
+EL   K7     337     609                                                        
CHEXA        728       2    1085    1090    1091    1086     609     337+EL   K8
+EL   K8     606     340                                                        
CHEXA        729       2    1086    1091    1092    1087     340     606+EL   K9
+EL   K9     607     339                                                        
CHEXA        730       2    1087    1092    1093    1088     339     607+EL   KA
+EL   KA     336     610                                                        
CHEXA        731       2    1088    1093    1094    1089     610     336+EL   KB
+EL   KB     608     338                                                        
CHEXA        732       2    1089    1094     703     704     338     608+EL   KC
+EL   KC     591     592                                                        
CHEXA        733       2     583     579    1095    1090     596     597+EL   KD
+EL   KD     335     337                                                        
CHEXA        734       2    1090    1095    1096    1091     337     335+EL   KE
+EL   KE     603     606                                                        
CHEXA        735       2    1091    1096    1097    1092     606     603+EL   KF
+EL   KF     604     607                                                        
CHEXA        736       2    1092    1097    1098    1093     607     604+EL   KG
+EL   KG     605     336                                                        
CHEXA        737       2    1093    1098    1099    1094     336     605+EL   KH
+EL   KH     334     608                                                        
CHEXA        738       2    1094    1099     702     703     608     334+EL   KI
+EL   KI     276     591                                                        
CHEXA        739       2     579     575    1100    1095     597     332+EL   KJ
+EL   KJ     598     335                                                        
CHEXA        740       2    1095    1100    1101    1096     335     598+EL   KK
+EL   KK     599     603                                                        
CHEXA        741       2    1096    1101    1102    1097     603     599+EL   KL
+EL   KL     600     604                                                        
CHEXA        742       2    1097    1102    1103    1098     604     600+EL   KM
+EL   KM     601     605                                                        
CHEXA        743       2    1098    1103    1104    1099     605     601+EL   KN
+EL   KN     602     334                                                        
CHEXA        744       2    1099    1104     701     702     334     602+EL   KO
+EL   KO     590     276                                                        
CHEXA        745       2     575     626     636    1100     332     567+EL   KP
+EL   KP     588     598                                                        
CHEXA        746       2    1100     636     640    1101     598     588+EL   KQ
+EL   KQ     635     599                                                        
CHEXA        747       2    1101     640     644    1102     599     635+EL   KR
+EL   KR     368     600                                                        
CHEXA        748       2    1102     644     648    1103     600     368+EL   KS
+EL   KS     331     601                                                        
CHEXA        749       2    1103     648     652    1104     601     331+EL   KT
+EL   KT     330     602                                                        
CHEXA        750       2    1104     652     688     701     602     330+EL   KU
+EL   KU     589     590                                                        
ENDDATA
