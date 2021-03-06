ID C:\Users\Dima\Lame3d2\Mpl,FEMAP
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
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\Ono_box\Solid\ono_box0.MOD
$   Date       : Sat Mar 29 11:32:47 2008
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
$ FEMAP Load Set 1 : Srcplane Vn -1.0 0.0
PLOAD4         1     255    111.                             389     225
$ FEMAP Load Set 2 : Endwall Absorb 1.0 0.0
PLOAD4         2     233    222.                              52     402
PLOAD4         2     244    222.                               6     378
PLOAD4         2     265    222.                             378     201
PLOAD4         2     276    222.                             378     332
PLOAD4         2     277    222.                             378     447
PLOAD4         2     278    222.                             402     448
$ FEMAP Property 1 : ONO BOX
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : INCLUSION
PSHELL         2       1      0.       1               1              0.
$ FEMAP Property 3 : Untitled
PSOLID         3       1       0        
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           4       0      0.    0.25     2.5       0        
GRID           6       0   0.175      0.     2.5       0        
GRID          11       0   0.175      0.      0.       0        
GRID          16       0      0.    0.25 2.04545       0        
GRID          17       0      0.    0.25 1.81818       0        
GRID          23       0      0.    0.25 0.45455       0        
GRID          24       0      0.    0.25 0.22727       0        
GRID          26       0      0.      0.      0.       0        
GRID          27       0      0.      0. 0.22727       0        
GRID          28       0      0.      0. 0.45455       0        
GRID          29       0      0.      0. 0.68182       0        
GRID          31       0      0.      0. 1.13636       0        
GRID          32       0      0.      0. 1.36364       0        
GRID          33       0      0.      0. 1.59091       0        
GRID          36       0      0.      0. 2.27273       0        
GRID          43       0    0.35      0. 0.90909       0        
GRID          46       0    0.35      0. 1.59091       0        
GRID          52       0      0.      0.     2.5       0        
GRID          54       0      0.      0. 2.04545       0        
GRID          55       0      0.      0. 1.81818       0        
GRID          59       0      0.      0. 0.90909       0        
GRID          63       0   0.175      0. 0.22727       0        
GRID          64       0   0.175      0. 0.45455       0        
GRID          65       0   0.175      0. 0.68182       0        
GRID          66       0   0.175      0. 0.90909       0        
GRID          67       0   0.175      0. 1.13636       0        
GRID          68       0   0.175      0. 1.36364       0        
GRID          69       0   0.175      0. 1.59091       0        
GRID          70       0   0.175      0. 1.81818       0        
GRID          71       0   0.175      0. 2.04545       0        
GRID          72       0   0.175      0. 2.27273       0        
GRID          74       0    0.35      0. 2.27273       0        
GRID          76       0    0.35      0. 1.81818       0        
GRID          78       0    0.35      0. 1.36364       0        
GRID          81       0    0.35      0. 0.68182       0        
GRID          87       0    0.35    0.25 0.45455       0        
GRID          88       0    0.35    0.25 0.68182       0        
GRID          91       0    0.35    0.25 1.36364       0        
GRID         100       0    0.35    0.25 2.27273       0        
GRID         101       0    0.35    0.25 2.04545       0        
GRID         111       0   0.175    0.25      0.       0        
GRID         112       0      0.    0.25      0.       0        
GRID         116       0      0.    0.25 0.90909       0        
GRID         117       0      0.    0.25 1.13636       0        
GRID         119       0      0.    0.25 1.59091       0        
GRID         124       0   0.175    0.25 2.04545       0        
GRID         126       0   0.175    0.25 1.59091       0        
GRID         128       0   0.175    0.25 1.13636       0        
GRID         130       0   0.175    0.25 0.68182       0        
GRID         131       0   0.175    0.25 0.45455       0        
GRID         132       0   0.175    0.25 0.22727       0        
GRID         141       0    0.35     0.5 0.90909       0        
GRID         151       0    0.35    0.25 1.13636       0        
GRID         154       0    0.35    0.25 1.81818       0        
GRID         159       0     0.6    0.25 0.22727       0        
GRID         167       0     0.6    0.25 2.04545       0        
GRID         177       0    0.35    0.25 0.90909       0        
GRID         185       0    0.35     0.5      0.       0        
GRID         189       0     0.6    0.25      0.       0        
GRID         192       0     0.6     0.5 0.45455       0        
GRID         194       0     0.6     0.5 0.90909       0        
GRID         195       0     0.6     0.5 1.13636       0        
GRID         197       0     0.6     0.5 1.59091       0        
GRID         198       0     0.6     0.5 1.81818       0        
GRID         200       0     0.6     0.5 2.27273       0        
GRID         201       0     0.6     0.5     2.5       0        
GRID         203       0     0.6    0.25 2.27273       0        
GRID         205       0     0.6    0.25 1.81818       0        
GRID         206       0     0.6    0.25 1.59091       0        
GRID         216       0     0.6     0.5 2.04545       0        
GRID         219       0     0.6     0.5 1.36364       0        
GRID         222       0     0.6     0.5 0.68182       0        
GRID         224       0     0.6     0.5 0.22727       0        
GRID         225       0     0.6     0.5      0.       0        
GRID         232       0    0.35     0.5 1.36364       0        
GRID         235       0    0.35     0.5 2.04545       0        
GRID         237       0     0.6    0.25     2.5       0        
GRID         242       0     0.6    0.25 1.36364       0        
GRID         243       0     0.6    0.25 1.13636       0        
GRID         244       0     0.6    0.25 0.90909       0        
GRID         264       0     0.6    0.25 0.45455       0        
GRID         265       0     0.6    0.25 0.68182       0        
GRID         276       0     0.6      0. 2.04545       0        
GRID         277       0     0.6      0. 1.81818       0        
GRID         279       0     0.6      0. 1.36364       0        
GRID         281       0     0.6      0. 0.90909       0        
GRID         283       0     0.6      0. 0.45455       0        
GRID         284       0     0.6      0. 0.22727       0        
GRID         301       0    0.35    0.25 1.59091       0        
GRID         309       0    0.35      0.     2.5       0        
GRID         311       0    0.35      0. 2.04545       0        
GRID         315       0    0.35      0. 1.13636       0        
GRID         318       0    0.35      0. 0.45455       0        
GRID         319       0    0.35      0. 0.22727       0        
GRID         320       0    0.35      0.      0.       0        
GRID         324       0     0.6      0. 0.68182       0        
GRID         326       0     0.6      0. 1.13636       0        
GRID         328       0     0.6      0. 1.59091       0        
GRID         331       0     0.6      0. 2.27273       0        
GRID         332       0     0.6      0.     2.5       0        
GRID         333       0     0.6      0.      0.       0        
GRID         344       0      0.     0.5 0.22727       0        
GRID         345       0      0.     0.5 0.45455       0        
GRID         348       0      0.     0.5 1.13636       0        
GRID         349       0      0.     0.5 1.36364       0        
GRID         350       0      0.     0.5 1.59091       0        
GRID         352       0      0.     0.5 2.04545       0        
GRID         357       0    0.35     0.5 2.27273       0        
GRID         359       0    0.35     0.5 1.81818       0        
GRID         360       0    0.35     0.5 1.59091       0        
GRID         362       0    0.35     0.5 1.13636       0        
GRID         364       0    0.35     0.5 0.68182       0        
GRID         365       0    0.35     0.5 0.45455       0        
GRID         366       0    0.35     0.5 0.22727       0        
GRID         367       0   0.175     0.5 0.22727       0        
GRID         368       0   0.175     0.5 0.45455       0        
GRID         369       0   0.175     0.5 0.68182       0        
GRID         370       0   0.175     0.5 0.90909       0        
GRID         371       0   0.175     0.5 1.13636       0        
GRID         372       0   0.175     0.5 1.36364       0        
GRID         373       0   0.175     0.5 1.59091       0        
GRID         374       0   0.175     0.5 1.81818       0        
GRID         375       0   0.175     0.5 2.04545       0        
GRID         376       0   0.175     0.5 2.27273       0        
GRID         378       0    0.35    0.25     2.5       0        
GRID         389       0    0.35    0.25      0.       0        
GRID         402       0   0.175    0.25     2.5       0        
GRID         404       0      0.    0.25 2.27273       0        
GRID         411       0      0.    0.25 0.68182       0        
GRID         417       0    0.35    0.25 0.22727       0        
GRID         427       0   0.175    0.25 2.27273       0        
GRID         429       0   0.175    0.25 1.81818       0        
GRID         431       0   0.175    0.25 1.36364       0        
GRID         433       0   0.175    0.25 0.90909       0        
GRID         442       0   0.175     0.5      0.       0        
GRID         446       0    0.35     0.5     2.5       0        
GRID         447       0   0.175     0.5     2.5       0        
GRID         448       0      0.     0.5     2.5       0        
GRID         455       0      0.    0.25 1.36364       0        
GRID         462       0      0.     0.5 2.27273       0        
GRID         464       0      0.     0.5 1.81818       0        
GRID         468       0      0.     0.5 0.90909       0        
GRID         469       0      0.     0.5 0.68182       0        
GRID         472       0      0.     0.5      0.       0        
CHEXA        233       3      52      36     404       4       6      72+EL   6H
+EL   6H     427     402                                                        
CHEXA        234       3      36      54      16     404      72      71+EL   6I
+EL   6I     124     427                                                        
CHEXA        235       3      54      55      17      16      71      70+EL   6J
+EL   6J     429     124                                                        
CHEXA        236       3      55      33     119      17      70      69+EL   6K
+EL   6K     126     429                                                        
CHEXA        237       3      33      32     455     119      69      68+EL   6L
+EL   6L     431     126                                                        
CHEXA        238       3      32      31     117     455      68      67+EL   6M
+EL   6M     128     431                                                        
CHEXA        239       3      31      59     116     117      67      66+EL   6N
+EL   6N     433     128                                                        
CHEXA        240       3      59      29     411     116      66      65+EL   6O
+EL   6O     130     433                                                        
CHEXA        241       3      29      28      23     411      65      64+EL   6P
+EL   6P     131     130                                                        
CHEXA        242       3      28      27      24      23      64      63+EL   6Q
+EL   6Q     132     131                                                        
CHEXA        243       3      27      26     112      24      63      11+EL   6R
+EL   6R     111     132                                                        
CHEXA        244       3       6      72     427     402     309      74+EL   6S
+EL   6S     100     378                                                        
CHEXA        245       3      72      71     124     427      74     311+EL   6T
+EL   6T     101     100                                                        
CHEXA        246       3      71      70     429     124     311      76+EL   6U
+EL   6U     154     101                                                        
CHEXA        247       3      70      69     126     429      76      46+EL   6V
+EL   6V     301     154                                                        
CHEXA        248       3      69      68     431     126      46      78+EL   6W
+EL   6W      91     301                                                        
CHEXA        249       3      68      67     128     431      78     315+EL   6X
+EL   6X     151      91                                                        
CHEXA        250       3      67      66     433     128     315      43+EL   6Y
+EL   6Y     177     151                                                        
CHEXA        251       3      66      65     130     433      43      81+EL   6Z
+EL   6Z      88     177                                                        
CHEXA        252       3      65      64     131     130      81     318+EL   70
+EL   70      87      88                                                        
CHEXA        253       3      64      63     132     131     318     319+EL   71
+EL   71     417      87                                                        
CHEXA        254       3      63      11     111     132     319     320+EL   72
+EL   72     389     417                                                        
CHEXA        255       3     389     417     159     189     185     366+EL   73
+EL   73     224     225                                                        
CHEXA        256       3     417      87     264     159     366     365+EL   74
+EL   74     192     224                                                        
CHEXA        257       3      87      88     265     264     365     364+EL   75
+EL   75     222     192                                                        
CHEXA        258       3      88     177     244     265     364     141+EL   76
+EL   76     194     222                                                        
CHEXA        259       3     177     151     243     244     141     362+EL   77
+EL   77     195     194                                                        
CHEXA        260       3     151      91     242     243     362     232+EL   78
+EL   78     219     195                                                        
CHEXA        261       3      91     301     206     242     232     360+EL   79
+EL   79     197     219                                                        
CHEXA        262       3     301     154     205     206     360     359+EL   7A
+EL   7A     198     197                                                        
CHEXA        263       3     154     101     167     205     359     235+EL   7B
+EL   7B     216     198                                                        
CHEXA        264       3     101     100     203     167     235     357+EL   7C
+EL   7C     200     216                                                        
CHEXA        265       3     100     378     237     203     357     446+EL   7D
+EL   7D     201     200                                                        
CHEXA        266       3     320     389     417     319     333     189+EL   7E
+EL   7E     159     284                                                        
CHEXA        267       3     319     417      87     318     284     159+EL   7F
+EL   7F     264     283                                                        
CHEXA        268       3     318      87      88      81     283     264+EL   7G
+EL   7G     265     324                                                        
CHEXA        269       3      81      88     177      43     324     265+EL   7H
+EL   7H     244     281                                                        
CHEXA        270       3      43     177     151     315     281     244+EL   7I
+EL   7I     243     326                                                        
CHEXA        271       3     315     151      91      78     326     243+EL   7J
+EL   7J     242     279                                                        
CHEXA        272       3      78      91     301      46     279     242+EL   7K
+EL   7K     206     328                                                        
CHEXA        273       3      46     301     154      76     328     206+EL   7L
+EL   7L     205     277                                                        
CHEXA        274       3      76     154     101     311     277     205+EL   7M
+EL   7M     167     276                                                        
CHEXA        275       3     311     101     100      74     276     167+EL   7N
+EL   7N     203     331                                                        
CHEXA        276       3      74     100     378     309     331     203+EL   7O
+EL   7O     237     332                                                        
CHEXA        277       3     378     100     427     402     446     357+EL   7P
+EL   7P     376     447                                                        
CHEXA        278       3     402     427     404       4     447     376+EL   7Q
+EL   7Q     462     448                                                        
CHEXA        279       3     100     101     124     427     357     235+EL   7R
+EL   7R     375     376                                                        
CHEXA        280       3     427     124      16     404     376     375+EL   7S
+EL   7S     352     462                                                        
CHEXA        281       3     101     154     429     124     235     359+EL   7T
+EL   7T     374     375                                                        
CHEXA        282       3     124     429      17      16     375     374+EL   7U
+EL   7U     464     352                                                        
CHEXA        283       3     154     301     126     429     359     360+EL   7V
+EL   7V     373     374                                                        
CHEXA        284       3     429     126     119      17     374     373+EL   7W
+EL   7W     350     464                                                        
CHEXA        285       3     301      91     431     126     360     232+EL   7X
+EL   7X     372     373                                                        
CHEXA        286       3     126     431     455     119     373     372+EL   7Y
+EL   7Y     349     350                                                        
CHEXA        287       3      91     151     128     431     232     362+EL   7Z
+EL   7Z     371     372                                                        
CHEXA        288       3     431     128     117     455     372     371+EL   80
+EL   80     348     349                                                        
CHEXA        289       3     151     177     433     128     362     141+EL   81
+EL   81     370     371                                                        
CHEXA        290       3     128     433     116     117     371     370+EL   82
+EL   82     468     348                                                        
CHEXA        291       3     177      88     130     433     141     364+EL   83
+EL   83     369     370                                                        
CHEXA        292       3     433     130     411     116     370     369+EL   84
+EL   84     469     468                                                        
CHEXA        293       3      88      87     131     130     364     365+EL   85
+EL   85     368     369                                                        
CHEXA        294       3     130     131      23     411     369     368+EL   86
+EL   86     345     469                                                        
CHEXA        295       3      87     417     132     131     365     366+EL   87
+EL   87     367     368                                                        
CHEXA        296       3     131     132      24      23     368     367+EL   88
+EL   88     344     345                                                        
CHEXA        297       3     417     389     111     132     366     185+EL   89
+EL   89     442     367                                                        
CHEXA        298       3     132     111     112      24     367     442+EL   8A
+EL   8A     472     344                                                        
ENDDATA
