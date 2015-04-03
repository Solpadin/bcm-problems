ID FEMAP,FEMAP
SOL SESTATICS
TIME 10000
CEND
  ECHO = NONE
  DISPLACEMENT(PLOT) = ALL
  SPCFORCE(PLOT) = ALL
  OLOAD(PLOT) = ALL
  FORCE(PLOT,CORNER) = ALL
  STRESS(PLOT,CORNER) = ALL
BEGIN BULK
$ ***************************************************************************
$   Written by : FEMAP
$   Version    : 7.00
$   Translator : MSC/NASTRAN
$   From Model : 
$   Date       : Fri Apr 09 10:36:35 2010
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
$ FEMAP Property 1 : Untitled
PSOLID         1       1       0        
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           1       0     0.4      0.      1.       0        
GRID           2       0     0.6      0.      1.       0        
GRID           3       0     0.8      0.      1.       0        
GRID           4       0      1.     0.4      1.       0        
GRID           5       0      1.      1.      1.       0        
GRID           6       0     0.6      1.      1.       0        
GRID           7       0      0.     0.6      1.       0        
GRID           8       0      0.     0.2      1.       0        
GRID           9       0     0.2     0.2      1.       0        
GRID          10       0     0.4     0.2      1.       0        
GRID          11       0     0.6     0.2      1.       0        
GRID          12       0     0.8     0.2      1.       0        
GRID          13       0     0.2     0.4      1.       0        
GRID          14       0     0.4     0.4      1.       0        
GRID          15       0     0.6     0.4      1.       0        
GRID          16       0     0.8     0.4      1.       0        
GRID          17       0     0.2     0.6      1.       0        
GRID          18       0     0.4     0.6      1.       0        
GRID          19       0     0.6     0.6      1.       0        
GRID          20       0     0.8     0.6      1.       0        
GRID          21       0     0.2     0.8      1.       0        
GRID          22       0     0.4     0.8      1.       0        
GRID          23       0     0.6     0.8      1.       0        
GRID          24       0     0.8     0.8      1.       0        
GRID          25       0      1.      0.     0.4       0        
GRID          26       0      1.      0.     0.6       0        
GRID          27       0      1.      0.      1.       0        
GRID          28       0     0.2      0.      1.       0        
GRID          29       0      0.      0.     0.8       0        
GRID          30       0      0.      0.     0.2       0        
GRID          31       0     0.2      0.      0.       0        
GRID          32       0     0.4      0.      0.       0        
GRID          33       0     0.8      0.      0.       0        
GRID          34       0     0.8      0.     0.2       0        
GRID          35       0     0.8      0.     0.4       0        
GRID          36       0     0.8      0.     0.6       0        
GRID          37       0     0.8      0.     0.8       0        
GRID          38       0     0.6      0.     0.2       0        
GRID          39       0     0.6      0.     0.4       0        
GRID          40       0     0.6      0.     0.6       0        
GRID          41       0     0.6      0.     0.8       0        
GRID          42       0     0.4      0.     0.2       0        
GRID          43       0     0.4      0.     0.4       0        
GRID          44       0     0.4      0.     0.6       0        
GRID          45       0     0.4      0.     0.8       0        
GRID          46       0     0.2      0.     0.2       0        
GRID          47       0     0.2      0.     0.4       0        
GRID          48       0     0.2      0.     0.6       0        
GRID          49       0     0.2      0.     0.8       0        
GRID          50       0      0.      0.      0.       0        
GRID          51       0      0.      0.     0.4       0        
GRID          52       0      0.      0.     0.6       0        
GRID          53       0      0.      0.      1.       0        
GRID          54       0      0.     0.4      1.       0        
GRID          55       0      0.     0.8      1.       0        
GRID          56       0      0.      1.      1.       0        
GRID          57       0      0.      1.     0.6       0        
GRID          58       0      0.     0.6      0.       0        
GRID          59       0      0.     0.2     0.2       0        
GRID          60       0      0.     0.2     0.4       0        
GRID          61       0      0.     0.2     0.6       0        
GRID          62       0      0.     0.2     0.8       0        
GRID          63       0      0.     0.4     0.2       0        
GRID          64       0      0.     0.4     0.4       0        
GRID          65       0      0.     0.4     0.6       0        
GRID          66       0      0.     0.4     0.8       0        
GRID          67       0      0.     0.6     0.2       0        
GRID          68       0      0.     0.6     0.4       0        
GRID          69       0      0.     0.6     0.6       0        
GRID          70       0      0.     0.6     0.8       0        
GRID          71       0      0.     0.8     0.2       0        
GRID          72       0      0.     0.8     0.4       0        
GRID          73       0      0.     0.8     0.6       0        
GRID          74       0      0.     0.8     0.8       0        
GRID          75       0      0.      1.      0.       0        
GRID          76       0      0.      1.     0.2       0        
GRID          77       0      0.      1.     0.4       0        
GRID          78       0      0.      1.     0.8       0        
GRID          79       0     0.2      1.      1.       0        
GRID          80       0     0.4      1.      1.       0        
GRID          81       0     0.8      1.      1.       0        
GRID          82       0      1.      1.     0.8       0        
GRID          83       0      1.      1.     0.6       0        
GRID          84       0      1.      1.     0.4       0        
GRID          85       0      1.      1.      0.       0        
GRID          86       0     0.8      1.      0.       0        
GRID          87       0     0.2      1.      0.       0        
GRID          88       0     0.2      1.     0.2       0        
GRID          89       0     0.2      1.     0.4       0        
GRID          90       0     0.2      1.     0.6       0        
GRID          91       0     0.2      1.     0.8       0        
GRID          92       0     0.4      1.     0.2       0        
GRID          93       0     0.4      1.     0.4       0        
GRID          94       0     0.4      1.     0.6       0        
GRID          95       0     0.4      1.     0.8       0        
GRID          96       0     0.6      1.     0.2       0        
GRID          97       0     0.6      1.     0.4       0        
GRID          98       0     0.6      1.     0.6       0        
GRID          99       0     0.6      1.     0.8       0        
GRID         100       0     0.8      1.     0.2       0        
GRID         101       0     0.8      1.     0.4       0        
GRID         102       0     0.8      1.     0.6       0        
GRID         103       0     0.8      1.     0.8       0        
GRID         104       0      1.     0.6      0.       0        
GRID         105       0      1.     0.2      0.       0        
GRID         106       0      1.      0.      0.       0        
GRID         107       0     0.6      0.      0.       0        
GRID         108       0      0.     0.2      0.       0        
GRID         109       0      0.     0.4      0.       0        
GRID         110       0      0.     0.8      0.       0        
GRID         111       0     0.4      1.      0.       0        
GRID         112       0     0.6      1.      0.       0        
GRID         113       0     0.8     0.8      0.       0        
GRID         114       0     0.8     0.6      0.       0        
GRID         115       0     0.8     0.4      0.       0        
GRID         116       0     0.8     0.2      0.       0        
GRID         117       0     0.6     0.8      0.       0        
GRID         118       0     0.6     0.6      0.       0        
GRID         119       0     0.6     0.4      0.       0        
GRID         120       0     0.6     0.2      0.       0        
GRID         121       0     0.4     0.8      0.       0        
GRID         122       0     0.4     0.6      0.       0        
GRID         123       0     0.4     0.4      0.       0        
GRID         124       0     0.4     0.2      0.       0        
GRID         125       0     0.2     0.8      0.       0        
GRID         126       0     0.2     0.6      0.       0        
GRID         127       0     0.2     0.4      0.       0        
GRID         128       0     0.2     0.2      0.       0        
GRID         129       0      1.     0.4      0.       0        
GRID         130       0      1.     0.8      0.       0        
GRID         131       0      1.      1.     0.2       0        
GRID         132       0      1.     0.8      1.       0        
GRID         133       0      1.     0.6      1.       0        
GRID         134       0      1.     0.2      1.       0        
GRID         135       0      1.      0.     0.8       0        
GRID         136       0      1.      0.     0.2       0        
GRID         137       0      1.     0.2     0.2       0        
GRID         138       0      1.     0.4     0.2       0        
GRID         139       0      1.     0.6     0.2       0        
GRID         140       0      1.     0.8     0.2       0        
GRID         141       0      1.     0.2     0.4       0        
GRID         142       0      1.     0.4     0.4       0        
GRID         143       0      1.     0.6     0.4       0        
GRID         144       0      1.     0.8     0.4       0        
GRID         145       0      1.     0.2     0.6       0        
GRID         146       0      1.     0.4     0.6       0        
GRID         147       0      1.     0.6     0.6       0        
GRID         148       0      1.     0.8     0.6       0        
GRID         149       0      1.     0.2     0.8       0        
GRID         150       0      1.     0.4     0.8       0        
GRID         151       0      1.     0.6     0.8       0        
GRID         152       0      1.     0.8     0.8       0        
GRID         153       0     0.2     0.2     0.8       0        
GRID         154       0     0.4     0.2     0.8       0        
GRID         155       0     0.6     0.2     0.8       0        
GRID         156       0     0.8     0.2     0.8       0        
GRID         157       0     0.2     0.4     0.8       0        
GRID         158       0     0.4     0.4     0.8       0        
GRID         159       0     0.6     0.4     0.8       0        
GRID         160       0     0.8     0.4     0.8       0        
GRID         161       0     0.2     0.6     0.8       0        
GRID         162       0     0.4     0.6     0.8       0        
GRID         163       0     0.6     0.6     0.8       0        
GRID         164       0     0.8     0.6     0.8       0        
GRID         165       0     0.2     0.8     0.8       0        
GRID         166       0     0.4     0.8     0.8       0        
GRID         167       0     0.6     0.8     0.8       0        
GRID         168       0     0.8     0.8     0.8       0        
GRID         169       0     0.2     0.2     0.6       0        
GRID         170       0     0.4     0.2     0.6       0        
GRID         171       0     0.6     0.2     0.6       0        
GRID         172       0     0.8     0.2     0.6       0        
GRID         173       0     0.2     0.4     0.6       0        
GRID         174       0     0.4     0.4     0.6       0        
GRID         175       0     0.6     0.4     0.6       0        
GRID         176       0     0.8     0.4     0.6       0        
GRID         177       0     0.2     0.6     0.6       0        
GRID         178       0     0.4     0.6     0.6       0        
GRID         179       0     0.6     0.6     0.6       0        
GRID         180       0     0.8     0.6     0.6       0        
GRID         181       0     0.2     0.8     0.6       0        
GRID         182       0     0.4     0.8     0.6       0        
GRID         183       0     0.6     0.8     0.6       0        
GRID         184       0     0.8     0.8     0.6       0        
GRID         185       0     0.2     0.2     0.4       0        
GRID         186       0     0.4     0.2     0.4       0        
GRID         187       0     0.6     0.2     0.4       0        
GRID         188       0     0.8     0.2     0.4       0        
GRID         189       0     0.2     0.4     0.4       0        
GRID         190       0     0.4     0.4     0.4       0        
GRID         191       0     0.6     0.4     0.4       0        
GRID         192       0     0.8     0.4     0.4       0        
GRID         193       0     0.2     0.6     0.4       0        
GRID         194       0     0.4     0.6     0.4       0        
GRID         195       0     0.6     0.6     0.4       0        
GRID         196       0     0.8     0.6     0.4       0        
GRID         197       0     0.2     0.8     0.4       0        
GRID         198       0     0.4     0.8     0.4       0        
GRID         199       0     0.6     0.8     0.4       0        
GRID         200       0     0.8     0.8     0.4       0        
GRID         201       0     0.2     0.2     0.2       0        
GRID         202       0     0.4     0.2     0.2       0        
GRID         203       0     0.6     0.2     0.2       0        
GRID         204       0     0.8     0.2     0.2       0        
GRID         205       0     0.2     0.4     0.2       0        
GRID         206       0     0.4     0.4     0.2       0        
GRID         207       0     0.6     0.4     0.2       0        
GRID         208       0     0.8     0.4     0.2       0        
GRID         209       0     0.2     0.6     0.2       0        
GRID         210       0     0.4     0.6     0.2       0        
GRID         211       0     0.6     0.6     0.2       0        
GRID         212       0     0.8     0.6     0.2       0        
GRID         213       0     0.2     0.8     0.2       0        
GRID         214       0     0.4     0.8     0.2       0        
GRID         215       0     0.6     0.8     0.2       0        
GRID         216       0     0.8     0.8     0.2       0        
CHEXA          1       1      53       8       9      28      29      62+EL    1
+EL    1     153      49                                                        
CHEXA          2       1      28       9      10       1      49     153+EL    2
+EL    2     154      45                                                        
CHEXA          3       1       1      10      11       2      45     154+EL    3
+EL    3     155      41                                                        
CHEXA          4       1       2      11      12       3      41     155+EL    4
+EL    4     156      37                                                        
CHEXA          5       1       3      12     134      27      37     156+EL    5
+EL    5     149     135                                                        
CHEXA          6       1       8      54      13       9      62      66+EL    6
+EL    6     157     153                                                        
CHEXA          7       1       9      13      14      10     153     157+EL    7
+EL    7     158     154                                                        
CHEXA          8       1      10      14      15      11     154     158+EL    8
+EL    8     159     155                                                        
CHEXA          9       1      11      15      16      12     155     159+EL    9
+EL    9     160     156                                                        
CHEXA         10       1      12      16       4     134     156     160+EL    A
+EL    A     150     149                                                        
CHEXA         11       1      54       7      17      13      66      70+EL    B
+EL    B     161     157                                                        
CHEXA         12       1      13      17      18      14     157     161+EL    C
+EL    C     162     158                                                        
CHEXA         13       1      14      18      19      15     158     162+EL    D
+EL    D     163     159                                                        
CHEXA         14       1      15      19      20      16     159     163+EL    E
+EL    E     164     160                                                        
CHEXA         15       1      16      20     133       4     160     164+EL    F
+EL    F     151     150                                                        
CHEXA         16       1       7      55      21      17      70      74+EL    G
+EL    G     165     161                                                        
CHEXA         17       1      17      21      22      18     161     165+EL    H
+EL    H     166     162                                                        
CHEXA         18       1      18      22      23      19     162     166+EL    I
+EL    I     167     163                                                        
CHEXA         19       1      19      23      24      20     163     167+EL    J
+EL    J     168     164                                                        
CHEXA         20       1      20      24     132     133     164     168+EL    K
+EL    K     152     151                                                        
CHEXA         21       1      55      56      79      21      74      78+EL    L
+EL    L      91     165                                                        
CHEXA         22       1      21      79      80      22     165      91+EL    M
+EL    M      95     166                                                        
CHEXA         23       1      22      80       6      23     166      95+EL    N
+EL    N      99     167                                                        
CHEXA         24       1      23       6      81      24     167      99+EL    O
+EL    O     103     168                                                        
CHEXA         25       1      24      81       5     132     168     103+EL    P
+EL    P      82     152                                                        
CHEXA         26       1      29      62     153      49      52      61+EL    Q
+EL    Q     169      48                                                        
CHEXA         27       1      49     153     154      45      48     169+EL    R
+EL    R     170      44                                                        
CHEXA         28       1      45     154     155      41      44     170+EL    S
+EL    S     171      40                                                        
CHEXA         29       1      41     155     156      37      40     171+EL    T
+EL    T     172      36                                                        
CHEXA         30       1      37     156     149     135      36     172+EL    U
+EL    U     145      26                                                        
CHEXA         31       1      62      66     157     153      61      65+EL    V
+EL    V     173     169                                                        
CHEXA         32       1     153     157     158     154     169     173+EL    W
+EL    W     174     170                                                        
CHEXA         33       1     154     158     159     155     170     174+EL    X
+EL    X     175     171                                                        
CHEXA         34       1     155     159     160     156     171     175+EL    Y
+EL    Y     176     172                                                        
CHEXA         35       1     156     160     150     149     172     176+EL    Z
+EL    Z     146     145                                                        
CHEXA         36       1      66      70     161     157      65      69+EL   10
+EL   10     177     173                                                        
CHEXA         37       1     157     161     162     158     173     177+EL   11
+EL   11     178     174                                                        
CHEXA         38       1     158     162     163     159     174     178+EL   12
+EL   12     179     175                                                        
CHEXA         39       1     159     163     164     160     175     179+EL   13
+EL   13     180     176                                                        
CHEXA         40       1     160     164     151     150     176     180+EL   14
+EL   14     147     146                                                        
CHEXA         41       1      70      74     165     161      69      73+EL   15
+EL   15     181     177                                                        
CHEXA         42       1     161     165     166     162     177     181+EL   16
+EL   16     182     178                                                        
CHEXA         43       1     162     166     167     163     178     182+EL   17
+EL   17     183     179                                                        
CHEXA         44       1     163     167     168     164     179     183+EL   18
+EL   18     184     180                                                        
CHEXA         45       1     164     168     152     151     180     184+EL   19
+EL   19     148     147                                                        
CHEXA         46       1      74      78      91     165      73      57+EL   1A
+EL   1A      90     181                                                        
CHEXA         47       1     165      91      95     166     181      90+EL   1B
+EL   1B      94     182                                                        
CHEXA         48       1     166      95      99     167     182      94+EL   1C
+EL   1C      98     183                                                        
CHEXA         49       1     167      99     103     168     183      98+EL   1D
+EL   1D     102     184                                                        
CHEXA         50       1     168     103      82     152     184     102+EL   1E
+EL   1E      83     148                                                        
CHEXA         51       1      52      61     169      48      51      60+EL   1F
+EL   1F     185      47                                                        
CHEXA         52       1      48     169     170      44      47     185+EL   1G
+EL   1G     186      43                                                        
CHEXA         53       1      44     170     171      40      43     186+EL   1H
+EL   1H     187      39                                                        
CHEXA         54       1      40     171     172      36      39     187+EL   1I
+EL   1I     188      35                                                        
CHEXA         55       1      36     172     145      26      35     188+EL   1J
+EL   1J     141      25                                                        
CHEXA         56       1      61      65     173     169      60      64+EL   1K
+EL   1K     189     185                                                        
CHEXA         57       1     169     173     174     170     185     189+EL   1L
+EL   1L     190     186                                                        
CHEXA         58       1     170     174     175     171     186     190+EL   1M
+EL   1M     191     187                                                        
CHEXA         59       1     171     175     176     172     187     191+EL   1N
+EL   1N     192     188                                                        
CHEXA         60       1     172     176     146     145     188     192+EL   1O
+EL   1O     142     141                                                        
CHEXA         61       1      65      69     177     173      64      68+EL   1P
+EL   1P     193     189                                                        
CHEXA         62       1     173     177     178     174     189     193+EL   1Q
+EL   1Q     194     190                                                        
CHEXA         63       1     174     178     179     175     190     194+EL   1R
+EL   1R     195     191                                                        
CHEXA         64       1     175     179     180     176     191     195+EL   1S
+EL   1S     196     192                                                        
CHEXA         65       1     176     180     147     146     192     196+EL   1T
+EL   1T     143     142                                                        
CHEXA         66       1      69      73     181     177      68      72+EL   1U
+EL   1U     197     193                                                        
CHEXA         67       1     177     181     182     178     193     197+EL   1V
+EL   1V     198     194                                                        
CHEXA         68       1     178     182     183     179     194     198+EL   1W
+EL   1W     199     195                                                        
CHEXA         69       1     179     183     184     180     195     199+EL   1X
+EL   1X     200     196                                                        
CHEXA         70       1     180     184     148     147     196     200+EL   1Y
+EL   1Y     144     143                                                        
CHEXA         71       1      73      57      90     181      72      77+EL   1Z
+EL   1Z      89     197                                                        
CHEXA         72       1     181      90      94     182     197      89+EL   20
+EL   20      93     198                                                        
CHEXA         73       1     182      94      98     183     198      93+EL   21
+EL   21      97     199                                                        
CHEXA         74       1     183      98     102     184     199      97+EL   22
+EL   22     101     200                                                        
CHEXA         75       1     184     102      83     148     200     101+EL   23
+EL   23      84     144                                                        
CHEXA         76       1      51      60     185      47      30      59+EL   24
+EL   24     201      46                                                        
CHEXA         77       1      47     185     186      43      46     201+EL   25
+EL   25     202      42                                                        
CHEXA         78       1      43     186     187      39      42     202+EL   26
+EL   26     203      38                                                        
CHEXA         79       1      39     187     188      35      38     203+EL   27
+EL   27     204      34                                                        
CHEXA         80       1      35     188     141      25      34     204+EL   28
+EL   28     137     136                                                        
CHEXA         81       1      60      64     189     185      59      63+EL   29
+EL   29     205     201                                                        
CHEXA         82       1     185     189     190     186     201     205+EL   2A
+EL   2A     206     202                                                        
CHEXA         83       1     186     190     191     187     202     206+EL   2B
+EL   2B     207     203                                                        
CHEXA         84       1     187     191     192     188     203     207+EL   2C
+EL   2C     208     204                                                        
CHEXA         85       1     188     192     142     141     204     208+EL   2D
+EL   2D     138     137                                                        
CHEXA         86       1      64      68     193     189      63      67+EL   2E
+EL   2E     209     205                                                        
CHEXA         87       1     189     193     194     190     205     209+EL   2F
+EL   2F     210     206                                                        
CHEXA         88       1     190     194     195     191     206     210+EL   2G
+EL   2G     211     207                                                        
CHEXA         89       1     191     195     196     192     207     211+EL   2H
+EL   2H     212     208                                                        
CHEXA         90       1     192     196     143     142     208     212+EL   2I
+EL   2I     139     138                                                        
CHEXA         91       1      68      72     197     193      67      71+EL   2J
+EL   2J     213     209                                                        
CHEXA         92       1     193     197     198     194     209     213+EL   2K
+EL   2K     214     210                                                        
CHEXA         93       1     194     198     199     195     210     214+EL   2L
+EL   2L     215     211                                                        
CHEXA         94       1     195     199     200     196     211     215+EL   2M
+EL   2M     216     212                                                        
CHEXA         95       1     196     200     144     143     212     216+EL   2N
+EL   2N     140     139                                                        
CHEXA         96       1      72      77      89     197      71      76+EL   2O
+EL   2O      88     213                                                        
CHEXA         97       1     197      89      93     198     213      88+EL   2P
+EL   2P      92     214                                                        
CHEXA         98       1     198      93      97     199     214      92+EL   2Q
+EL   2Q      96     215                                                        
CHEXA         99       1     199      97     101     200     215      96+EL   2R
+EL   2R     100     216                                                        
CHEXA        100       1     200     101      84     144     216     100+EL   2S
+EL   2S     131     140                                                        
CHEXA        101       1      30      59     201      46      50     108+EL   2T
+EL   2T     128      31                                                        
CHEXA        102       1      46     201     202      42      31     128+EL   2U
+EL   2U     124      32                                                        
CHEXA        103       1      42     202     203      38      32     124+EL   2V
+EL   2V     120     107                                                        
CHEXA        104       1      38     203     204      34     107     120+EL   2W
+EL   2W     116      33                                                        
CHEXA        105       1      34     204     137     136      33     116+EL   2X
+EL   2X     105     106                                                        
CHEXA        106       1      59      63     205     201     108     109+EL   2Y
+EL   2Y     127     128                                                        
CHEXA        107       1     201     205     206     202     128     127+EL   2Z
+EL   2Z     123     124                                                        
CHEXA        108       1     202     206     207     203     124     123+EL   30
+EL   30     119     120                                                        
CHEXA        109       1     203     207     208     204     120     119+EL   31
+EL   31     115     116                                                        
CHEXA        110       1     204     208     138     137     116     115+EL   32
+EL   32     129     105                                                        
CHEXA        111       1      63      67     209     205     109      58+EL   33
+EL   33     126     127                                                        
CHEXA        112       1     205     209     210     206     127     126+EL   34
+EL   34     122     123                                                        
CHEXA        113       1     206     210     211     207     123     122+EL   35
+EL   35     118     119                                                        
CHEXA        114       1     207     211     212     208     119     118+EL   36
+EL   36     114     115                                                        
CHEXA        115       1     208     212     139     138     115     114+EL   37
+EL   37     104     129                                                        
CHEXA        116       1      67      71     213     209      58     110+EL   38
+EL   38     125     126                                                        
CHEXA        117       1     209     213     214     210     126     125+EL   39
+EL   39     121     122                                                        
CHEXA        118       1     210     214     215     211     122     121+EL   3A
+EL   3A     117     118                                                        
CHEXA        119       1     211     215     216     212     118     117+EL   3B
+EL   3B     113     114                                                        
CHEXA        120       1     212     216     140     139     114     113+EL   3C
+EL   3C     130     104                                                        
CHEXA        121       1      71      76      88     213     110      75+EL   3D
+EL   3D      87     125                                                        
CHEXA        122       1     213      88      92     214     125      87+EL   3E
+EL   3E     111     121                                                        
CHEXA        123       1     214      92      96     215     121     111+EL   3F
+EL   3F     112     117                                                        
CHEXA        124       1     215      96     100     216     117     112+EL   3G
+EL   3G      86     113                                                        
CHEXA        125       1     216     100     131     140     113      86+EL   3H
+EL   3H      85     130                                                        
ENDDATA