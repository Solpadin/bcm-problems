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
  SPC = 1
  LOAD = 1
BEGIN BULK
$ ***************************************************************************
$   Written by : FEMAP
$   Version    : 8.00
$   Translator : NE/Nastran
$   From Model : 
$   Date       : Fri Feb 27 18:20:49 2004
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
$ FEMAP Load Set 1 : L
FORCE          1     241       0      1.   1000.      0.      0.
FORCE          1     242       0      1.   1000.      0.      0.
FORCE          1     243       0      1.   1000.      0.      0.
FORCE          1     244       0      1.   1000.      0.      0.
$ FEMAP Constraint Set 1 : NASTRAN SPC 1
SPC            1       1     123      0.
SPC            1       2     123      0.
SPC            1       3     123      0.
SPC            1       4     123      0.
$ FEMAP Property 1 : Untitled
PSOLID         1       1       0        
$ FEMAP Material 1 : Untitled
MAT1           1 7.2E+10 6.4E+10     0.3      0.      0.      0.        
$ FEMAP Material 2 : Untitled
MAT1           2 7.2E+10 6.4E+10     0.3      0.      0.      0.        
GRID           1       0      0.      0.      0.       0        
GRID           2       0      0.      0.      1.       0        
GRID           3       0      1.      0.      1.       0        
GRID           4       0      1.      0.      0.       0        
GRID           5       0      0.      1.      0.       0        
GRID           6       0      0.      1.      1.       0        
GRID           7       0      1.      1.      1.       0        
GRID           8       0      1.      1.      0.       0        
GRID          13       0      0.      2.      0.       0        
GRID          14       0      0.      2.      1.       0        
GRID          15       0      1.      2.      1.       0        
GRID          16       0      1.      2.      0.       0        
GRID          17       0      0.      3.      0.       0        
GRID          18       0      0.      3.      1.       0        
GRID          19       0      1.      3.      1.       0        
GRID          20       0      1.      3.      0.       0        
GRID          24       0      1.      4.      0.       0        
GRID          25       0      0.      4.      0.       0        
GRID          26       0      0.      4.      1.       0        
GRID          27       0      1.      4.      1.       0        
GRID          29       0      0.      5.      0.       0        
GRID          30       0      0.      5.      1.       0        
GRID          31       0      1.      5.      1.       0        
GRID          32       0      1.      5.      0.       0        
GRID          37       0      0.      6.      0.       0        
GRID          38       0      0.      6.      1.       0        
GRID          39       0      1.      6.      1.       0        
GRID          40       0      1.      6.      0.       0        
GRID          41       0      0.      7.      0.       0        
GRID          42       0      0.      7.      1.       0        
GRID          43       0      1.      7.      1.       0        
GRID          44       0      1.      7.      0.       0        
GRID          49       0      0.      8.      0.       0        
GRID          50       0      0.      8.      1.       0        
GRID          51       0      1.      8.      1.       0        
GRID          52       0      1.      8.      0.       0        
GRID          53       0      0.      9.      0.       0        
GRID          54       0      0.      9.      1.       0        
GRID          55       0      1.      9.      1.       0        
GRID          56       0      1.      9.      0.       0        
GRID          61       0      0.     10.      0.       0        
GRID          62       0      0.     10.      1.       0        
GRID          63       0      1.     10.      1.       0        
GRID          64       0      1.     10.      0.       0        
GRID          65       0      0.     11.      0.       0        
GRID          66       0      0.     11.      1.       0        
GRID          67       0      1.     11.      1.       0        
GRID          68       0      1.     11.      0.       0        
GRID          73       0      0.     12.      0.       0        
GRID          74       0      0.     12.      1.       0        
GRID          75       0      1.     12.      1.       0        
GRID          76       0      1.     12.      0.       0        
GRID          77       0      0.     13.      0.       0        
GRID          78       0      0.     13.      1.       0        
GRID          79       0      1.     13.      1.       0        
GRID          80       0      1.     13.      0.       0        
GRID          85       0      0.     14.      0.       0        
GRID          86       0      0.     14.      1.       0        
GRID          87       0      1.     14.      1.       0        
GRID          88       0      1.     14.      0.       0        
GRID          89       0      0.     15.      0.       0        
GRID          90       0      0.     15.      1.       0        
GRID          91       0      1.     15.      1.       0        
GRID          92       0      1.     15.      0.       0        
GRID          96       0      1.     16.      0.       0        
GRID          97       0      0.     16.      0.       0        
GRID          98       0      0.     16.      1.       0        
GRID          99       0      1.     16.      1.       0        
GRID         101       0      0.     17.      0.       0        
GRID         102       0      0.     17.      1.       0        
GRID         103       0      1.     17.      1.       0        
GRID         104       0      1.     17.      0.       0        
GRID         109       0      0.     18.      0.       0        
GRID         110       0      0.     18.      1.       0        
GRID         111       0      1.     18.      1.       0        
GRID         112       0      1.     18.      0.       0        
GRID         113       0      0.     19.      0.       0        
GRID         114       0      0.     19.      1.       0        
GRID         115       0      1.     19.      1.       0        
GRID         116       0      1.     19.      0.       0        
GRID         121       0      0.     20.      0.       0        
GRID         122       0      0.     20.      1.       0        
GRID         123       0      1.     20.      1.       0        
GRID         124       0      1.     20.      0.       0        
GRID         125       0      0.     21.      0.       0        
GRID         126       0      0.     21.      1.       0        
GRID         127       0      1.     21.      1.       0        
GRID         128       0      1.     21.      0.       0        
GRID         132       0      1.     22.      0.       0        
GRID         133       0      0.     22.      0.       0        
GRID         134       0      0.     22.      1.       0        
GRID         135       0      1.     22.      1.       0        
GRID         137       0      0.     23.      0.       0        
GRID         138       0      0.     23.      1.       0        
GRID         139       0      1.     23.      1.       0        
GRID         140       0      1.     23.      0.       0        
GRID         145       0      0.     24.      0.       0        
GRID         146       0      0.     24.      1.       0        
GRID         147       0      1.     24.      1.       0        
GRID         148       0      1.     24.      0.       0        
GRID         149       0      0.     25.      0.       0        
GRID         150       0      0.     25.      1.       0        
GRID         151       0      1.     25.      1.       0        
GRID         152       0      1.     25.      0.       0        
GRID         157       0      0.     26.      0.       0        
GRID         158       0      0.     26.      1.       0        
GRID         159       0      1.     26.      1.       0        
GRID         160       0      1.     26.      0.       0        
GRID         161       0      0.     27.      0.       0        
GRID         162       0      0.     27.      1.       0        
GRID         163       0      1.     27.      1.       0        
GRID         164       0      1.     27.      0.       0        
GRID         169       0      0.     28.      0.       0        
GRID         170       0      0.     28.      1.       0        
GRID         171       0      1.     28.      1.       0        
GRID         172       0      1.     28.      0.       0        
GRID         173       0      0.     29.      0.       0        
GRID         174       0      0.     29.      1.       0        
GRID         175       0      1.     29.      1.       0        
GRID         176       0      1.     29.      0.       0        
GRID         181       0      0.     30.      0.       0        
GRID         182       0      0.     30.      1.       0        
GRID         183       0      1.     30.      1.       0        
GRID         184       0      1.     30.      0.       0        
GRID         185       0      0.     31.      0.       0        
GRID         186       0      0.     31.      1.       0        
GRID         187       0      1.     31.      1.       0        
GRID         188       0      1.     31.      0.       0        
GRID         193       0      0.     32.      0.       0        
GRID         194       0      0.     32.      1.       0        
GRID         195       0      1.     32.      1.       0        
GRID         196       0      1.     32.      0.       0        
GRID         197       0      0.     33.      0.       0        
GRID         198       0      0.     33.      1.       0        
GRID         199       0      1.     33.      1.       0        
GRID         200       0      1.     33.      0.       0        
GRID         204       0      1.     34.      0.       0        
GRID         205       0      0.     34.      0.       0        
GRID         206       0      0.     34.      1.       0        
GRID         207       0      1.     34.      1.       0        
GRID         209       0      0.     35.      0.       0        
GRID         210       0      0.     35.      1.       0        
GRID         211       0      1.     35.      1.       0        
GRID         212       0      1.     35.      0.       0        
GRID         217       0      0.     36.      0.       0        
GRID         218       0      0.     36.      1.       0        
GRID         219       0      1.     36.      1.       0        
GRID         220       0      1.     36.      0.       0        
GRID         221       0      0.     37.      0.       0        
GRID         222       0      0.     37.      1.       0        
GRID         223       0      1.     37.      1.       0        
GRID         224       0      1.     37.      0.       0        
GRID         229       0      0.     38.      0.       0        
GRID         230       0      0.     38.      1.       0        
GRID         231       0      1.     38.      1.       0        
GRID         232       0      1.     38.      0.       0        
GRID         233       0      0.     39.      0.       0        
GRID         234       0      0.     39.      1.       0        
GRID         235       0      1.     39.      1.       0        
GRID         236       0      1.     39.      0.       0        
GRID         241       0      0.     40.      0.       0        
GRID         242       0      0.     40.      1.       0        
GRID         243       0      1.     40.      1.       0        
GRID         244       0      1.     40.      0.       0        
CTETRA         1       1       1       2       4       5                        
CTETRA         2       1       2       3       4       7                        
CTETRA         3       1       5       2       7       6                        
CTETRA         4       1       4       8       5       7                        
CTETRA         5       1       5       4       7       2                        
CTETRA         6       1       5       6       7      14                        
CTETRA         7       1       5      16      13      14                        
CTETRA         8       1       7      16      14      15                        
CTETRA         9       1      14       5      16       7                        
CTETRA        10       1      16       7       5       8                        
CTETRA        11       1      13      14      16      17                        
CTETRA        12       1      14      15      16      19                        
CTETRA        13       1      17      14      19      18                        
CTETRA        14       1      16      20      17      19                        
CTETRA        15       1      17      16      19      14                        
CTETRA        16       1      17      18      19      26                        
CTETRA        17       1      17      24      25      26                        
CTETRA        18       1      19      24      26      27                        
CTETRA        19       1      26      17      24      19                        
CTETRA        20       1      24      19      17      20                        
CTETRA        21       1      25      26      24      29                        
CTETRA        22       1      26      27      24      31                        
CTETRA        23       1      29      26      31      30                        
CTETRA        24       1      24      32      29      31                        
CTETRA        25       1      29      24      31      26                        
CTETRA        26       1      29      30      31      38                        
CTETRA        27       1      29      40      37      38                        
CTETRA        28       1      31      40      38      39                        
CTETRA        29       1      38      29      40      31                        
CTETRA        30       1      40      31      29      32                        
CTETRA        31       1      37      38      40      41                        
CTETRA        32       1      38      39      40      43                        
CTETRA        33       1      41      38      43      42                        
CTETRA        34       1      40      44      41      43                        
CTETRA        35       1      41      40      43      38                        
CTETRA        36       1      41      42      43      50                        
CTETRA        37       1      41      52      49      50                        
CTETRA        38       1      43      52      50      51                        
CTETRA        39       1      50      41      52      43                        
CTETRA        40       1      52      43      41      44                        
CTETRA        41       1      49      50      52      53                        
CTETRA        42       1      50      51      52      55                        
CTETRA        43       1      53      50      55      54                        
CTETRA        44       1      52      56      53      55                        
CTETRA        45       1      53      52      55      50                        
CTETRA        46       1      53      54      55      62                        
CTETRA        47       1      53      64      61      62                        
CTETRA        48       1      55      64      62      63                        
CTETRA        49       1      62      53      64      55                        
CTETRA        50       1      64      55      53      56                        
CTETRA        51       1      61      62      64      65                        
CTETRA        52       1      62      63      64      67                        
CTETRA        53       1      65      62      67      66                        
CTETRA        54       1      64      68      65      67                        
CTETRA        55       1      65      64      67      62                        
CTETRA        56       1      65      66      67      74                        
CTETRA        57       1      65      76      73      74                        
CTETRA        58       1      67      76      74      75                        
CTETRA        59       1      74      65      76      67                        
CTETRA        60       1      76      67      65      68                        
CTETRA        61       1      73      74      76      77                        
CTETRA        62       1      74      75      76      79                        
CTETRA        63       1      77      74      79      78                        
CTETRA        64       1      76      80      77      79                        
CTETRA        65       1      77      76      79      74                        
CTETRA        66       1      77      78      79      86                        
CTETRA        67       1      77      88      85      86                        
CTETRA        68       1      79      88      86      87                        
CTETRA        69       1      86      77      88      79                        
CTETRA        70       1      88      79      77      80                        
CTETRA        71       1      85      86      88      89                        
CTETRA        72       1      86      87      88      91                        
CTETRA        73       1      89      86      91      90                        
CTETRA        74       1      88      92      89      91                        
CTETRA        75       1      89      88      91      86                        
CTETRA        76       1      89      90      91      98                        
CTETRA        77       1      89      96      97      98                        
CTETRA        78       1      91      96      98      99                        
CTETRA        79       1      98      89      96      91                        
CTETRA        80       1      96      91      89      92                        
CTETRA        81       1      97      98      96     101                        
CTETRA        82       1      98      99      96     103                        
CTETRA        83       1     101      98     103     102                        
CTETRA        84       1      96     104     101     103                        
CTETRA        85       1     101      96     103      98                        
CTETRA        86       1     101     102     103     110                        
CTETRA        87       1     101     112     109     110                        
CTETRA        88       1     103     112     110     111                        
CTETRA        89       1     110     101     112     103                        
CTETRA        90       1     112     103     101     104                        
CTETRA        91       1     109     110     112     113                        
CTETRA        92       1     110     111     112     115                        
CTETRA        93       1     113     110     115     114                        
CTETRA        94       1     112     116     113     115                        
CTETRA        95       1     113     112     115     110                        
CTETRA        96       1     113     114     115     122                        
CTETRA        97       1     113     124     121     122                        
CTETRA        98       1     115     124     122     123                        
CTETRA        99       1     122     113     124     115                        
CTETRA       100       1     124     115     113     116                        
CTETRA       101       1     121     122     124     125                        
CTETRA       102       1     122     123     124     127                        
CTETRA       103       1     125     122     127     126                        
CTETRA       104       1     124     128     125     127                        
CTETRA       105       1     125     124     127     122                        
CTETRA       106       1     125     126     127     134                        
CTETRA       107       1     125     132     133     134                        
CTETRA       108       1     127     132     134     135                        
CTETRA       109       1     134     125     132     127                        
CTETRA       110       1     132     127     125     128                        
CTETRA       111       1     133     134     132     137                        
CTETRA       112       1     134     135     132     139                        
CTETRA       113       1     137     134     139     138                        
CTETRA       114       1     132     140     137     139                        
CTETRA       115       1     137     132     139     134                        
CTETRA       116       1     137     138     139     146                        
CTETRA       117       1     137     148     145     146                        
CTETRA       118       1     139     148     146     147                        
CTETRA       119       1     146     137     148     139                        
CTETRA       120       1     148     139     137     140                        
CTETRA       121       1     145     146     148     149                        
CTETRA       122       1     146     147     148     151                        
CTETRA       123       1     149     146     151     150                        
CTETRA       124       1     148     152     149     151                        
CTETRA       125       1     149     148     151     146                        
CTETRA       126       1     149     150     151     158                        
CTETRA       127       1     149     160     157     158                        
CTETRA       128       1     151     160     158     159                        
CTETRA       129       1     158     149     160     151                        
CTETRA       130       1     160     151     149     152                        
CTETRA       131       1     157     158     160     161                        
CTETRA       132       1     158     159     160     163                        
CTETRA       133       1     161     158     163     162                        
CTETRA       134       1     160     164     161     163                        
CTETRA       135       1     161     160     163     158                        
CTETRA       136       1     161     162     163     170                        
CTETRA       137       1     161     172     169     170                        
CTETRA       138       1     163     172     170     171                        
CTETRA       139       1     170     161     172     163                        
CTETRA       140       1     172     163     161     164                        
CTETRA       141       1     169     170     172     173                        
CTETRA       142       1     170     171     172     175                        
CTETRA       143       1     173     170     175     174                        
CTETRA       144       1     172     176     173     175                        
CTETRA       145       1     173     172     175     170                        
CTETRA       146       1     173     174     175     182                        
CTETRA       147       1     173     184     181     182                        
CTETRA       148       1     175     184     182     183                        
CTETRA       149       1     182     173     184     175                        
CTETRA       150       1     184     175     173     176                        
CTETRA       151       1     181     182     184     185                        
CTETRA       152       1     182     183     184     187                        
CTETRA       153       1     185     182     187     186                        
CTETRA       154       1     184     188     185     187                        
CTETRA       155       1     185     184     187     182                        
CTETRA       156       1     185     186     187     194                        
CTETRA       157       1     185     196     193     194                        
CTETRA       158       1     187     196     194     195                        
CTETRA       159       1     194     185     196     187                        
CTETRA       160       1     196     187     185     188                        
CTETRA       161       1     193     194     196     197                        
CTETRA       162       1     194     195     196     199                        
CTETRA       163       1     197     194     199     198                        
CTETRA       164       1     196     200     197     199                        
CTETRA       165       1     197     196     199     194                        
CTETRA       166       1     197     198     199     206                        
CTETRA       167       1     197     204     205     206                        
CTETRA       168       1     199     204     206     207                        
CTETRA       169       1     206     197     204     199                        
CTETRA       170       1     204     199     197     200                        
CTETRA       171       1     205     206     204     209                        
CTETRA       172       1     206     207     204     211                        
CTETRA       173       1     209     206     211     210                        
CTETRA       174       1     204     212     209     211                        
CTETRA       175       1     209     204     211     206                        
CTETRA       176       1     209     210     211     218                        
CTETRA       177       1     209     220     217     218                        
CTETRA       178       1     211     220     218     219                        
CTETRA       179       1     218     209     220     211                        
CTETRA       180       1     220     211     209     212                        
CTETRA       181       1     217     218     220     221                        
CTETRA       182       1     218     219     220     223                        
CTETRA       183       1     221     218     223     222                        
CTETRA       184       1     220     224     221     223                        
CTETRA       185       1     221     220     223     218                        
CTETRA       186       1     221     222     223     230                        
CTETRA       187       1     221     232     229     230                        
CTETRA       188       1     223     232     230     231                        
CTETRA       189       1     230     221     232     223                        
CTETRA       190       1     232     223     221     224                        
CTETRA       191       1     229     230     232     233                        
CTETRA       192       1     230     231     232     235                        
CTETRA       193       1     233     230     235     234                        
CTETRA       194       1     232     236     233     235                        
CTETRA       195       1     233     232     235     230                        
CTETRA       196       1     233     234     235     242                        
CTETRA       197       1     233     244     241     242                        
CTETRA       198       1     235     244     242     243                        
CTETRA       199       1     242     233     244     235                        
CTETRA       200       1     244     235     233     236                        
ENDDATA
