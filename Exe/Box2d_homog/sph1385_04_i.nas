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
$   Date       : Wed Jun 13 20:49:12 2007
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
$ FEMAP Load Set 111 : Inclusion BSOURCE 1.50965 1.50965 1.385 1.385 0.
$ FEMAP Property 1 : Untitled
PSHELL         1       1      0.       1      1.       1 0.83333      0.
$ FEMAP Property 2 : Untitled
PSHELL         2       1      0.       1      1.       1 0.83333      0.
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           2       0 0.13321 1.87847      0.       0        
GRID           6       0 1.14083 2.88609      0.       0        
GRID           8       0 1.50965  3.0193      0.       0        
GRID           9       0 1.13224  3.0193      0.       0        
GRID          10       0 0.75482  3.0193      0.       0        
GRID          11       0 0.37741  3.0193      0.       0        
GRID          12       0      0.  3.0193      0.       0        
GRID          13       0      0. 2.64189      0.       0        
GRID          14       0      0. 2.26448      0.       0        
GRID          15       0      0. 1.88706      0.       0        
GRID          17       0 1.50965 2.89465      0.       0        
GRID          20       0 0.79715 2.74374      0.       0        
GRID          21       0 0.50202 2.51728      0.       0        
GRID          22       0 0.27556 2.22215      0.       0        
GRID          26       0 0.17184 1.86811      0.       0        
GRID          27       0  0.3102 2.20215      0.       0        
GRID          29       0 0.81715  2.7091      0.       0        
GRID          30       0 1.15119 2.84746      0.       0        
GRID          32       0 1.50965  1.8559      0.       0        
GRID          38       0 0.53031 2.48899      0.       0        
GRID          42       0  0.4709 1.50965      0.       0        
GRID          45       0  1.1768 1.84215      0.       0        
GRID          46       0 1.18722 2.17089      0.       0        
GRID          47       0 1.18346 2.49774      0.       0        
GRID          48       0 0.84853  1.8322      0.       0        
GRID          49       0 0.52119 1.83546      0.       0        
GRID          50       0 0.87875  2.1398      0.       0        
GRID          51       0  0.5952 2.13193      0.       0        
GRID          52       0  0.8875 2.42424      0.       0        
GRID          53       0  0.6704 2.34837      0.       0        
GRID          56       0 0.27556 0.79715      0.       0        
GRID          57       0 0.50202 0.50202      0.       0        
GRID          60       0 1.50965 0.08465      0.       0        
GRID          61       0 1.50965      0.      0.       0        
GRID          62       0 1.13224      0.      0.       0        
GRID          63       0 0.75482      0.      0.       0        
GRID          64       0 0.37741      0.      0.       0        
GRID          65       0      0.      0.      0.       0        
GRID          66       0      0. 0.37741      0.       0        
GRID          67       0      0. 0.75482      0.       0        
GRID          68       0      0. 1.13224      0.       0        
GRID          69       0      0. 1.50965      0.       0        
GRID          72       0 1.14083 0.13321      0.       0        
GRID          73       0 0.79715 0.27556      0.       0        
GRID          76       0 0.13321 1.14083      0.       0        
GRID          77       0 0.08465 1.50965      0.       0        
GRID          79       0 0.17184 1.15119      0.       0        
GRID          82       0 0.81715  0.3102      0.       0        
GRID          83       0 1.15119 0.17184      0.       0        
GRID          84       0 1.50965 1.50965      0.       0        
GRID          85       0 1.50965  1.1634      0.       0        
GRID          86       0 1.50965 0.81715      0.       0        
GRID          87       0 1.50965  0.4709      0.       0        
GRID          88       0 1.50965 0.12465      0.       0        
GRID          91       0 0.53031 0.53031      0.       0        
GRID          92       0  0.3102 0.81715      0.       0        
GRID          94       0 0.12465 1.50965      0.       0        
GRID          96       0 0.81715 1.50965      0.       0        
GRID          97       0  1.1634 1.50965      0.       0        
GRID          98       0  1.1768 1.17715      0.       0        
GRID          99       0 1.18722 0.84841      0.       0        
GRID         100       0 1.18346 0.52156      0.       0        
GRID         101       0 0.84853  1.1871      0.       0        
GRID         102       0 0.52119 1.18384      0.       0        
GRID         103       0 0.87875  0.8795      0.       0        
GRID         104       0  0.5952 0.88737      0.       0        
GRID         105       0  0.8875 0.59506      0.       0        
GRID         106       0  0.6704 0.67093      0.       0        
GRID         110       0 2.51728 2.51728      0.       0        
GRID         111       0 2.22215 2.74374      0.       0        
GRID         115       0 1.88706  3.0193      0.       0        
GRID         116       0 2.26448  3.0193      0.       0        
GRID         117       0 2.64189  3.0193      0.       0        
GRID         118       0  3.0193  3.0193      0.       0        
GRID         119       0  3.0193 2.64189      0.       0        
GRID         120       0  3.0193 2.26448      0.       0        
GRID         121       0  3.0193 1.88706      0.       0        
GRID         122       0  3.0193 1.50965      0.       0        
GRID         124       0 1.50965 2.93465      0.       0        
GRID         125       0 1.87847 2.88609      0.       0        
GRID         128       0 2.74374 2.22215      0.       0        
GRID         129       0 2.88609 1.87847      0.       0        
GRID         133       0  2.7091 2.20215      0.       0        
GRID         135       0 2.20215  2.7091      0.       0        
GRID         139       0 1.50965 2.20215      0.       0        
GRID         140       0 1.50965  2.5484      0.       0        
GRID         142       0 1.86811 2.84746      0.       0        
GRID         144       0 2.48899 2.48899      0.       0        
GRID         146       0 2.84746 1.86811      0.       0        
GRID         147       0 2.89465 1.50965      0.       0        
GRID         148       0  2.5484 1.50965      0.       0        
GRID         149       0 2.20215 1.50965      0.       0        
GRID         150       0  1.8559 1.50965      0.       0        
GRID         151       0  1.8425 1.84215      0.       0        
GRID         152       0 1.83208 2.17089      0.       0        
GRID         153       0 1.83584 2.49774      0.       0        
GRID         154       0 2.17077  1.8322      0.       0        
GRID         155       0 2.49811 1.83546      0.       0        
GRID         156       0 2.14055  2.1398      0.       0        
GRID         157       0  2.4241 2.13193      0.       0        
GRID         158       0  2.1318 2.42424      0.       0        
GRID         159       0  2.3489 2.34837      0.       0        
GRID         164       0 2.22215 0.27556      0.       0        
GRID         165       0 1.87847 0.13321      0.       0        
GRID         168       0 1.88706      0.      0.       0        
GRID         169       0 2.26448      0.      0.       0        
GRID         170       0 2.64189      0.      0.       0        
GRID         171       0  3.0193      0.      0.       0        
GRID         172       0  3.0193 0.37741      0.       0        
GRID         173       0  3.0193 0.75482      0.       0        
GRID         174       0  3.0193 1.13224      0.       0        
GRID         180       0 2.51728 0.50202      0.       0        
GRID         181       0 2.74374 0.79715      0.       0        
GRID         182       0 2.88609 1.14083      0.       0        
GRID         183       0 2.93465 1.50965      0.       0        
GRID         186       0  2.7091 0.81715      0.       0        
GRID         188       0 2.20215  0.3102      0.       0        
GRID         195       0 1.86811 0.17184      0.       0        
GRID         197       0 2.48899 0.53031      0.       0        
GRID         199       0 2.84746 1.15119      0.       0        
GRID         204       0  1.8425 1.17715      0.       0        
GRID         205       0 1.83208 0.84841      0.       0        
GRID         206       0 1.83584 0.52156      0.       0        
GRID         207       0 2.17077  1.1871      0.       0        
GRID         208       0 2.49811 1.18384      0.       0        
GRID         209       0 2.14055  0.8795      0.       0        
GRID         210       0  2.4241 0.88737      0.       0        
GRID         211       0  2.1318 0.59506      0.       0        
GRID         212       0  2.3489 0.67093      0.       0        
CQUAD4         1       1       6     124       8       9                
CQUAD4         2       1      20       6       9      10                
CQUAD4         3       1      21      20      10      11                
CQUAD4         5       1      22      21      13      14                
CQUAD4         6       1       2      22      14      15                
CQUAD4         7       1      77       2      15      69                
CTRIA3         8       1      13      11      12                        
CTRIA3         9       1      21      11      13                        
CQUAD4        10       1      17     124       6      30                
CQUAD4        11       1      30       6      20      29                
CQUAD4        12       1      29      20      21      38                
CQUAD4        13       1      38      21      22      27                
CQUAD4        14       1      27      22       2      26                
CQUAD4        15       1      26       2      77      94                
CQUAD4        16       2      96      97      45      48                
CQUAD4        17       2      42      96      48      49                
CQUAD4        18       2      26      94      42      49                
CQUAD4        19       2      48      45      46      50                
CQUAD4        20       2      49      48      50      51                
CQUAD4        21       2      27      26      49      51                
CQUAD4        22       2      38      27      51      53                
CQUAD4        23       2      53      51      50      52                
CQUAD4        24       2      29      38      53      52                
CQUAD4        25       2      52      50      46      47                
CQUAD4        26       2      30      29      52      47                
CQUAD4        27       2     140      17      30      47                
CQUAD4        28       2     139     140      47      46                
CQUAD4        29       2      32     139      46      45                
CQUAD4        30       2      84      32      45      97                
CQUAD4        31       1      60      72      62      61                
CQUAD4        32       1      72      73      63      62                
CQUAD4        33       1      73      57      64      63                
CQUAD4        34       1      57      56      67      66                
CQUAD4        35       1      56      76      68      67                
CQUAD4        36       1      76      77      69      68                
CTRIA3        37       1      64      66      65                        
CTRIA3        38       1      57      66      64                        
CQUAD4        39       1      60      88      83      72                
CQUAD4        40       1      72      83      82      73                
CQUAD4        41       1      73      82      91      57                
CQUAD4        42       1      57      91      92      56                
CQUAD4        43       1      56      92      79      76                
CQUAD4        44       1      76      79      94      77                
CQUAD4        45       2      97      96     101      98                
CQUAD4        46       2      96      42     102     101                
CQUAD4        47       2      94      79     102      42                
CQUAD4        48       2      98     101     103      99                
CQUAD4        49       2     101     102     104     103                
CQUAD4        50       2      79      92     104     102                
CQUAD4        51       2      92      91     106     104                
CQUAD4        52       2     104     106     105     103                
CQUAD4        53       2      91      82     105     106                
CQUAD4        54       2     103     105     100      99                
CQUAD4        55       2      82      83     100     105                
CQUAD4        56       2      88      87     100      83                
CQUAD4        57       2      87      86      99     100                
CQUAD4        58       2      86      85      98      99                
CQUAD4        59       2      85      84      97      98                
CQUAD4        60       1     124     125     115       8                
CQUAD4        61       1     125     111     116     115                
CQUAD4        62       1     111     110     117     116                
CQUAD4        63       1     110     128     120     119                
CQUAD4        64       1     128     129     121     120                
CQUAD4        65       1     129     183     122     121                
CTRIA3        66       1     117     119     118                        
CTRIA3        67       1     110     119     117                        
CQUAD4        68       1     124      17     142     125                
CQUAD4        69       1     125     142     135     111                
CQUAD4        70       1     111     135     144     110                
CQUAD4        71       1     110     144     133     128                
CQUAD4        72       1     128     133     146     129                
CQUAD4        73       1     129     146     147     183                
CQUAD4        74       2     150     149     154     151                
CQUAD4        75       2     149     148     155     154                
CQUAD4        76       2     147     146     155     148                
CQUAD4        77       2     151     154     156     152                
CQUAD4        78       2     154     155     157     156                
CQUAD4        79       2     146     133     157     155                
CQUAD4        80       2     133     144     159     157                
CQUAD4        81       2     157     159     158     156                
CQUAD4        82       2     144     135     158     159                
CQUAD4        83       2     156     158     153     152                
CQUAD4        84       2     135     142     153     158                
CQUAD4        85       2      17     140     153     142                
CQUAD4        86       2     140     139     152     153                
CQUAD4        87       2     139      32     151     152                
CQUAD4        88       2      32      84     150     151                
CQUAD4        89       1     165      60      61     168                
CQUAD4        90       1     164     165     168     169                
CQUAD4        91       1     180     164     169     170                
CQUAD4        92       1     181     180     172     173                
CQUAD4        93       1     182     181     173     174                
CQUAD4        94       1     183     182     174     122                
CTRIA3        95       1     172     170     171                        
CTRIA3        96       1     180     170     172                        
CQUAD4        97       1      88      60     165     195                
CQUAD4        98       1     195     165     164     188                
CQUAD4        99       1     188     164     180     197                
CQUAD4       100       1     197     180     181     186                
CQUAD4       101       1     186     181     182     199                
CQUAD4       102       1     199     182     183     147                
CQUAD4       103       2     149     150     204     207                
CQUAD4       104       2     148     149     207     208                
CQUAD4       105       2     199     147     148     208                
CQUAD4       106       2     207     204     205     209                
CQUAD4       107       2     208     207     209     210                
CQUAD4       108       2     186     199     208     210                
CQUAD4       109       2     197     186     210     212                
CQUAD4       110       2     212     210     209     211                
CQUAD4       111       2     188     197     212     211                
CQUAD4       112       2     211     209     205     206                
CQUAD4       113       2     195     188     211     206                
CQUAD4       114       2      87      88     195     206                
CQUAD4       115       2      86      87     206     205                
CQUAD4       116       2      85      86     205     204                
CQUAD4       117       2      84      85     204     150                
ENDDATA