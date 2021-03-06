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
$   Date       : Wed Feb 10 23:03:45 2010
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
$ FEMAP Load Set 111 : Inclusion BSOURCE 0.5 0.5 0.25 0.25 0.
$ FEMAP Property 1 : Untitled
PSHELL         1       1      0.       1      1.       1 0.83333      0.
$ FEMAP Property 2 : Untitled
PSHELL         2       1      0.       1      1.       1 0.83333      0.
$ FEMAP Material 1 : Untitled
MAT1           1     7.8     2.3     0.3      0.      0.      0.        
GRID           1       0     0.5     0.5      0.       0        
GRID           2       0     0.5 0.58333      0.       0        
GRID           3       0     0.5    0.75      0.       0        
GRID           4       0 0.32322 0.67678      0.       0        
GRID           5       0 0.42242 0.65335      0.       0        
GRID           6       0 0.42122 0.57844      0.       0        
GRID           7       0 0.34693 0.57786      0.       0        
GRID           8       0 0.36382 0.63573      0.       0        
GRID           9       0      0.      1.      0.       0        
GRID          10       0      0.     0.9      0.       0        
GRID          11       0      0.     0.8      0.       0        
GRID          12       0      0.     0.7      0.       0        
GRID          13       0      0.     0.6      0.       0        
GRID          14       0 0.16667     0.5      0.       0        
GRID          15       0    0.25     0.5      0.       0        
GRID          16       0 0.26903 0.59567      0.       0        
GRID          17       0 0.40433 0.73097      0.       0        
GRID          18       0     0.4      1.      0.       0        
GRID          19       0     0.3      1.      0.       0        
GRID          20       0     0.2      1.      0.       0        
GRID          21       0     0.1      1.      0.       0        
GRID          22       0 0.39203 0.81257      0.       0        
GRID          23       0 0.30428 0.81746      0.       0        
GRID          24       0  0.2162   0.786      0.       0        
GRID          25       0 0.10707 0.77368      0.       0        
GRID          26       0 0.18026 0.60406      0.       0        
GRID          27       0 0.12491 0.64751      0.       0        
GRID          28       00.072436 0.68198      0.       0        
GRID          29       0 0.10346 0.59581      0.       0        
GRID          30       00.092168 0.54796      0.       0        
GRID          31       0 0.13278 0.54291      0.       0        
GRID          32       0 0.13859 0.58146      0.       0        
GRID          33       00.056232 0.60776      0.       0        
GRID          34       00.049721 0.55171      0.       0        
GRID          35       0 0.20019 0.69689      0.       0        
GRID          36       0 0.30789 0.76353      0.       0        
GRID          37       0 0.31187 0.71791      0.       0        
GRID          38       0 0.35632 0.74041      0.       0        
GRID          39       0 0.35179 0.77266      0.       0        
GRID          40       0 0.26108 0.68914      0.       0        
GRID          41       0 0.26132 0.74693      0.       0        
GRID          42       0 0.14467 0.70567      0.       0        
GRID          43       0 0.20525 0.89634      0.       0        
GRID          44       0 0.39848 0.90887      0.       0        
GRID          45       0 0.30228 0.90556      0.       0        
GRID          46       0 0.10329  0.8924      0.       0        
GRID          47       0 0.66667     0.5      0.       0        
GRID          48       0     0.5 0.66667      0.       0        
GRID          49       0 0.67678 0.67678      0.       0        
GRID          50       0 0.57758 0.65335      0.       0        
GRID          51       0 0.57878 0.57844      0.       0        
GRID          52       0 0.65307 0.57786      0.       0        
GRID          53       0 0.63618 0.63573      0.       0        
GRID          54       0      1.      1.      0.       0        
GRID          55       0      1.     0.9      0.       0        
GRID          56       0      1.     0.8      0.       0        
GRID          57       0      1.     0.7      0.       0        
GRID          58       0      1.     0.6      0.       0        
GRID          59       0      1.     0.5      0.       0        
GRID          60       0 0.83333     0.5      0.       0        
GRID          61       0    0.75     0.5      0.       0        
GRID          62       0 0.73097 0.59567      0.       0        
GRID          63       0 0.59567 0.73097      0.       0        
GRID          64       0     0.5 0.83333      0.       0        
GRID          65       0     0.5 0.91667      0.       0        
GRID          66       0     0.5      1.      0.       0        
GRID          67       0     0.6      1.      0.       0        
GRID          68       0     0.7      1.      0.       0        
GRID          69       0     0.8      1.      0.       0        
GRID          70       0     0.9      1.      0.       0        
GRID          71       0 0.60797 0.81257      0.       0        
GRID          72       0 0.69572 0.81746      0.       0        
GRID          73       0  0.7838   0.786      0.       0        
GRID          74       0 0.89293 0.77368      0.       0        
GRID          75       0 0.81974 0.60406      0.       0        
GRID          76       0 0.87509 0.64751      0.       0        
GRID          77       0 0.92756 0.68198      0.       0        
GRID          78       0 0.89654 0.59581      0.       0        
GRID          79       0 0.90783 0.54796      0.       0        
GRID          80       0 0.86722 0.54291      0.       0        
GRID          81       0 0.86141 0.58146      0.       0        
GRID          82       0 0.94377 0.60776      0.       0        
GRID          83       0 0.95028 0.55171      0.       0        
GRID          84       0 0.79981 0.69689      0.       0        
GRID          85       0 0.69211 0.76353      0.       0        
GRID          86       0 0.68813 0.71791      0.       0        
GRID          87       0 0.64368 0.74041      0.       0        
GRID          88       0 0.64821 0.77266      0.       0        
GRID          89       0 0.73892 0.68914      0.       0        
GRID          90       0 0.73868 0.74693      0.       0        
GRID          91       0 0.85533 0.70567      0.       0        
GRID          92       0 0.79475 0.89634      0.       0        
GRID          93       0 0.60152 0.90887      0.       0        
GRID          94       0 0.69772 0.90556      0.       0        
GRID          95       0 0.89671  0.8924      0.       0        
GRID          96       0 0.33333     0.5      0.       0        
GRID          97       0 0.41667     0.5      0.       0        
GRID          98       0 0.32322 0.32322      0.       0        
GRID          99       0 0.26903 0.40433      0.       0        
GRID         100       0 0.42242 0.34665      0.       0        
GRID         101       0 0.42122 0.42156      0.       0        
GRID         102       0 0.34693 0.42214      0.       0        
GRID         103       0 0.36382 0.36427      0.       0        
GRID         104       0      0.      0.      0.       0        
GRID         105       0      0.     0.1      0.       0        
GRID         106       0      0.     0.2      0.       0        
GRID         107       0      0.     0.3      0.       0        
GRID         108       0      0.     0.4      0.       0        
GRID         109       0      0.     0.5      0.       0        
GRID         110       00.083333     0.5      0.       0        
GRID         111       0 0.40433 0.26903      0.       0        
GRID         112       0     0.5 0.16667      0.       0        
GRID         113       0     0.50.083333      0.       0        
GRID         114       0     0.5      0.      0.       0        
GRID         115       0     0.4      0.      0.       0        
GRID         116       0     0.3      0.      0.       0        
GRID         117       0     0.2      0.      0.       0        
GRID         118       0     0.1      0.      0.       0        
GRID         119       0 0.39203 0.18743      0.       0        
GRID         120       0 0.30428 0.18254      0.       0        
GRID         121       0  0.2162   0.214      0.       0        
GRID         122       0 0.10707 0.22632      0.       0        
GRID         123       0 0.18026 0.39594      0.       0        
GRID         124       0 0.12491 0.35249      0.       0        
GRID         125       00.072436 0.31802      0.       0        
GRID         126       0 0.10346 0.40419      0.       0        
GRID         127       00.092168 0.45204      0.       0        
GRID         128       0 0.13278 0.45709      0.       0        
GRID         129       0 0.13859 0.41854      0.       0        
GRID         130       00.056232 0.39224      0.       0        
GRID         131       00.049721 0.44829      0.       0        
GRID         132       0 0.20019 0.30311      0.       0        
GRID         133       0 0.30789 0.23647      0.       0        
GRID         134       0 0.31187 0.28209      0.       0        
GRID         135       0 0.35632 0.25959      0.       0        
GRID         136       0 0.35179 0.22734      0.       0        
GRID         137       0 0.26108 0.31086      0.       0        
GRID         138       0 0.26132 0.25307      0.       0        
GRID         139       0 0.14467 0.29433      0.       0        
GRID         140       0 0.20525 0.10366      0.       0        
GRID         141       0 0.398480.091129      0.       0        
GRID         142       0 0.30228 0.09444      0.       0        
GRID         143       0 0.10329  0.1076      0.       0        
GRID         144       0 0.58333     0.5      0.       0        
GRID         145       0     0.5 0.41667      0.       0        
GRID         146       0     0.5 0.33333      0.       0        
GRID         147       0 0.57758 0.34665      0.       0        
GRID         148       0 0.57878 0.42156      0.       0        
GRID         149       0 0.65307 0.42214      0.       0        
GRID         150       0 0.63618 0.36427      0.       0        
GRID         151       0      1.      0.      0.       0        
GRID         152       0      1.     0.1      0.       0        
GRID         153       0      1.     0.2      0.       0        
GRID         154       0      1.     0.3      0.       0        
GRID         155       0      1.     0.4      0.       0        
GRID         156       0 0.91667     0.5      0.       0        
GRID         157       0 0.73097 0.40433      0.       0        
GRID         158       0 0.67678 0.32322      0.       0        
GRID         159       0 0.59567 0.26903      0.       0        
GRID         160       0     0.5    0.25      0.       0        
GRID         161       0     0.6      0.      0.       0        
GRID         162       0     0.7      0.      0.       0        
GRID         163       0     0.8      0.      0.       0        
GRID         164       0     0.9      0.      0.       0        
GRID         165       0 0.60797 0.18743      0.       0        
GRID         166       0 0.69572 0.18254      0.       0        
GRID         167       0  0.7838   0.214      0.       0        
GRID         168       0 0.89293 0.22632      0.       0        
GRID         169       0 0.81974 0.39594      0.       0        
GRID         170       0 0.87509 0.35249      0.       0        
GRID         171       0 0.92756 0.31802      0.       0        
GRID         172       0 0.89654 0.40419      0.       0        
GRID         173       0 0.90783 0.45204      0.       0        
GRID         174       0 0.86722 0.45709      0.       0        
GRID         175       0 0.86141 0.41854      0.       0        
GRID         176       0 0.94377 0.39224      0.       0        
GRID         177       0 0.95028 0.44829      0.       0        
GRID         178       0 0.79981 0.30311      0.       0        
GRID         179       0 0.69211 0.23647      0.       0        
GRID         180       0 0.68813 0.28209      0.       0        
GRID         181       0 0.64368 0.25959      0.       0        
GRID         182       0 0.64821 0.22734      0.       0        
GRID         183       0 0.73892 0.31086      0.       0        
GRID         184       0 0.73868 0.25307      0.       0        
GRID         185       0 0.85533 0.29433      0.       0        
GRID         186       0 0.79475 0.10366      0.       0        
GRID         187       0 0.601520.091129      0.       0        
GRID         188       0 0.69772 0.09444      0.       0        
GRID         189       0 0.89671  0.1076      0.       0        
CQUAD4         1       2      48       3      17       5                
CQUAD4         2       2       2      48       5       6                
CQUAD4         3       2      97       1       2       6                
CQUAD4         4       2       5      17       4       8                
CQUAD4         5       2       8       4      16       7                
CQUAD4         6       2       6       5       8       7                
CQUAD4         7       2      96      97       6       7                
CQUAD4         8       2      15      96       7      16                
CQUAD4         9       1      14      15      16      26                
CQUAD4        10       1      26      27      29      32                
CQUAD4        11       1      32      29      30      31                
CQUAD4        12       1      14      26      32      31                
CQUAD4        13       1     110      14      31      30                
CQUAD4        14       1      30      29      33      34                
CQUAD4        15       1     109     110      30      34                
CQUAD4        16       1      13     109      34      33                
CQUAD4        17       1      33      29      27      28                
CQUAD4        18       1      12      13      33      28                
CQUAD4        19       1      17       3      64      22                
CQUAD4        20       1      22      23      36      39                
CQUAD4        21       1      39      36      37      38                
CQUAD4        22       1      17      22      39      38                
CQUAD4        23       1       4      17      38      37                
CQUAD4        24       1      24      35      40      41                
CQUAD4        25       1      36      23      24      41                
CQUAD4        26       1      37      36      41      40                
CQUAD4        27       1      16       4      37      40                
CQUAD4        28       1      26      16      40      35                
CQUAD4        29       1      27      26      35      42                
CQUAD4        30       1      42      35      24      25                
CQUAD4        31       1      28      27      42      25                
CQUAD4        32       1      11      12      28      25                
CQUAD4        33       1      65      66      18      44                
CQUAD4        34       1      22      64      65      44                
CQUAD4        35       1      44      18      19      45                
CQUAD4        36       1      23      22      44      45                
CQUAD4        37       1      45      19      20      43                
CQUAD4        38       1      24      23      45      43                
CQUAD4        39       1      43      20      21      46                
CQUAD4        40       1      25      24      43      46                
CQUAD4        41       1      10      11      25      46                
CQUAD4        42       1       9      10      46      21                
CQUAD4        43       2       3      48      50      63                
CQUAD4        44       2      48       2      51      50                
CQUAD4        45       2       1     144      51       2                
CQUAD4        46       2      63      50      53      49                
CQUAD4        47       2      49      53      52      62                
CQUAD4        48       2      50      51      52      53                
CQUAD4        49       2     144      47      52      51                
CQUAD4        50       2      47      61      62      52                
CQUAD4        51       1      61      60      75      62                
CQUAD4        52       1      76      75      81      78                
CQUAD4        53       1      78      81      80      79                
CQUAD4        54       1      75      60      80      81                
CQUAD4        55       1      60     156      79      80                
CQUAD4        56       1      78      79      83      82                
CQUAD4        57       1     156      59      83      79                
CQUAD4        58       1      59      58      82      83                
CQUAD4        59       1      78      82      77      76                
CQUAD4        60       1      58      57      77      82                
CQUAD4        61       1       3      63      71      64                
CQUAD4        62       1      72      71      88      85                
CQUAD4        63       1      85      88      87      86                
CQUAD4        64       1      71      63      87      88                
CQUAD4        65       1      63      49      86      87                
CQUAD4        66       1      84      73      90      89                
CQUAD4        67       1      72      85      90      73                
CQUAD4        68       1      85      86      89      90                
CQUAD4        69       1      49      62      89      86                
CQUAD4        70       1      62      75      84      89                
CQUAD4        71       1      75      76      91      84                
CQUAD4        72       1      84      91      74      73                
CQUAD4        73       1      76      77      74      91                
CQUAD4        74       1      57      56      74      77                
CQUAD4        75       1      66      65      93      67                
CQUAD4        76       1      64      71      93      65                
CQUAD4        77       1      67      93      94      68                
CQUAD4        78       1      71      72      94      93                
CQUAD4        79       1      68      94      92      69                
CQUAD4        80       1      72      73      92      94                
CQUAD4        81       1      69      92      95      70                
CQUAD4        82       1      73      74      95      92                
CQUAD4        83       1      56      55      95      74                
CQUAD4        84       1      55      54      70      95                
CQUAD4        85       2     160     146     100     111                
CQUAD4        86       2     146     145     101     100                
CQUAD4        87       2       1      97     101     145                
CQUAD4        88       2     111     100     103      98                
CQUAD4        89       2      98     103     102      99                
CQUAD4        90       2     100     101     102     103                
CQUAD4        91       2      97      96     102     101                
CQUAD4        92       2      96      15      99     102                
CQUAD4        93       1      15      14     123      99                
CQUAD4        94       1     124     123     129     126                
CQUAD4        95       1     126     129     128     127                
CQUAD4        96       1     123      14     128     129                
CQUAD4        97       1      14     110     127     128                
CQUAD4        98       1     126     127     131     130                
CQUAD4        99       1     110     109     131     127                
CQUAD4       100       1     109     108     130     131                
CQUAD4       101       1     126     130     125     124                
CQUAD4       102       1     108     107     125     130                
CQUAD4       103       1     160     111     119     112                
CQUAD4       104       1     120     119     136     133                
CQUAD4       105       1     133     136     135     134                
CQUAD4       106       1     119     111     135     136                
CQUAD4       107       1     111      98     134     135                
CQUAD4       108       1     132     121     138     137                
CQUAD4       109       1     120     133     138     121                
CQUAD4       110       1     133     134     137     138                
CQUAD4       111       1      98      99     137     134                
CQUAD4       112       1      99     123     132     137                
CQUAD4       113       1     123     124     139     132                
CQUAD4       114       1     132     139     122     121                
CQUAD4       115       1     124     125     122     139                
CQUAD4       116       1     107     106     122     125                
CQUAD4       117       1     114     113     141     115                
CQUAD4       118       1     112     119     141     113                
CQUAD4       119       1     115     141     142     116                
CQUAD4       120       1     119     120     142     141                
CQUAD4       121       1     116     142     140     117                
CQUAD4       122       1     120     121     140     142                
CQUAD4       123       1     117     140     143     118                
CQUAD4       124       1     121     122     143     140                
CQUAD4       125       1     106     105     143     122                
CQUAD4       126       1     105     104     118     143                
CQUAD4       127       2     146     160     159     147                
CQUAD4       128       2     145     146     147     148                
CQUAD4       129       2     144       1     145     148                
CQUAD4       130       2     147     159     158     150                
CQUAD4       131       2     150     158     157     149                
CQUAD4       132       2     148     147     150     149                
CQUAD4       133       2      47     144     148     149                
CQUAD4       134       2      61      47     149     157                
CQUAD4       135       1      60      61     157     169                
CQUAD4       136       1     169     170     172     175                
CQUAD4       137       1     175     172     173     174                
CQUAD4       138       1      60     169     175     174                
CQUAD4       139       1     156      60     174     173                
CQUAD4       140       1     173     172     176     177                
CQUAD4       141       1      59     156     173     177                
CQUAD4       142       1     155      59     177     176                
CQUAD4       143       1     176     172     170     171                
CQUAD4       144       1     154     155     176     171                
CQUAD4       145       1     159     160     112     165                
CQUAD4       146       1     165     166     179     182                
CQUAD4       147       1     182     179     180     181                
CQUAD4       148       1     159     165     182     181                
CQUAD4       149       1     158     159     181     180                
CQUAD4       150       1     167     178     183     184                
CQUAD4       151       1     179     166     167     184                
CQUAD4       152       1     180     179     184     183                
CQUAD4       153       1     157     158     180     183                
CQUAD4       154       1     169     157     183     178                
CQUAD4       155       1     170     169     178     185                
CQUAD4       156       1     185     178     167     168                
CQUAD4       157       1     171     170     185     168                
CQUAD4       158       1     153     154     171     168                
CQUAD4       159       1     113     114     161     187                
CQUAD4       160       1     165     112     113     187                
CQUAD4       161       1     187     161     162     188                
CQUAD4       162       1     166     165     187     188                
CQUAD4       163       1     188     162     163     186                
CQUAD4       164       1     167     166     188     186                
CQUAD4       165       1     186     163     164     189                
CQUAD4       166       1     168     167     186     189                
CQUAD4       167       1     152     153     168     189                
CQUAD4       168       1     151     152     189     164                
ENDDATA
