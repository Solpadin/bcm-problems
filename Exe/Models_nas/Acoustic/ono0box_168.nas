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
$   Date       : Tue Feb 05 21:42:03 2008
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
PSHELL         1       1      0.       1      1.       1 0.83333      0.
$ FEMAP Property 2 : Untitled
PSHELL         2       1      0.       1      1.       1 0.83333      0.
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           5       0 0.76667     0.5      0.       0        
GRID           7       0    0.85 0.33333      0.       0        
GRID           8       0 0.76667    0.25      0.       0        
GRID           9       0 0.68333    0.25      0.       0        
GRID          10       0 0.68333 0.33333      0.       0        
GRID          11       0 0.68333 0.41667      0.       0        
GRID          12       0 0.76667 0.33333      0.       0        
GRID          13       0 0.76667 0.41667      0.       0        
GRID          14       0     1.2      0.      0.       0        
GRID          15       0 1.11429      0.      0.       0        
GRID          16       0 1.02857      0.      0.       0        
GRID          17       0 0.94286      0.      0.       0        
GRID          18       0 0.85714      0.      0.       0        
GRID          19       0 0.77143      0.      0.       0        
GRID          20       0 0.68571      0.      0.       0        
GRID          22       0     0.60.083333      0.       0        
GRID          25       0    0.85    0.25      0.       0        
GRID          26       0    0.85 0.41667      0.       0        
GRID          27       0  0.9375     0.5      0.       0        
GRID          28       0   1.025     0.5      0.       0        
GRID          29       0  1.1125     0.5      0.       0        
GRID          30       0     1.2     0.5      0.       0        
GRID          31       0     1.2 0.41667      0.       0        
GRID          32       0     1.2 0.33333      0.       0        
GRID          33       0     1.2    0.25      0.       0        
GRID          34       0     1.2 0.16667      0.       0        
GRID          35       0     1.20.083333      0.       0        
GRID          36       0  0.9383 0.41633      0.       0        
GRID          37       0 0.93897 0.33276      0.       0        
GRID          38       0 0.93946 0.24946      0.       0        
GRID          39       0 0.94084 0.16628      0.       0        
GRID          40       0 0.941890.083217      0.       0        
GRID          41       0 0.85358 0.16687      0.       0        
GRID          42       0 0.855570.083587      0.       0        
GRID          43       0 0.68412 0.16703      0.       0        
GRID          44       0 0.68491  0.0837      0.       0        
GRID          45       0 0.76853 0.16709      0.       0        
GRID          46       0 0.769970.083793      0.       0        
GRID          47       0 1.11429  0.2496      0.       0        
GRID          48       0 1.02757 0.24926      0.       0        
GRID          49       0 1.11399 0.33295      0.       0        
GRID          50       0 1.02703 0.33267      0.       0        
GRID          51       0 1.11338 0.41643      0.       0        
GRID          52       0 1.02623 0.41627      0.       0        
GRID          53       0 1.02802 0.16616      0.       0        
GRID          54       0 1.028310.083093      0.       0        
GRID          55       0 1.114410.083162      0.       0        
GRID          56       0 1.11438 0.16634      0.       0        
GRID          60       0 0.68333     0.5      0.       0        
GRID          62       0    0.85     0.5      0.       0        
GRID          63       0    0.85 0.66667      0.       0        
GRID          64       0 0.76667    0.75      0.       0        
GRID          65       0 0.68333    0.75      0.       0        
GRID          66       0 0.68333 0.66667      0.       0        
GRID          67       0 0.68333 0.58333      0.       0        
GRID          68       0 0.76667 0.66667      0.       0        
GRID          69       0 0.76667 0.58333      0.       0        
GRID          70       0     1.2      1.      0.       0        
GRID          71       0 1.11429      1.      0.       0        
GRID          72       0 1.02857      1.      0.       0        
GRID          73       0 0.94286      1.      0.       0        
GRID          74       0 0.85714      1.      0.       0        
GRID          75       0 0.77143      1.      0.       0        
GRID          76       0 0.68571      1.      0.       0        
GRID          79       0     0.6 0.83333      0.       0        
GRID          81       0    0.85    0.75      0.       0        
GRID          82       0    0.85 0.58333      0.       0        
GRID          87       0     1.2 0.58333      0.       0        
GRID          88       0     1.2 0.66667      0.       0        
GRID          89       0     1.2    0.75      0.       0        
GRID          90       0     1.2 0.83333      0.       0        
GRID          91       0     1.2 0.91667      0.       0        
GRID          92       0  0.9383 0.58367      0.       0        
GRID          93       0 0.93897 0.66724      0.       0        
GRID          94       0 0.93946 0.75054      0.       0        
GRID          95       0 0.94084 0.83372      0.       0        
GRID          96       0 0.94189 0.91678      0.       0        
GRID          97       0 0.85358 0.83313      0.       0        
GRID          98       0 0.85557 0.91641      0.       0        
GRID          99       0 0.68412 0.83297      0.       0        
GRID         100       0 0.68491  0.9163      0.       0        
GRID         101       0 0.76853 0.83291      0.       0        
GRID         102       0 0.76997 0.91621      0.       0        
GRID         103       0 1.11429  0.7504      0.       0        
GRID         104       0 1.02757 0.75074      0.       0        
GRID         105       0 1.11399 0.66705      0.       0        
GRID         106       0 1.02703 0.66733      0.       0        
GRID         107       0 1.11338 0.58357      0.       0        
GRID         108       0 1.02623 0.58373      0.       0        
GRID         109       0 1.02802 0.83384      0.       0        
GRID         110       0 1.02831 0.91691      0.       0        
GRID         111       0 1.11441 0.91684      0.       0        
GRID         112       0 1.11438 0.83366      0.       0        
GRID         113       0     0.6 0.66667      0.       0        
GRID         114       0     0.6 0.58333      0.       0        
GRID         117       0 0.43333     0.5      0.       0        
GRID         118       0    0.35     0.5      0.       0        
GRID         119       0    0.35 0.66667      0.       0        
GRID         120       0 0.43333    0.75      0.       0        
GRID         121       0 0.51667    0.75      0.       0        
GRID         122       0 0.51667 0.66667      0.       0        
GRID         123       0 0.51667 0.58333      0.       0        
GRID         124       0 0.43333 0.66667      0.       0        
GRID         125       0 0.43333 0.58333      0.       0        
GRID         126       0      0.      1.      0.       0        
GRID         127       00.085714      1.      0.       0        
GRID         128       0 0.17143      1.      0.       0        
GRID         129       0 0.25714      1.      0.       0        
GRID         130       0 0.34286      1.      0.       0        
GRID         131       0 0.42857      1.      0.       0        
GRID         132       0 0.51429      1.      0.       0        
GRID         133       0     0.6      1.      0.       0        
GRID         134       0     0.6 0.91667      0.       0        
GRID         136       0     0.6    0.75      0.       0        
GRID         137       0    0.35    0.75      0.       0        
GRID         138       0    0.35 0.58333      0.       0        
GRID         139       0  0.2625     0.5      0.       0        
GRID         141       0  0.0875     0.5      0.       0        
GRID         143       0      0. 0.58333      0.       0        
GRID         144       0      0. 0.66667      0.       0        
GRID         145       0      0.    0.75      0.       0        
GRID         146       0      0. 0.83333      0.       0        
GRID         147       0      0. 0.91667      0.       0        
GRID         148       0  0.2617 0.58367      0.       0        
GRID         149       0 0.26103 0.66724      0.       0        
GRID         150       0 0.26054 0.75054      0.       0        
GRID         151       0 0.25916 0.83372      0.       0        
GRID         152       0 0.25811 0.91678      0.       0        
GRID         153       0     0.6 0.33333      0.       0        
GRID         154       0     0.6 0.41667      0.       0        
GRID         155       0     0.6     0.5      0.       0        
GRID         156       0 0.51667     0.5      0.       0        
GRID         159       0 0.34642 0.83313      0.       0        
GRID         160       0    0.35 0.33333      0.       0        
GRID         161       0 0.34443 0.91641      0.       0        
GRID         162       0 0.43333    0.25      0.       0        
GRID         163       0 0.51667    0.25      0.       0        
GRID         164       0 0.51667 0.33333      0.       0        
GRID         165       0 0.51667 0.41667      0.       0        
GRID         166       0 0.43333 0.33333      0.       0        
GRID         167       0 0.43333 0.41667      0.       0        
GRID         168       0      0.      0.      0.       0        
GRID         169       00.085714      0.      0.       0        
GRID         170       0 0.17143      0.      0.       0        
GRID         171       0 0.25714      0.      0.       0        
GRID         172       0 0.34286      0.      0.       0        
GRID         173       0 0.42857      0.      0.       0        
GRID         174       0 0.51429      0.      0.       0        
GRID         175       0     0.6      0.      0.       0        
GRID         177       0     0.6 0.16667      0.       0        
GRID         178       0     0.6    0.25      0.       0        
GRID         179       0 0.51588 0.83297      0.       0        
GRID         180       0 0.51509  0.9163      0.       0        
GRID         181       0    0.35    0.25      0.       0        
GRID         182       0 0.43147 0.83291      0.       0        
GRID         183       0    0.35 0.41667      0.       0        
GRID         184       0 0.43003 0.91621      0.       0        
GRID         186       0   0.175     0.5      0.       0        
GRID         188       0      0.     0.5      0.       0        
GRID         189       0      0. 0.41667      0.       0        
GRID         190       0      0. 0.33333      0.       0        
GRID         191       0      0.    0.25      0.       0        
GRID         192       0      0. 0.16667      0.       0        
GRID         193       0      0.0.083333      0.       0        
GRID         194       0  0.2617 0.41633      0.       0        
GRID         195       0 0.26103 0.33276      0.       0        
GRID         196       0 0.26054 0.24946      0.       0        
GRID         197       0 0.25916 0.16628      0.       0        
GRID         198       0 0.258110.083217      0.       0        
GRID         199       0 0.34642 0.16687      0.       0        
GRID         200       0 0.344430.083587      0.       0        
GRID         201       0 0.51588 0.16703      0.       0        
GRID         202       0 0.51509  0.0837      0.       0        
GRID         203       0 0.43147 0.16709      0.       0        
GRID         204       0 0.430030.083793      0.       0        
GRID         205       0 0.08571  0.2496      0.       0        
GRID         206       0 0.17243 0.24926      0.       0        
GRID         207       00.086007 0.33295      0.       0        
GRID         208       0 0.17297 0.33267      0.       0        
GRID         209       00.086621 0.41643      0.       0        
GRID         210       0 0.17377 0.41627      0.       0        
GRID         211       0 0.17198 0.16616      0.       0        
GRID         212       0 0.171690.083093      0.       0        
GRID         213       00.0855860.083162      0.       0        
GRID         214       00.085623 0.16634      0.       0        
GRID         215       0 0.08571  0.7504      0.       0        
GRID         216       0 0.17243 0.75074      0.       0        
GRID         217       00.086007 0.66705      0.       0        
GRID         218       0 0.17297 0.66733      0.       0        
GRID         219       00.086621 0.58357      0.       0        
GRID         220       0 0.17377 0.58373      0.       0        
GRID         221       0 0.17198 0.83384      0.       0        
GRID         222       0 0.17169 0.91691      0.       0        
GRID         223       00.085586 0.91684      0.       0        
GRID         224       00.085623 0.83366      0.       0        
CQUAD4         1       2     153     178       9      10                
CQUAD4         2       2     154     153      10      11                
CQUAD4         3       2     155     154      11      60                
CQUAD4         4       2      10       9       8      12                
CQUAD4         5       2      11      10      12      13                
CQUAD4         6       2      60      11      13       5                
CQUAD4         7       2      12       8      25       7                
CQUAD4         8       2      13      12       7      26                
CQUAD4         9       2       5      13      26      62                
CQUAD4        10       1     178     177      43       9                
CQUAD4        11       1     177      22      44      43                
CQUAD4        12       1     175      20      44      22                
CQUAD4        13       1       9      43      45       8                
CQUAD4        14       1      43      44      46      45                
CQUAD4        15       1      20      19      46      44                
CQUAD4        16       1       8      45      41      25                
CQUAD4        17       1      45      46      42      41                
CQUAD4        18       1      19      18      42      46                
CQUAD4        19       1      62      26      36      27                
CQUAD4        20       1      26       7      37      36                
CQUAD4        21       1       7      25      38      37                
CQUAD4        22       1      25      41      39      38                
CQUAD4        23       1      41      42      40      39                
CQUAD4        24       1      18      17      40      42                
CQUAD4        25       1      30      29      51      31                
CQUAD4        26       1      29      28      52      51                
CQUAD4        27       1      27      36      52      28                
CQUAD4        28       1      31      51      49      32                
CQUAD4        29       1      51      52      50      49                
CQUAD4        30       1      36      37      50      52                
CQUAD4        31       1      32      49      47      33                
CQUAD4        32       1      49      50      48      47                
CQUAD4        33       1      37      38      48      50                
CQUAD4        34       1      38      39      53      48                
CQUAD4        35       1      39      40      54      53                
CQUAD4        36       1      17      16      54      40                
CQUAD4        37       1      35      34      56      55                
CQUAD4        38       1      33      47      56      34                
CQUAD4        39       1      48      53      56      47                
CQUAD4        40       1      53      54      55      56                
CQUAD4        41       1      16      15      55      54                
CQUAD4        42       1      15      14      35      55                
CQUAD4        43       2     136     113      66      65                
CQUAD4        44       2     113     114      67      66                
CQUAD4        45       2     114     155      60      67                
CQUAD4        46       2      65      66      68      64                
CQUAD4        47       2      66      67      69      68                
CQUAD4        48       2      67      60       5      69                
CQUAD4        49       2      64      68      63      81                
CQUAD4        50       2      68      69      82      63                
CQUAD4        51       2      69       5      62      82                
CQUAD4        52       1      79     136      65      99                
CQUAD4        53       1     134      79      99     100                
CQUAD4        54       1      76     133     134     100                
CQUAD4        55       1      99      65      64     101                
CQUAD4        56       1     100      99     101     102                
CQUAD4        57       1      75      76     100     102                
CQUAD4        58       1     101      64      81      97                
CQUAD4        59       1     102     101      97      98                
CQUAD4        60       1      74      75     102      98                
CQUAD4        61       1      82      62      27      92                
CQUAD4        62       1      63      82      92      93                
CQUAD4        63       1      81      63      93      94                
CQUAD4        64       1      97      81      94      95                
CQUAD4        65       1      98      97      95      96                
CQUAD4        66       1      73      74      98      96                
CQUAD4        67       1      29      30      87     107                
CQUAD4        68       1      28      29     107     108                
CQUAD4        69       1      92      27      28     108                
CQUAD4        70       1     107      87      88     105                
CQUAD4        71       2     178     153     164     163                
CQUAD4        72       2     153     154     165     164                
CQUAD4        73       2     154     155     156     165                
CQUAD4        74       2     163     164     166     162                
CQUAD4        75       2     164     165     167     166                
CQUAD4        76       2     165     156     117     167                
CQUAD4        77       2     162     166     160     181                
CQUAD4        78       2     166     167     183     160                
CQUAD4        79       2     167     117     118     183                
CQUAD4        80       1     177     178     163     201                
CQUAD4        81       1      22     177     201     202                
CQUAD4        82       1     174     175      22     202                
CQUAD4        83       1     201     163     162     203                
CQUAD4        84       1     202     201     203     204                
CQUAD4        85       1     173     174     202     204                
CQUAD4        86       1     203     162     181     199                
CQUAD4        87       1     204     203     199     200                
CQUAD4        88       1     172     173     204     200                
CQUAD4        89       1     183     118     139     194                
CQUAD4        90       1     160     183     194     195                
CQUAD4        91       1     181     160     195     196                
CQUAD4        92       1     199     181     196     197                
CQUAD4        93       1     200     199     197     198                
CQUAD4        94       1     171     172     200     198                
CQUAD4        95       1     141     188     189     209                
CQUAD4        96       1     186     141     209     210                
CQUAD4        97       1     194     139     186     210                
CQUAD4        98       1     209     189     190     207                
CQUAD4        99       1     210     209     207     208                
CQUAD4       100       1     195     194     210     208                
CQUAD4       101       1     207     190     191     205                
CQUAD4       102       1     208     207     205     206                
CQUAD4       103       1     196     195     208     206                
CQUAD4       104       1     197     196     206     211                
CQUAD4       105       1     198     197     211     212                
CQUAD4       106       1     170     171     198     212                
CQUAD4       107       1     192     193     213     214                
CQUAD4       108       1     205     191     192     214                
CQUAD4       109       1     211     206     205     214                
CQUAD4       110       1     212     211     214     213                
CQUAD4       111       1     169     170     212     213                
CQUAD4       112       1     168     169     213     193                
CQUAD4       113       1     108     107     105     106                
CQUAD4       114       1      93      92     108     106                
CQUAD4       115       1     105      88      89     103                
CQUAD4       116       1     106     105     103     104                
CQUAD4       117       1      94      93     106     104                
CQUAD4       118       1      95      94     104     109                
CQUAD4       119       1      96      95     109     110                
CQUAD4       120       1      72      73      96     110                
CQUAD4       121       1      90      91     111     112                
CQUAD4       122       1     103      89      90     112                
CQUAD4       123       1     109     104     103     112                
CQUAD4       124       1     110     109     112     111                
CQUAD4       125       1      71      72     110     111                
CQUAD4       126       1      70      71     111      91                
CQUAD4       127       2     113     136     121     122                
CQUAD4       128       2     114     113     122     123                
CQUAD4       129       2     155     114     123     156                
CQUAD4       130       2     122     121     120     124                
CQUAD4       131       2     123     122     124     125                
CQUAD4       132       2     156     123     125     117                
CQUAD4       133       2     124     120     137     119                
CQUAD4       134       2     125     124     119     138                
CQUAD4       135       2     117     125     138     118                
CQUAD4       136       1     136      79     179     121                
CQUAD4       137       1      79     134     180     179                
CQUAD4       138       1     133     132     180     134                
CQUAD4       139       1     121     179     182     120                
CQUAD4       140       1     179     180     184     182                
CQUAD4       141       1     132     131     184     180                
CQUAD4       142       1     120     182     159     137                
CQUAD4       143       1     182     184     161     159                
CQUAD4       144       1     131     130     161     184                
CQUAD4       145       1     118     138     148     139                
CQUAD4       146       1     138     119     149     148                
CQUAD4       147       1     119     137     150     149                
CQUAD4       148       1     137     159     151     150                
CQUAD4       149       1     159     161     152     151                
CQUAD4       150       1     130     129     152     161                
CQUAD4       151       1     188     141     219     143                
CQUAD4       152       1     141     186     220     219                
CQUAD4       153       1     139     148     220     186                
CQUAD4       154       1     143     219     217     144                
CQUAD4       155       1     219     220     218     217                
CQUAD4       156       1     148     149     218     220                
CQUAD4       157       1     144     217     215     145                
CQUAD4       158       1     217     218     216     215                
CQUAD4       159       1     149     150     216     218                
CQUAD4       160       1     150     151     221     216                
CQUAD4       161       1     151     152     222     221                
CQUAD4       162       1     129     128     222     152                
CQUAD4       163       1     147     146     224     223                
CQUAD4       164       1     145     215     224     146                
CQUAD4       165       1     216     221     224     215                
CQUAD4       166       1     221     222     223     224                
CQUAD4       167       1     128     127     223     222                
CQUAD4       168       1     127     126     147     223                
ENDDATA
