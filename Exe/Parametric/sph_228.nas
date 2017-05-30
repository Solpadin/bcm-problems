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
BEGIN BULK
$ ***************************************************************************
$   Written by : FEMAP
$   Version    : 7.00
$   Translator : MSC/NASTRAN
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\Parametric\Solid\sphere1.MOD
$   Date       : Sat Nov 19 16:48:13 2005
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
$ FEMAP Property 1 : Sphere surface
PSHELL         1       1      0.       1               1              0.
$ FEMAP Material 1 : Mat1
MAT1           1                              0.      0.      0.        
GRID          57       0     -1.      0.      0.       0        
GRID          58       0-0.96593      0. 0.25882       0        
GRID          60       0-0.70711      0. 0.70711       0        
GRID          61       0    -0.5      0. 0.86603       0        
GRID          62       0-0.25882      0. 0.96593       0        
GRID          64       0      0. 0.25882 0.96593       0        
GRID          65       0      0.     0.5 0.86603       0        
GRID          66       0      0. 0.70711 0.70711       0        
GRID          67       0      0. 0.86603     0.5       0        
GRID          68       0      0. 0.96593 0.25882       0        
GRID          70       0-0.25882 0.96593      0.       0        
GRID          74       0-0.96593 0.25882      0.       0        
GRID          75       0-0.51617 0.78226 0.34877       0        
GRID          76       0-0.51705 0.60525 0.60525       0        
GRID          77       0-0.51617 0.34877 0.78226       0        
GRID          78       0-0.35688 0.66054 0.66054       0        
GRID          79       0-0.18394 0.69504 0.69504       0        
GRID          80       0-0.35631 0.82346 0.44155       0        
GRID          81       0 -0.1874 0.85516 0.48331       0        
GRID          82       0-0.20195 0.94737 0.24842       0        
GRID          83       0-0.35607 0.90551 0.23078       0        
GRID          84       0-0.35631 0.44155 0.82346       0        
GRID          85       0 -0.1874 0.48331 0.85516       0        
GRID          86       0-0.20195 0.24842 0.94737       0        
GRID          87       0-0.35607 0.23078 0.90551       0        
GRID          88       0-0.67141 0.68259 0.28858       0        
GRID          89       0-0.66467  0.5283  0.5283       0        
GRID          90       0-0.67141 0.28858 0.68259       0        
GRID          91       0-0.86847 0.35054 0.35054       0        
GRID          92       0-0.77713 0.44501 0.44501       0        
GRID          93       0-0.77802 0.57984 0.24183       0        
GRID          94       0-0.77802 0.24183 0.57984       0        
GRID          95       0-0.94023  0.2408  0.2408       0        
GRID          99       0-0.70711      0.-0.70711       0        
GRID         101       0-0.25882      0.-0.96593       0        
GRID         103       0      0.-0.25882-0.96593       0        
GRID         104       0      0.    -0.5-0.86603       0        
GRID         105       0      0.-0.70711-0.70711       0        
GRID         106       0      0.-0.86603    -0.5       0        
GRID         107       0      0.-0.96593-0.25882       0        
GRID         109       0-0.25882-0.96593      0.       0        
GRID         111       0-0.70711-0.70711      0.       0        
GRID         112       0-0.86603    -0.5      0.       0        
GRID         113       0-0.96593-0.25882      0.       0        
GRID         114       0-0.31123-0.86764-0.38772       0        
GRID         115       0-0.32272-0.66927-0.66927       0        
GRID         116       0-0.31123-0.38772-0.86764       0        
GRID         117       0-0.15795-0.69823-0.69823       0        
GRID         118       0-0.15268-0.86945-0.46983       0        
GRID         119       0-0.13785-0.95907-0.24736       0        
GRID         120       0-0.15268-0.46983-0.86945       0        
GRID         121       0-0.13785-0.24736-0.95907       0        
GRID         122       0-0.66544-0.68905-0.28706       0        
GRID         123       0-0.65689-0.53315-0.53315       0        
GRID         124       0-0.66544-0.28706-0.68905       0        
GRID         125       0-0.49916-0.61271-0.61271       0        
GRID         126       0-0.49942-0.79793-0.33746       0        
GRID         127       0-0.49942-0.33746-0.79793       0        
GRID         128       0-0.86793-0.35122-0.35122       0        
GRID         129       0-0.77434-0.44744-0.44744       0        
GRID         130       0-0.77546-0.58317-0.24203       0        
GRID         131       0-0.77546-0.24203-0.58317       0        
GRID         132       0-0.94013-0.24099-0.24099       0        
GRID         137       0    -0.5-0.86603      0.       0        
GRID         139       0      0.     -1.      0.       0        
GRID         140       0      0.-0.96593 0.25882       0        
GRID         141       0      0.-0.86603     0.5       0        
GRID         142       0      0.-0.70711 0.70711       0        
GRID         143       0      0.    -0.5 0.86603       0        
GRID         144       0      0.-0.25882 0.96593       0        
GRID         145       0      0.      0.      1.       0        
GRID         149       0-0.86603      0.     0.5       0        
GRID         151       0-0.31123-0.38772 0.86764       0        
GRID         152       0-0.32272-0.66927 0.66927       0        
GRID         153       0-0.31123-0.86764 0.38772       0        
GRID         154       0-0.15795-0.69823 0.69823       0        
GRID         155       0-0.15268-0.46983 0.86945       0        
GRID         156       0-0.13785-0.24736 0.95907       0        
GRID         157       0-0.15268-0.86945 0.46983       0        
GRID         158       0-0.13785-0.95907 0.24736       0        
GRID         159       0-0.49942-0.33746 0.79793       0        
GRID         160       0-0.49916-0.61271 0.61271       0        
GRID         161       0-0.49942-0.79793 0.33746       0        
GRID         162       0-0.66544-0.28706 0.68905       0        
GRID         163       0-0.65689-0.53315 0.53315       0        
GRID         164       0-0.66544-0.68905 0.28706       0        
GRID         165       0-0.86793-0.35122 0.35122       0        
GRID         166       0-0.77434-0.44744 0.44744       0        
GRID         167       0-0.77546-0.24203 0.58317       0        
GRID         168       0-0.77546-0.58317 0.24203       0        
GRID         169       0-0.94013-0.24099 0.24099       0        
GRID         172       0-0.86603     0.5      0.       0        
GRID         173       0-0.70711 0.70711      0.       0        
GRID         174       0    -0.5 0.86603      0.       0        
GRID         176       0      0.      1.      0.       0        
GRID         177       0      0. 0.96593-0.25882       0        
GRID         178       0      0. 0.86603    -0.5       0        
GRID         179       0      0. 0.70711-0.70711       0        
GRID         180       0      0.     0.5-0.86603       0        
GRID         181       0      0. 0.25882-0.96593       0        
GRID         182       0      0.      0.     -1.       0        
GRID         184       0    -0.5      0.-0.86603       0        
GRID         186       0-0.86603      0.    -0.5       0        
GRID         187       0-0.96593      0.-0.25882       0        
GRID         188       0-0.51617 0.34877-0.78226       0        
GRID         189       0-0.51705 0.60525-0.60525       0        
GRID         190       0-0.51617 0.78226-0.34877       0        
GRID         191       0-0.35688 0.66054-0.66054       0        
GRID         192       0-0.18394 0.69504-0.69504       0        
GRID         193       0-0.35631 0.44155-0.82346       0        
GRID         194       0 -0.1874 0.48331-0.85516       0        
GRID         195       0-0.20195 0.24842-0.94737       0        
GRID         196       0-0.35607 0.23078-0.90551       0        
GRID         197       0-0.35631 0.82346-0.44155       0        
GRID         198       0 -0.1874 0.85516-0.48331       0        
GRID         199       0-0.20195 0.94737-0.24842       0        
GRID         200       0-0.35607 0.90551-0.23078       0        
GRID         201       0-0.67141 0.28858-0.68259       0        
GRID         202       0-0.66467  0.5283 -0.5283       0        
GRID         203       0-0.67141 0.68259-0.28858       0        
GRID         204       0-0.86847 0.35054-0.35054       0        
GRID         205       0-0.77713 0.44501-0.44501       0        
GRID         206       0-0.77802 0.24183-0.57984       0        
GRID         207       0-0.77802 0.57984-0.24183       0        
GRID         208       0-0.94023  0.2408 -0.2408       0        
GRID         222       0 0.96593      0.-0.25882       0        
GRID         225       0     0.5      0.-0.86603       0        
GRID         227       0 0.51617 0.34877-0.78226       0        
GRID         228       0 0.51705 0.60525-0.60525       0        
GRID         229       0 0.51617 0.78226-0.34877       0        
GRID         230       0 0.67141 0.28858-0.68259       0        
GRID         231       0 0.66467  0.5283 -0.5283       0        
GRID         232       0 0.67141 0.68259-0.28858       0        
GRID         233       0 0.86847 0.35054-0.35054       0        
GRID         234       0 0.94023  0.2408 -0.2408       0        
GRID         235       0 0.77713 0.44501-0.44501       0        
GRID         236       0 0.77802 0.24183-0.57984       0        
GRID         237       0 0.77802 0.57984-0.24183       0        
GRID         238       0 0.35688 0.66054-0.66054       0        
GRID         239       0 0.18394 0.69504-0.69504       0        
GRID         240       0 0.35631 0.82346-0.44155       0        
GRID         241       0  0.1874 0.85516-0.48331       0        
GRID         242       0 0.20195 0.94737-0.24842       0        
GRID         243       0 0.35607 0.90551-0.23078       0        
GRID         244       0 0.35631 0.44155-0.82346       0        
GRID         245       0  0.1874 0.48331-0.85516       0        
GRID         246       0 0.20195 0.24842-0.94737       0        
GRID         247       0 0.35607 0.23078-0.90551       0        
GRID         256       0     0.5-0.86603      0.       0        
GRID         257       0 0.70711-0.70711      0.       0        
GRID         259       0 0.96593-0.25882      0.       0        
GRID         260       0      1.      0.      0.       0        
GRID         261       0 0.96593      0. 0.25882       0        
GRID         262       0 0.86603      0.     0.5       0        
GRID         264       0     0.5      0. 0.86603       0        
GRID         266       0 0.51617-0.34877 0.78226       0        
GRID         267       0 0.51705-0.60525 0.60525       0        
GRID         268       0 0.51617-0.78226 0.34877       0        
GRID         269       0 0.67141-0.28858 0.68259       0        
GRID         270       0 0.66467 -0.5283  0.5283       0        
GRID         271       0 0.67141-0.68259 0.28858       0        
GRID         272       0 0.86847-0.35054 0.35054       0        
GRID         273       0 0.94023 -0.2408  0.2408       0        
GRID         274       0 0.77713-0.44501 0.44501       0        
GRID         275       0 0.77802-0.24183 0.57984       0        
GRID         276       0 0.77802-0.57984 0.24183       0        
GRID         277       0 0.35688-0.66054 0.66054       0        
GRID         278       0 0.18394-0.69504 0.69504       0        
GRID         279       0 0.35631-0.82346 0.44155       0        
GRID         280       0  0.1874-0.85516 0.48331       0        
GRID         281       0 0.20195-0.94737 0.24842       0        
GRID         282       0 0.35607-0.90551 0.23078       0        
GRID         283       0 0.35631-0.44155 0.82346       0        
GRID         284       0  0.1874-0.48331 0.85516       0        
GRID         285       0 0.20195-0.24842 0.94737       0        
GRID         286       0 0.35607-0.23078 0.90551       0        
GRID         294       0 0.25882      0.-0.96593       0        
GRID         296       0 0.70711      0.-0.70711       0        
GRID         297       0 0.86603      0.    -0.5       0        
GRID         301       0 0.86603    -0.5      0.       0        
GRID         304       0 0.25882-0.96593      0.       0        
GRID         305       0 0.51617-0.78226-0.34877       0        
GRID         306       0 0.51705-0.60525-0.60525       0        
GRID         307       0 0.51617-0.34877-0.78226       0        
GRID         308       0 0.67141-0.68259-0.28858       0        
GRID         309       0 0.66467 -0.5283 -0.5283       0        
GRID         310       0 0.67141-0.28858-0.68259       0        
GRID         311       0 0.86847-0.35054-0.35054       0        
GRID         312       0 0.94023 -0.2408 -0.2408       0        
GRID         313       0 0.77713-0.44501-0.44501       0        
GRID         314       0 0.77802-0.57984-0.24183       0        
GRID         315       0 0.77802-0.24183-0.57984       0        
GRID         316       0 0.35688-0.66054-0.66054       0        
GRID         317       0 0.18394-0.69504-0.69504       0        
GRID         318       0 0.35631-0.44155-0.82346       0        
GRID         319       0  0.1874-0.48331-0.85516       0        
GRID         320       0 0.20195-0.24842-0.94737       0        
GRID         321       0 0.35607-0.23078-0.90551       0        
GRID         322       0 0.35631-0.82346-0.44155       0        
GRID         323       0  0.1874-0.85516-0.48331       0        
GRID         324       0 0.20195-0.94737-0.24842       0        
GRID         325       0 0.35607-0.90551-0.23078       0        
GRID         333       0 0.25882      0. 0.96593       0        
GRID         335       0 0.70711      0. 0.70711       0        
GRID         339       0 0.96593 0.25882      0.       0        
GRID         340       0 0.86603     0.5      0.       0        
GRID         341       0 0.70711 0.70711      0.       0        
GRID         342       0     0.5 0.86603      0.       0        
GRID         343       0 0.25882 0.96593      0.       0        
GRID         344       0 0.51617 0.78226 0.34877       0        
GRID         345       0 0.51705 0.60525 0.60525       0        
GRID         346       0 0.51617 0.34877 0.78226       0        
GRID         347       0 0.67141 0.68259 0.28858       0        
GRID         348       0 0.66467  0.5283  0.5283       0        
GRID         349       0 0.67141 0.28858 0.68259       0        
GRID         350       0 0.86847 0.35054 0.35054       0        
GRID         351       0 0.94023  0.2408  0.2408       0        
GRID         352       0 0.77713 0.44501 0.44501       0        
GRID         353       0 0.77802 0.57984 0.24183       0        
GRID         354       0 0.77802 0.24183 0.57984       0        
GRID         355       0 0.35688 0.66054 0.66054       0        
GRID         356       0 0.18394 0.69504 0.69504       0        
GRID         357       0 0.35631 0.44155 0.82346       0        
GRID         358       0  0.1874 0.48331 0.85516       0        
GRID         359       0 0.20195 0.24842 0.94737       0        
GRID         360       0 0.35607 0.23078 0.90551       0        
GRID         361       0 0.35631 0.82346 0.44155       0        
GRID         362       0  0.1874 0.85516 0.48331       0        
GRID         363       0 0.20195 0.94737 0.24842       0        
GRID         364       0 0.35607 0.90551 0.23078       0        
CQUAD4        25       1      80      81      82      83                
CQUAD4        26       1     174      75      80      83                
CQUAD4        27       1      70     174      83      82                
CQUAD4        28       1      68     176      70      82                
CQUAD4        29       1      67      68      82      81                
CQUAD4        30       1      80      75      76      78                
CQUAD4        31       1      81      80      78      79                
CQUAD4        32       1      66      67      81      79                
CQUAD4        33       1      78      76      77      84                
CQUAD4        34       1      79      78      84      85                
CQUAD4        35       1      65      66      79      85                
CQUAD4        36       1      64      65      85      86                
CQUAD4        37       1      62     145      64      86                
CQUAD4        38       1      86      85      84      87                
CQUAD4        39       1      61      62      86      87                
CQUAD4        40       1      61      87      84      77                
CQUAD4        41       1      75     174     173      88                
CQUAD4        42       1      76      75      88      89                
CQUAD4        43       1      77      76      89      90                
CQUAD4        44       1      60      61      77      90                
CQUAD4        45       1      88     173     172      93                
CQUAD4        46       1      93     172      91      92                
CQUAD4        47       1      89      88      93      92                
CQUAD4        48       1     149      60      90      94                
CQUAD4        49       1      94      90      89      92                
CQUAD4        50       1     149      94      92      91                
CQUAD4        51       1      91     172      74      95                
CQUAD4        52       1      58     149      91      95                
CQUAD4        53       1      57      58      95      74                
CQUAD4        54       1     107     139     109     119                
CQUAD4        55       1     119     109     114     118                
CQUAD4        56       1     106     107     119     118                
CQUAD4        57       1     118     114     115     117                
CQUAD4        58       1     105     106     118     117                
CQUAD4        59       1     117     115     116     120                
CQUAD4        60       1     104     105     117     120                
CQUAD4        61       1     101     182     103     121                
CQUAD4        62       1     121     103     104     120                
CQUAD4        63       1     101     121     120     116                
CQUAD4        64       1     137     111     122     126                
CQUAD4        65       1     114     109     137     126                
CQUAD4        66       1     126     122     123     125                
CQUAD4        67       1     115     114     126     125                
CQUAD4        68       1     125     123     124     127                
CQUAD4        69       1     116     115     125     127                
CQUAD4        70       1     184     101     116     127                
CQUAD4        71       1      99     184     127     124                
CQUAD4        72       1     122     111     112     130                
CQUAD4        73       1     130     112     128     129                
CQUAD4        74       1     123     122     130     129                
CQUAD4        75       1     186      99     124     131                
CQUAD4        76       1     131     124     123     129                
CQUAD4        77       1     186     131     129     128                
CQUAD4        78       1     128     112     113     132                
CQUAD4        79       1     187     186     128     132                
CQUAD4        80       1      57     187     132     113                
CQUAD4        81       1     144     145      62     156                
CQUAD4        82       1     156      62     151     155                
CQUAD4        83       1     143     144     156     155                
CQUAD4        84       1     155     151     152     154                
CQUAD4        85       1     142     143     155     154                
CQUAD4        86       1     154     152     153     157                
CQUAD4        87       1     141     142     154     157                
CQUAD4        88       1     109     139     140     158                
CQUAD4        89       1     158     140     141     157                
CQUAD4        90       1     109     158     157     153                
CQUAD4        91       1     151      62      61     159                
CQUAD4        92       1     152     151     159     160                
CQUAD4        93       1     153     152     160     161                
CQUAD4        94       1     137     109     153     161                
CQUAD4        95       1     159      61      60     162                
CQUAD4        96       1     160     159     162     163                
CQUAD4        97       1     161     160     163     164                
CQUAD4        98       1     111     137     161     164                
CQUAD4        99       1     162      60     149     167                
CQUAD4       100       1     167     149     165     166                
CQUAD4       101       1     163     162     167     166                
CQUAD4       102       1     112     111     164     168                
CQUAD4       103       1     168     164     163     166                
CQUAD4       104       1     112     168     166     165                
CQUAD4       105       1     165     149      58     169                
CQUAD4       106       1     113     112     165     169                
CQUAD4       107       1      57     113     169      58                
CQUAD4       108       1     193     194     195     196                
CQUAD4       109       1     184     188     193     196                
CQUAD4       110       1     101     184     196     195                
CQUAD4       111       1     181     182     101     195                
CQUAD4       112       1     180     181     195     194                
CQUAD4       113       1     193     188     189     191                
CQUAD4       114       1     194     193     191     192                
CQUAD4       115       1     179     180     194     192                
CQUAD4       116       1     191     189     190     197                
CQUAD4       117       1     192     191     197     198                
CQUAD4       118       1     178     179     192     198                
CQUAD4       119       1     177     178     198     199                
CQUAD4       120       1      70     176     177     199                
CQUAD4       121       1     199     198     197     200                
CQUAD4       122       1     174      70     199     200                
CQUAD4       123       1     174     200     197     190                
CQUAD4       124       1     188     184      99     201                
CQUAD4       125       1     189     188     201     202                
CQUAD4       126       1     190     189     202     203                
CQUAD4       127       1     173     174     190     203                
CQUAD4       128       1     201      99     186     206                
CQUAD4       129       1     206     186     204     205                
CQUAD4       130       1     202     201     206     205                
CQUAD4       131       1     172     173     203     207                
CQUAD4       132       1     207     203     202     205                
CQUAD4       133       1     172     207     205     204                
CQUAD4       134       1     204     186     187     208                
CQUAD4       135       1      74     172     204     208                
CQUAD4       136       1      57      74     208     187                
CQUAD4       137       1     339     260     222     234                
CQUAD4       138       1     234     222     297     233                
CQUAD4       139       1     340     339     234     233                
CQUAD4       140       1     297     296     230     236                
CQUAD4       141       1     236     230     231     235                
CQUAD4       142       1     233     297     236     235                
CQUAD4       143       1     235     231     232     237                
CQUAD4       144       1     340     233     235     237                
CQUAD4       145       1     341     340     237     232                
CQUAD4       146       1     230     296     225     227                
CQUAD4       147       1     231     230     227     228                
CQUAD4       148       1     232     231     228     229                
CQUAD4       149       1     342     341     232     229                
CQUAD4       150       1     240     241     242     243                
CQUAD4       151       1     342     229     240     243                
CQUAD4       152       1     343     342     243     242                
CQUAD4       153       1     177     176     343     242                
CQUAD4       154       1     178     177     242     241                
CQUAD4       155       1     240     229     228     238                
CQUAD4       156       1     241     240     238     239                
CQUAD4       157       1     179     178     241     239                
CQUAD4       158       1     238     228     227     244                
CQUAD4       159       1     239     238     244     245                
CQUAD4       160       1     180     179     239     245                
CQUAD4       161       1     225     294     246     247                
CQUAD4       162       1     244     227     225     247                
CQUAD4       163       1     245     244     247     246                
CQUAD4       164       1     181     180     245     246                
CQUAD4       165       1     182     181     246     294                
CQUAD4       166       1     259     260     261     273                
CQUAD4       167       1     273     261     262     272                
CQUAD4       168       1     301     259     273     272                
CQUAD4       169       1     262     335     269     275                
CQUAD4       170       1     275     269     270     274                
CQUAD4       171       1     272     262     275     274                
CQUAD4       172       1     274     270     271     276                
CQUAD4       173       1     301     272     274     276                
CQUAD4       174       1     257     301     276     271                
CQUAD4       175       1     269     335     264     266                
CQUAD4       176       1     270     269     266     267                
CQUAD4       177       1     271     270     267     268                
CQUAD4       178       1     256     257     271     268                
CQUAD4       179       1     279     280     281     282                
CQUAD4       180       1     256     268     279     282                
CQUAD4       181       1     304     256     282     281                
CQUAD4       182       1     140     139     304     281                
CQUAD4       183       1     141     140     281     280                
CQUAD4       184       1     279     268     267     277                
CQUAD4       185       1     280     279     277     278                
CQUAD4       186       1     142     141     280     278                
CQUAD4       187       1     277     267     266     283                
CQUAD4       188       1     278     277     283     284                
CQUAD4       189       1     143     142     278     284                
CQUAD4       190       1     264     333     285     286                
CQUAD4       191       1     283     266     264     286                
CQUAD4       192       1     284     283     286     285                
CQUAD4       193       1     144     143     284     285                
CQUAD4       194       1     145     144     285     333                
CQUAD4       195       1     222     260     259     312                
CQUAD4       196       1     312     259     301     311                
CQUAD4       197       1     297     222     312     311                
CQUAD4       198       1     301     257     308     314                
CQUAD4       199       1     314     308     309     313                
CQUAD4       200       1     311     301     314     313                
CQUAD4       201       1     313     309     310     315                
CQUAD4       202       1     297     311     313     315                
CQUAD4       203       1     296     297     315     310                
CQUAD4       204       1     308     257     256     305                
CQUAD4       205       1     309     308     305     306                
CQUAD4       206       1     310     309     306     307                
CQUAD4       207       1     225     296     310     307                
CQUAD4       208       1     318     319     320     321                
CQUAD4       209       1     225     307     318     321                
CQUAD4       210       1     294     225     321     320                
CQUAD4       211       1     103     182     294     320                
CQUAD4       212       1     104     103     320     319                
CQUAD4       213       1     318     307     306     316                
CQUAD4       214       1     319     318     316     317                
CQUAD4       215       1     105     104     319     317                
CQUAD4       216       1     316     306     305     322                
CQUAD4       217       1     317     316     322     323                
CQUAD4       218       1     106     105     317     323                
CQUAD4       219       1     256     304     324     325                
CQUAD4       220       1     322     305     256     325                
CQUAD4       221       1     323     322     325     324                
CQUAD4       222       1     107     106     323     324                
CQUAD4       223       1     139     107     324     304                
CQUAD4       224       1     261     260     339     351                
CQUAD4       225       1     351     339     340     350                
CQUAD4       226       1     262     261     351     350                
CQUAD4       227       1     340     341     347     353                
CQUAD4       228       1     353     347     348     352                
CQUAD4       229       1     350     340     353     352                
CQUAD4       230       1     352     348     349     354                
CQUAD4       231       1     262     350     352     354                
CQUAD4       232       1     335     262     354     349                
CQUAD4       233       1     347     341     342     344                
CQUAD4       234       1     348     347     344     345                
CQUAD4       235       1     349     348     345     346                
CQUAD4       236       1     264     335     349     346                
CQUAD4       237       1     357     358     359     360                
CQUAD4       238       1     264     346     357     360                
CQUAD4       239       1     333     264     360     359                
CQUAD4       240       1      64     145     333     359                
CQUAD4       241       1      65      64     359     358                
CQUAD4       242       1     357     346     345     355                
CQUAD4       243       1     358     357     355     356                
CQUAD4       244       1      66      65     358     356                
CQUAD4       245       1     355     345     344     361                
CQUAD4       246       1     356     355     361     362                
CQUAD4       247       1      67      66     356     362                
CQUAD4       248       1     342     343     363     364                
CQUAD4       249       1     361     344     342     364                
CQUAD4       250       1     362     361     364     363                
CQUAD4       251       1      68      67     362     363                
CQUAD4       252       1     176      68     363     343                
ENDDATA
