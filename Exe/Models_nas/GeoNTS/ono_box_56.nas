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
$   Date       : Sat Feb 23 18:35:13 2008
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
GRID          69       0      0.     0.5      0.       0        
GRID          70       0   0.175     0.5      0.       0        
GRID          73       0    0.43    0.75      0.       0        
GRID          76       0     0.6      1.      0.       0        
GRID          77       0    0.45      1.      0.       0        
GRID          78       0     0.3      1.      0.       0        
GRID          79       0    0.15      1.      0.       0        
GRID          80       0      0.      1.      0.       0        
GRID          81       0      0. 0.83333      0.       0        
GRID          82       0      0. 0.66667      0.       0        
GRID          83       0 0.15934 0.83638      0.       0        
GRID          84       0 0.31485 0.84358      0.       0        
GRID          85       0 0.44891 0.86718      0.       0        
GRID          86       0 0.17127 0.66837      0.       0        
GRID          90       0     1.2 0.66667      0.       0        
GRID          91       0     1.2 0.83333      0.       0        
GRID          92       0     1.2      1.      0.       0        
GRID          93       0    1.05      1.      0.       0        
GRID          94       0     0.9      1.      0.       0        
GRID          95       0    0.75      1.      0.       0        
GRID          97       0     0.6   0.875      0.       0        
GRID          98       0     0.6    0.75      0.       0        
GRID         101       0 0.75109 0.86718      0.       0        
GRID         102       0 0.88515 0.84358      0.       0        
GRID         103       0 1.04066 0.83638      0.       0        
GRID         104       0 1.02873 0.66837      0.       0        
GRID         105       0     1.2      0.      0.       0        
GRID         106       0     1.2 0.16667      0.       0        
GRID         107       0     1.2 0.33333      0.       0        
GRID         108       0     1.2     0.5      0.       0        
GRID         109       0   1.025     0.5      0.       0        
GRID         110       0    0.85     0.5      0.       0        
GRID         113       0     0.6    0.25      0.       0        
GRID         115       0     0.6      0.      0.       0        
GRID         116       0    0.75      0.      0.       0        
GRID         117       0     0.9      0.      0.       0        
GRID         118       0    1.05      0.      0.       0        
GRID         119       0 0.75109 0.13282      0.       0        
GRID         120       0 0.88515 0.15642      0.       0        
GRID         121       0 1.04066 0.16362      0.       0        
GRID         122       0 1.02873 0.33163      0.       0        
GRID         126       0      0. 0.33333      0.       0        
GRID         127       0      0. 0.16667      0.       0        
GRID         128       0      0.      0.      0.       0        
GRID         129       0    0.15      0.      0.       0        
GRID         130       0     0.3      0.      0.       0        
GRID         131       0    0.45      0.      0.       0        
GRID         133       0     0.6   0.125      0.       0        
GRID         135       0    0.43    0.25      0.       0        
GRID         137       0 0.44891 0.13282      0.       0        
GRID         138       0 0.31485 0.15642      0.       0        
GRID         139       0 0.15934 0.16362      0.       0        
GRID         140       0 0.17127 0.33163      0.       0        
GRID         141       0    0.35     0.5      0.       0        
GRID         147       0    0.35    0.67      0.       0        
GRID         152       0    0.85    0.67      0.       0        
GRID         153       0    0.77    0.75      0.       0        
GRID         162       0    0.77    0.25      0.       0        
GRID         163       0    0.85    0.33      0.       0        
GRID         170       0    0.35    0.33      0.       0        
CQUAD4        76       1      79      80      81      83                
CQUAD4        77       1      78      79      83      84                
CQUAD4        78       1      77      78      84      85                
CQUAD4        79       1      97      76      77      85                
CQUAD4        80       1      83      81      82      86                
CQUAD4        81       1      73      98      97      85                
CQUAD4        82       1     147      73      85      84                
CQUAD4        83       1     147      84      83      86                
CQUAD4        84       1      70     141     147      86                
CQUAD4        85       1      69      70      86      82                
CQUAD4        86       1      95      76      97     101                
CQUAD4        87       1      94      95     101     102                
CQUAD4        88       1      93      94     102     103                
CQUAD4        89       1      91      92      93     103                
CQUAD4        90       1     101      97      98     153                
CQUAD4        91       1     102     101     153     152                
CQUAD4        92       1     103     102     152     104                
CQUAD4        93       1      90      91     103     104                
CQUAD4        94       1     109     108      90     104                
CQUAD4        95       1     110     109     104     152                
CQUAD4        96       1     109     110     163     122                
CQUAD4        97       1     107     108     109     122                
CQUAD4        98       1     162     113     133     119                
CQUAD4        99       1     163     162     119     120                
CQUAD4       100       1     122     163     120     121                
CQUAD4       101       1     106     107     122     121                
CQUAD4       102       1     119     133     115     116                
CQUAD4       103       1     120     119     116     117                
CQUAD4       104       1     121     120     117     118                
CQUAD4       105       1     105     106     121     118                
CQUAD4       106       1     131     115     133     137                
CQUAD4       107       1     130     131     137     138                
CQUAD4       108       1     129     130     138     139                
CQUAD4       109       1     127     128     129     139                
CQUAD4       110       1     137     133     113     135                
CQUAD4       111       1     138     137     135     170                
CQUAD4       112       1     139     138     170     140                
CQUAD4       113       1     126     127     139     140                
CQUAD4       114       1      70      69     126     140                
CQUAD4       115       1     141      70     140     170                
ENDDATA
