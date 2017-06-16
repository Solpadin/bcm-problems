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
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\Acoustic\Solid\ono_box0.MOD
$   Date       : Wed Mar 14 19:17:18 2007
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
$ FEMAP Property 1 : ONO BOX
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : INCLUSION
PSHELL         2       1      0.       1               1              0.
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID         118       0     0.6   0.375     0.0       0        
GRID         119       0     0.6     0.5     0.0       0        
GRID         120       0   0.475     0.5     0.0       0        
GRID         124       0   0.475    0.25     0.0       0        
GRID         125       0   0.475   0.375     0.0       0        
GRID         126       0      0.      0.     0.0       0        
GRID         127       0    0.12      0.     0.0       0        
GRID         128       0    0.24      0.     0.0       0        
GRID         129       0    0.36      0.     0.0       0        
GRID         130       0    0.48      0.     0.0       0        
GRID         131       0     0.6      0.     0.0       0        
GRID         132       0     0.6   0.125     0.0       0        
GRID         133       0     0.6    0.25     0.0       0        
GRID         135       0    0.35    0.25     0.0       0        
GRID         136       0    0.35   0.375     0.0       0        
GRID         137       0    0.35     0.5     0.0       0        
GRID         138       0 0.23333     0.5     0.0       0        
GRID         139       0 0.11667     0.5     0.0       0        
GRID         140       0      0.     0.5     0.0       0        
GRID         141       0      0.   0.375     0.0       0        
GRID         142       0      0.    0.25     0.0       0        
GRID         143       0      0.   0.125     0.0       0        
GRID         144       0 0.23406 0.37562     0.0       0        
GRID         145       0  0.2351 0.25072     0.0       0        
GRID         146       0 0.23776 0.12557     0.0       0        
GRID         147       0 0.35631 0.12494     0.0       0        
GRID         148       0 0.47787 0.12495     0.0       0        
GRID         149       0 0.11734 0.37562     0.0       0        
GRID         150       0 0.11801 0.25088     0.0       0        
GRID         151       0 0.11901 0.12546     0.0       0        
CQUAD4        51       2     133     118     125     124                
CQUAD4        52       2     118     119     120     125                
CQUAD4        53       2     124     125     136     135                
CQUAD4        54       2     125     120     137     136                
CQUAD4        55       1     132     133     124     148                
CQUAD4        56       1     130     131     132     148                
CQUAD4        57       1     148     124     135     147                
CQUAD4        58       1     129     130     148     147                
CQUAD4        59       1     136     137     138     144                
CQUAD4        60       1     135     136     144     145                
CQUAD4        61       1     147     135     145     146                
CQUAD4        62       1     128     129     147     146                
CQUAD4        63       1     144     138     139     149                
CQUAD4        64       1     145     144     149     150                
CQUAD4        65       1     146     145     150     151                
CQUAD4        66       1     127     128     146     151                
CQUAD4        67       1     149     139     140     141                
CQUAD4        68       1     150     149     141     142                
CQUAD4        69       1     151     150     142     143                
CQUAD4        70       1     126     127     151     143                
ENDDATA
