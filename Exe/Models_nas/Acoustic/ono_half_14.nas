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
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\Acoustic\Solid\ono_box.MOD
$   Date       : Fri Dec 15 07:12:17 2006
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
$ FEMAP Property 1 : Ono box
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Inclusion
PSHELL         2       1      0.       1               1              0.
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID         174       0     0.6   0.375     0.0       0        
GRID         175       0     0.6     0.5     0.0       0        
GRID         176       0   0.475     0.5     0.0       0        
GRID         177       0    0.35     0.5     0.0       0        
GRID         178       0    0.35    0.33     0.0       0        
GRID         179       0    0.43    0.25     0.0       0        
GRID         180       0 0.46375 0.36375     0.0       0        
GRID         182       0   0.175     0.5     0.0       0        
GRID         183       0      0.     0.5     0.0       0        
GRID         184       0      0. 0.33333     0.0       0        
GRID         185       0      0. 0.16667     0.0       0        
GRID         186       0      0.      0.     0.0       0        
GRID         187       0    0.15      0.     0.0       0        
GRID         188       0     0.3      0.     0.0       0        
GRID         189       0    0.45      0.     0.0       0        
GRID         190       0     0.6      0.     0.0       0        
GRID         191       0     0.6   0.125     0.0       0        
GRID         192       0     0.6    0.25     0.0       0        
GRID         195       0 0.44891 0.13282     0.0       0        
GRID         196       0 0.31485 0.15642     0.0       0        
GRID         197       0 0.15934 0.16362     0.0       0        
GRID         198       0 0.17127 0.33163     0.0       0        
CTRIA3       132       2     178     179     180                        
CQUAD4       133       2     176     177     178     180                
CQUAD4       134       2     174     175     176     180                
CQUAD4       135       2     192     174     180     179                
CQUAD4       136       1     189     190     191     195                
CQUAD4       137       1     188     189     195     196                
CQUAD4       138       1     187     188     196     197                
CQUAD4       139       1     185     186     187     197                
CQUAD4       140       1     195     191     192     179                
CQUAD4       141       1     196     195     179     178                
CQUAD4       142       1     197     196     178     198                
CQUAD4       143       1     184     185     197     198                
CQUAD4       144       1     182     183     184     198                
CQUAD4       145       1     177     182     198     178                
ENDDATA
