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
$   Date       : Mon Mar 12 22:41:04 2007
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
GRID          66       0     0.6   0.375     0.0       0        
GRID          67       0     0.6     0.5     0.0       0        
GRID          68       0   0.475     0.5     0.0       0        
GRID          69       0    0.35     0.5     0.0       0        
GRID          70       0    0.35   0.375     0.0       0        
GRID          73       0   0.475   0.375     0.0       0        
GRID          75       0   0.175    0.25     0.0       0        
GRID          76       0    0.35    0.25     0.0       0        
GRID          79       0   0.175     0.5     0.0       0        
GRID          80       0      0.     0.5     0.0       0        
GRID          81       0      0.   0.375     0.0       0        
GRID          82       0   0.175   0.375     0.0       0        
GRID          83       0      0.      0.     0.0       0        
GRID          84       0   0.175      0.     0.0       0        
GRID          85       0    0.35      0.     0.0       0        
GRID          89       0      0.    0.25     0.0       0        
GRID          90       0      0.   0.125     0.0       0        
GRID          91       0   0.175   0.125     0.0       0        
GRID          93       0    0.35   0.125     0.0       0        
GRID          95       0   0.475      0.     0.0       0        
GRID          96       0     0.6      0.     0.0       0        
GRID          97       0     0.6   0.125     0.0       0        
GRID          98       0     0.6    0.25     0.0       0        
GRID          99       0   0.475    0.25     0.0       0        
GRID         100       0   0.475   0.125     0.0       0        
CQUAD4        26       2      98      66      73      99                
CQUAD4        27       2      66      67      68      73                
CQUAD4        28       2      99      73      70      76                
CQUAD4        29       2      73      68      69      70                
CQUAD4        30       1      89      75      82      81                
CQUAD4        31       1      75      76      70      82                
CQUAD4        32       1      81      82      79      80                
CQUAD4        33       1      82      70      69      79                
CQUAD4        34       1      83      84      91      90                
CQUAD4        35       1      84      85      93      91                
CQUAD4        36       1      90      91      75      89                
CQUAD4        37       1      91      93      76      75                
CQUAD4        38       1      76      93     100      99                
CQUAD4        39       1      93      85      95     100                
CQUAD4        40       1      99     100      97      98                
CQUAD4        41       1     100      95      96      97                
ENDDATA
