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
$   Date       : Mon Mar 12 22:39:02 2007
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
GRID          49       0   0.175     0.5     0.0       0        
GRID          50       0      0.     0.5     0.0       0        
GRID          51       0      0.      0.     0.0       0        
GRID          52       0   0.175      0.     0.0       0        
GRID          53       0    0.35      0.     0.0       0        
GRID          54       0    0.35    0.25     0.0       0        
GRID          55       0   0.175    0.25     0.0       0        
GRID          56       0      0.    0.25     0.0       0        
GRID          59       0     0.6      0.     0.0       0        
GRID          60       0     0.6    0.25     0.0       0        
GRID          62       0     0.6     0.5     0.0       0        
GRID          63       0    0.35     0.5     0.0       0        
CQUAD4        20       1      56      55      49      50                
CQUAD4        21       1      55      54      63      49                
CQUAD4        22       1      51      52      55      56                
CQUAD4        23       1      52      53      54      55                
CQUAD4        24       1      54      53      59      60                
CQUAD4        25       2      60      62      63      54                
ENDDATA
