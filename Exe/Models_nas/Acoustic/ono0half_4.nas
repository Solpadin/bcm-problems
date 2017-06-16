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
$   Date       : Sun Mar 11 19:14:23 2007
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
GRID          32       0      0.     0.5     0.0       0        
GRID          33       0      0.      0.     0.0       0        
GRID          36       0      0.    0.25     0.0       0        
GRID          37       0    0.35    0.25     0.0       0        
GRID          38       0    0.35      0.     0.0       0        
GRID          39       0     0.6      0.     0.0       0        
GRID          41       0     0.6    0.25     0.0       0        
GRID          42       0     0.6     0.5     0.0       0        
GRID          43       0    0.35     0.5     0.0       0        
CQUAD4        16       1      36      37      43      32                
CQUAD4        17       1      33      38      37      36                
CQUAD4        18       1      37      38      39      41                
CQUAD4        19       2      41      42      43      37                
ENDDATA
