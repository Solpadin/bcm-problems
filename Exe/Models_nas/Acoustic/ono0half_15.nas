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
$   Date       : Sun Mar 11 16:45:47 2007
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
GRID           4       0   0.175     0.5     0.0       0        
GRID           5       0      0.     0.5     0.0       0        
GRID           6       0      0. 0.33333     0.0       0        
GRID           7       0      0. 0.16667     0.0       0        
GRID           8       0      0.      0.     0.0       0        
GRID           9       0    0.15      0.     0.0       0        
GRID          10       0     0.3      0.     0.0       0        
GRID          11       0    0.45      0.     0.0       0        
GRID          12       0     0.6      0.     0.0       0        
GRID          13       0     0.6   0.125     0.0       0        
GRID          15       0   0.475    0.25     0.0       0        
GRID          16       0 0.31994 0.13321     0.0       0        
GRID          17       0 0.16964 0.15573     0.0       0        
GRID          18       0 0.20878 0.32281     0.0       0        
GRID          19       0 0.46121  0.1271     0.0       0        
GRID          20       0     0.6    0.25     0.0       0        
GRID          21       0     0.6   0.375     0.0       0        
GRID          22       0     0.6     0.5     0.0       0        
GRID          23       0   0.475     0.5     0.0       0        
GRID          24       0    0.35     0.5     0.0       0        
GRID          25       0    0.35   0.375     0.0       0        
GRID          26       0    0.35    0.25     0.0       0        
GRID          28       0   0.475   0.375     0.0       0        
CQUAD4         1       1       7       8       9      17                
CQUAD4         2       1       6       7      17      18                
CQUAD4         3       1       4       5       6      18                
CQUAD4         4       1      17       9      10      16                
CQUAD4         5       1      25      24       4      18                
CTRIA3         6       1      26      25      18                        
CQUAD4         7       1      26      18      17      16                
CQUAD4         8       1      13      20      15      19                
CQUAD4         9       1      11      12      13      19                
CQUAD4        10       1      16      10      11      19                
CQUAD4        11       1      26      16      19      15                
CQUAD4        12       2      20      21      28      15                
CQUAD4        13       2      21      22      23      28                
CQUAD4        14       2      15      28      25      26                
CQUAD4        15       2      28      23      24      25                
ENDDATA
