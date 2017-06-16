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
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\Cylinder\Solid\cylinder.MOD
$   Date       : Tue Dec 12 15:32:46 2006
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
$ FEMAP Property 1 : Cylinder
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Inclusion
PSHELL         2       1      0.       1               1              0.
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           3       0   0.125-0.21651     0.0       0        
GRID           4       0    0.25      0.     0.0       0        
GRID           5       0   0.125 0.21651     0.0       0        
GRID           6       0  -0.125 0.21651     0.0       0        
GRID           7       0      0.      0.     0.0       0        
GRID           8       0   -0.55      0.     0.0       0        
GRID           9       0-0.47631  -0.275     0.0       0        
GRID          10       0  -0.275-0.47631     0.0       0        
GRID          11       0      0.   -0.55     0.0       0        
GRID          12       0   0.275-0.47631     0.0       0        
GRID          13       0 0.47631  -0.275     0.0       0        
GRID          14       0    0.55      0.     0.0       0        
GRID          15       0 0.47631   0.275     0.0       0        
GRID          16       0   0.275 0.47631     0.0       0        
GRID          17       0      0.    0.55     0.0       0        
GRID          18       0  -0.275 0.47631     0.0       0        
GRID          19       0-0.47631   0.275     0.0       0        
GRID          22       0  -0.125-0.21651     0.0       0        
GRID          23       0   -0.25      0.     0.0       0        
CTRIA3         1       2       3       4       7                        
CTRIA3         2       2      22       3       7                        
CTRIA3         3       2      23      22       7                        
CTRIA3         4       2       7       4       5                        
CTRIA3         5       2       7       5       6                        
CTRIA3         6       2      23       7       6                        
CTRIA3         7       1      13      14       4                        
CTRIA3         8       1      13       4       3                        
CTRIA3         9       1      12      13       3                        
CTRIA3        10       1      11      12       3                        
CTRIA3        11       1      11       3      22                        
CTRIA3        12       1      10      11      22                        
CTRIA3        13       1       9      10      22                        
CTRIA3        14       1       9      22      23                        
CTRIA3        15       1       8       9      23                        
CTRIA3        16       1       4      14      15                        
CTRIA3        17       1       5       4      15                        
CTRIA3        18       1       5      15      16                        
CTRIA3        19       1       5      16      17                        
CTRIA3        20       1       6       5      17                        
CTRIA3        21       1       6      17      18                        
CTRIA3        22       1       6      18      19                        
CTRIA3        23       1      23       6      19                        
CTRIA3        24       1       8      23      19                        
ENDDATA
