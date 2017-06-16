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
$   Date       : Tue Dec 12 15:50:15 2006
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
GRID           3       0    0.25      0.     0.0       0        
GRID           5       0   -0.55      0.     0.0       0        
GRID           6       0-0.38891-0.38891     0.0       0        
GRID           7       0      0.   -0.55     0.0       0        
GRID           8       0 0.38891-0.38891     0.0       0        
GRID           9       0    0.55      0.     0.0       0        
GRID          10       0 0.38891 0.38891     0.0       0        
GRID          11       0      0.    0.55     0.0       0        
GRID          12       0-0.38891 0.38891     0.0       0        
GRID          14       0      0.   -0.25     0.0       0        
GRID          15       0   -0.25      0.     0.0       0        
GRID          16       0      0.    0.25     0.0       0        
CTRIA3         1       2      15      14       3                        
CTRIA3         2       2      15       3      16                        
CTRIA3         3       1       8       9       3                        
CTRIA3         4       1       8       3      14                        
CTRIA3         5       1       7       8      14                        
CTRIA3         6       1       6       7      14                        
CTRIA3         7       1       6      14      15                        
CTRIA3         8       1       5       6      15                        
CTRIA3         9       1       3       9      10                        
CTRIA3        10       1      16       3      10                        
CTRIA3        11       1      16      10      11                        
CTRIA3        12       1      16      11      12                        
CTRIA3        13       1      15      16      12                        
CTRIA3        14       1       5      15      12                        
ENDDATA
