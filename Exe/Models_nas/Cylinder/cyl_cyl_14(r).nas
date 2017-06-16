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
$   Date       : Fri Nov 24 10:10:01 2006
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
GRID           1       0    0.55      0.      0.       0        
GRID           2       0 0.38891 0.38891      0.       0        
GRID           3       0      0.    0.55      0.       0        
GRID           4       0-0.38891 0.38891      0.       0        
GRID           5       0   -0.55      0.      0.       0        
GRID           6       0-0.38891-0.38891      0.       0        
GRID           7       0      0.   -0.55      0.       0        
GRID           8       0 0.38891-0.38891      0.       0        
GRID           9       0   -0.25      0.      0.       0        
GRID          10       0      0.    0.25      0.       0        
GRID          11       0    0.25      0.      0.       0        
GRID          12       0      0.   -0.25      0.       0        
CTRIA3         1       1       4       5       9                        
CTRIA3         2       1       4       9      10                        
CTRIA3         3       1       3       4      10                        
CTRIA3         4       1       2       3      10                        
CTRIA3         5       1       2      10      11                        
CTRIA3         6       1       1       2      11                        
CTRIA3         7       1       9       5       6                        
CTRIA3         8       1      12       9       6                        
CTRIA3         9       1      12       6       7                        
CTRIA3        10       1      12       7       8                        
CTRIA3        11       1      11      12       8                        
CTRIA3        12       1       1      11       8                        
CTRIA3        13       2      10       9      12                        
CTRIA3        14       2      11      10      12                        
ENDDATA
