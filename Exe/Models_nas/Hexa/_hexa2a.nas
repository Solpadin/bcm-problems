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
$   Date       : Sun Jan 13 10:02:58 2008
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
PSOLID         1       1       0        
$ FEMAP Property 2 : Another Phase
PSOLID         2       1       0        
$ FEMAP Material 1 : Untitled
MAT1           1 7.2E+102.77E+10     0.3      0.      0.      0.        
GRID           1       0      1.      1.      0.       0        
GRID           2       0      0.      1.      0.       0        
GRID           3       0      0.      0.      0.       0        
GRID           4       0      1.      0.      0.       0        
GRID           5       0      1.      1.      1.       0        
GRID           6       0      0.      0.     1.1       0        
GRID           7       0      1.      1.     2.1       0        
GRID           8       0      0.      1.     2.1       0        
GRID           9       0      0.      0.     2.1       0        
GRID          10       0      1.      0.     2.1       0        
GRID          11       0      0.      1.      1.       0        
GRID          12       0      0.      0.      1.       0        
GRID          13       0      1.      0.      1.       0        
GRID          14       0      1.      1.     1.1       0        
GRID          15       0      0.      1.     1.1       0        
GRID          16       0      1.      0.     1.1       0        
CHEXA          1       1       1       2       3       4       5      11+EL    1
+EL    1      12      13                                                        
CHEXA          2       2       5      11      12      13      14      15+EL    2
+EL    2       6      16                                                        
CHEXA          3       1      14      15       6      16       7       8+EL    3
+EL    3       9      10                                                        
ENDDATA
