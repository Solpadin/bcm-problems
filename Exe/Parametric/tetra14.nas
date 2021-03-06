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
$   Date       : Sun Feb 25 08:10:21 2007
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
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           1       0      1.      0.      0.       0        
GRID           2       0      0.      0.      1.       0        
GRID           3       0      0.      1.      0.       0        
GRID           4       0      0.      0.      0.       0        
GRID           5       0     0.5      0.     0.5       0        
GRID           6       0     0.5     0.5      0.       0        
GRID           7       0     0.5      0.      0.       0        
GRID           8       0      0.     0.5     0.5       0        
GRID           9       0      0.      0.     0.5       0        
GRID          10       0      0.     0.5      0.       0        
GRID          11       0 0.33333 0.33333 0.33333       0        
GRID          12       0 0.33333 0.33333      0.       0        
GRID          13       0 0.33333      0. 0.33333       0        
GRID          14       0      0. 0.33333 0.33333       0        
CTETRA         1       1       2       3       1       4       8       6+EL    1
+EL    1       5       9      10       7
ENDDATA
