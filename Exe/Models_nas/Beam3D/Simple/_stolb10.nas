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
$   Date       : Sat Jul 15 17:49:38 2006
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
PBEAM          1       1      0.      0.      0.      0.      0.      0.+PR    1
+PR    1      0.      0.      0.      0.      0.      0.      0.      0.+PA    1
+PA    1    YESA      1.                                                +PC    1
+PC    1                                                                        
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           2       0      0.    1.25      0.       0        
GRID           3       0      0.     1.5      0.       0        
GRID           4       0    0.25      0.      0.       0        
GRID           5       0  0.1875    0.25      0.       0        
GRID           6       0   0.125     0.5      0.       0        
GRID           7       0  0.0625    0.75      0.       0        
GRID           9       0   -0.25      0.      0.       0        
GRID          10       0 -0.1875    0.25      0.       0        
GRID          11       0  -0.125     0.5      0.       0        
GRID          12       0 -0.0625    0.75      0.       0        
GRID          13       0      0.      1.      0.       0        
CBEAM          1       1      13       2      1.      0.      0.
CBEAM          2       1       2       3      1.      0.      0.
CBEAM          3       1       4       5      1.      0.      0.
CBEAM          4       1       5       6      1.      0.      0.
CBEAM          5       1       6       7      1.      0.      0.
CBEAM          6       1       7      13      1.      0.      0.
CBEAM          7       1       9      10      1.      0.      0.
CBEAM          8       1      10      11      1.      0.      0.
CBEAM          9       1      11      12      1.      0.      0.
CBEAM         10       1      12      13      1.      0.      0.
ENDDATA
