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
$   Date       : Mon Oct 05 01:12:22 2009
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
PSHELL         1       1      0.       1      1.       1 0.83333      0.
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           1       0      0.      0.      0.       0        
GRID           2       0      1.      0.      0.       0        
GRID           3       0      1.     0.5      0.       0        
GRID           4       0      0.    0.25      0.       0        
GRID           5       0      1.    0.25      0.       0        
GRID           6       0      0.    0.75      0.       0        
GRID           7       0      0.     0.5      0.       0        
GRID           8       0      1.      1.      0.       0        
GRID           9       0      0.      1.      0.       0        
GRID          10       0      1.    0.75      0.       0        
GRID          11       0      1.    1.25      0.       0        
GRID          12       0      0.    1.25      0.       0        
GRID          13       0      1.     1.5      0.       0        
GRID          14       0      0.     1.5      0.       0        
GRID          15       0      1.      2.      0.       0        
GRID          16       0      0.      2.      0.       0        
GRID          17       0      0.    1.75      0.       0        
GRID          18       0      1.    1.75      0.       0        
CQUAD4         1       1       5       4       1       2                
CQUAD4         2       1       3       7       4       5                
CQUAD4         3       1      10       6       7       3                
CQUAD4         4       1       8       9       6      10                
CQUAD4         5       1      11      12       9       8                
CQUAD4         6       1      13      14      12      11                
CQUAD4         7       1      18      17      14      13                
CQUAD4         8       1      15      16      17      18                
ENDDATA
