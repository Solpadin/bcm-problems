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
  SPC = 1
  LOAD = 1
BEGIN BULK
$ ***************************************************************************
$   Written by : FEMAP
$   Version    : 8.00
$   Translator : NE/Nastran
$   From Model : 
$   Date       : Fri Feb 27 17:56:10 2004
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
$ FEMAP Load Set 1 : L
FORCE          1       5       0      1.   1000.      0.      0.
FORCE          1       6       0      1.   1000.      0.      0.
FORCE          1       7       0      1.   1000.      0.      0.
FORCE          1       8       0      1.   1000.      0.      0.
$ FEMAP Constraint Set 1 : NASTRAN SPC 1
SPC            1       1     123      0.
SPC            1       2     123      0.
SPC            1       3     123      0.
SPC            1       4     123      0.
$ FEMAP Property 1 : Untitled
PSOLID         1       1       0        
$ FEMAP Material 1 : Untitled
MAT1           1 7.2E+10 6.4E+10     0.3      0.      0.      0.        
$ FEMAP Material 2 : Untitled
MAT1           2 7.2E+10 6.4E+10     0.3      0.      0.      0.        
GRID           1       0      0.      0.      0.       0        
GRID           2       0      0.      0.      1.       0        
GRID           3       0      1.      0.      1.       0        
GRID           4       0      1.      0.      0.       0        
GRID           5       0      0.      1.      0.       0        
GRID           6       0      0.      1.      1.       0        
GRID           7       0      1.      1.      1.       0        
GRID           8       0      1.      1.      0.       0        
CTETRA         1       1       1       2       4       5                        
CTETRA         2       1       2       3       4       7                        
CTETRA         3       1       5       2       7       6                        
CTETRA         4       1       4       8       5       7                        
CTETRA         5       1       5       4       7       2                        
ENDDATA
