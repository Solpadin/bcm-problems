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
$   Date       : Fri Feb 27 18:00:37 2004
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
FORCE          1      29       0      1.   1000.      0.      0.
FORCE          1      30       0      1.   1000.      0.      0.
FORCE          1      31       0      1.   1000.      0.      0.
FORCE          1      32       0      1.   1000.      0.      0.
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
GRID          13       0      0.      2.      0.       0        
GRID          14       0      0.      2.      1.       0        
GRID          15       0      1.      2.      1.       0        
GRID          16       0      1.      2.      0.       0        
GRID          17       0      0.      3.      0.       0        
GRID          18       0      0.      3.      1.       0        
GRID          19       0      1.      3.      1.       0        
GRID          20       0      1.      3.      0.       0        
GRID          24       0      1.      4.      0.       0        
GRID          25       0      0.      4.      0.       0        
GRID          26       0      0.      4.      1.       0        
GRID          27       0      1.      4.      1.       0        
GRID          29       0      0.      5.      0.       0        
GRID          30       0      0.      5.      1.       0        
GRID          31       0      1.      5.      1.       0        
GRID          32       0      1.      5.      0.       0        
CTETRA         1       1       1       2       4       5                        
CTETRA         2       1       2       3       4       7                        
CTETRA         3       1       5       2       7       6                        
CTETRA         4       1       4       8       5       7                        
CTETRA         5       1       5       4       7       2                        
CTETRA         6       1       5       6       7      14                        
CTETRA         7       1       5      16      13      14                        
CTETRA         8       1       7      16      14      15                        
CTETRA         9       1      14       5      16       7                        
CTETRA        10       1      16       7       5       8                        
CTETRA        11       1      13      14      16      17                        
CTETRA        12       1      14      15      16      19                        
CTETRA        13       1      17      14      19      18                        
CTETRA        14       1      16      20      17      19                        
CTETRA        15       1      17      16      19      14                        
CTETRA        16       1      17      18      19      26                        
CTETRA        17       1      17      24      25      26                        
CTETRA        18       1      19      24      26      27                        
CTETRA        19       1      26      17      24      19                        
CTETRA        20       1      24      19      17      20                        
CTETRA        21       1      25      26      24      29                        
CTETRA        22       1      26      27      24      31                        
CTETRA        23       1      29      26      31      30                        
CTETRA        24       1      24      32      29      31                        
CTETRA        25       1      29      24      31      26                        
ENDDATA
