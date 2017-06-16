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
$   Date       : Mon Mar 01 18:02:09 2004
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
FORCE          1      57       0      1.   1000.      0.      0.
FORCE          1      58       0      1.   1000.      0.      0.
FORCE          1      59       0      1.   1000.      0.      0.
FORCE          1      60       0      1.   1000.      0.      0.
$ FEMAP Constraint Set 1 : C
SPC            1       1     123      0.
SPC            1       2     123      0.
SPC            1       3     123      0.
SPC            1       4     123      0.
$ FEMAP Property 1 : solid
PSOLID         1       1       0        
$ FEMAP Material 1 : amg
MAT1           1 7.2E+10             0.3      0.      0.      0.        
GRID           1       0      0.      0.      0.       0        
GRID           2       0      0.      0.      1.       0        
GRID           3       0      1.      0.      1.       0        
GRID           4       0      1.      0.      0.       0        
GRID           5       0      0.      1.      0.       0        
GRID           7       0      1.      1.      1.       0        
GRID          10       0      0.      1.      1.       0        
GRID          12       0      1.      1.      0.       0        
GRID          17       0      0.      2.      0.       0        
GRID          18       0      0.      2.      1.       0        
GRID          19       0      1.      2.      1.       0        
GRID          20       0      1.      2.      0.       0        
GRID          25       0      0.      3.      0.       0        
GRID          26       0      0.      3.      1.       0        
GRID          27       0      1.      3.      1.       0        
GRID          28       0      1.      3.      0.       0        
GRID          33       0      0.      4.      0.       0        
GRID          34       0      0.      4.      1.       0        
GRID          35       0      1.      4.      1.       0        
GRID          36       0      1.      4.      0.       0        
GRID          41       0      0.      5.      0.       0        
GRID          42       0      0.      5.      1.       0        
GRID          43       0      1.      5.      1.       0        
GRID          44       0      1.      5.      0.       0        
GRID          49       0      0.      6.      0.       0        
GRID          50       0      0.      6.      1.       0        
GRID          51       0      1.      6.      1.       0        
GRID          52       0      1.      6.      0.       0        
GRID          57       0      0.      7.      0.       0        
GRID          58       0      0.      7.      1.       0        
GRID          59       0      1.      7.      1.       0        
GRID          60       0      1.      7.      0.       0        
CHEXA          1       1       1       2       3       4       5      10+EL    1
+EL    1       7      12                                                        
CHEXA          2       1       5      10       7      12      17      18+EL    2
+EL    2      19      20                                                        
CHEXA          3       1      17      18      19      20      25      26+EL    3
+EL    3      27      28                                                        
CHEXA          4       1      25      26      27      28      33      34+EL    4
+EL    4      35      36                                                        
CHEXA          5       1      33      34      35      36      41      42+EL    5
+EL    5      43      44                                                        
CHEXA          6       1      41      42      43      44      49      50+EL    6
+EL    6      51      52                                                        
CHEXA          7       1      49      50      51      52      57      58+EL    7
+EL    7      59      60                                                        
ENDDATA
