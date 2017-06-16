ID D:\Users\Dima\My_document,FEMAP
SOL SESTATICS
TIME 10000
CEND
  ECHO = NONE
  DISPLACEMENT(PLOT) = ALL
  SPCFORCE(PLOT) = ALL
  OLOAD(PLOT) = ALL
  FORCE(PLOT,CORNER) = ALL
  STRESS(PLOT,CORNER) = ALL
  LOAD = 1
BEGIN BULK
$ ***************************************************************************
$   Written by : FEMAP
$   Version    : 7.00
$   Translator : MSC/NASTRAN
$   From Model : D:\Users\Dima\My_documents\PAVT_2009\models\Pipe\Solid\Pipe_model.MOD
$   Date       : Sat Feb 14 20:34:53 2009
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
$ FEMAP Load Set 1 : Source Displn -.001
PLOAD4         1       3      1.                               8      10
$ FEMAP Load Set 2 : Absorbing Uptake 1.0
PLOAD4         2       8      2.                              27      32
$ FEMAP Load Set 3 : Boundary Skews
PLOAD4         3       3      3.                               8      12
PLOAD4         3       4      3.                               9      20
PLOAD4         3       1      3.                               1      37
PLOAD4         3       1      3.                               1       2
PLOAD4         3       2      3.                               5      35
PLOAD4         3       7      3.                               4      28
PLOAD4         3       8      3.                              29      32
PLOAD4         3       9      3.                               4      38
PLOAD4         3       3      3.                              14      13
PLOAD4         3       4      3.                              15      17
PLOAD4         3       3      3.                              11      13
PLOAD4         3       4      3.                              12      17
PLOAD4         3       3      3.                               8      15
PLOAD4         3       4      3.                               9      16
PLOAD4         3       1      3.                               3      39
PLOAD4         3       2      3.                               2      31
PLOAD4         3       1      3.                               1       7
PLOAD4         3       2      3.                               5       4
PLOAD4         3       5      3.                              16      18
PLOAD4         3       6      3.                              23      19
PLOAD4         3       5      3.                              20      18
PLOAD4         3       6      3.                              22      19
PLOAD4         3       5      3.                              40      23
PLOAD4         3       6      3.                              21      24
PLOAD4         3       6      3.                              36      24
PLOAD4         3       7      3.                              33      30
PLOAD4         3       8      3.                              34      25
PLOAD4         3       7      3.                              30      31
PLOAD4         3       8      3.                              25      28
PLOAD4         3       7      3.                              33      29
PLOAD4         3       8      3.                              34      26
PLOAD4         3       9      3.                              31      21
PLOAD4         3      10      3.                              39      40
PLOAD4         3       9      3.                               4      22
PLOAD4         3      10      3.                               7      20
PLOAD4         3      10      3.                              37      20
$ FEMAP Property 1 : Acoustic space
PSOLID         1       1       0        
$ FEMAP Property 2 : Rigid body
PSOLID         2       1       0        
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           1       0     1.2     -2.     0.5       0        
GRID           2       0     1.8     -2.    -0.5       0        
GRID           3       0     1.2     -2.    -0.5       0        
GRID           4       0     2.4     -1.     0.5       0        
GRID           5       0     1.8     -2.     0.5       0        
GRID           6       0     1.2     -1.     0.5       0        
GRID           7       0     1.8     -1.     0.5       0        
GRID           8       0      0.      0.    -0.5       0        
GRID           9       0     0.6      0.    -0.5       0        
GRID          10       0      0.      1.     0.5       0        
GRID          11       0      0.      0.     0.5       0        
GRID          12       0     0.6      0.     0.5       0        
GRID          13       0     0.6      1.     0.5       0        
GRID          14       0      0.      1.    -0.5       0        
GRID          15       0     0.6      1.    -0.5       0        
GRID          16       0     1.2      1.    -0.5       0        
GRID          17       0     1.2      1.     0.5       0        
GRID          18       0     1.8      1.     0.5       0        
GRID          19       0     2.4      1.     0.5       0        
GRID          20       0     1.2      0.     0.5       0        
GRID          21       0     1.8      0.    -0.5       0        
GRID          22       0     1.8      0.     0.5       0        
GRID          23       0     1.8      1.    -0.5       0        
GRID          24       0     2.4      1.    -0.5       0        
GRID          25       0     3.6     -2.    -0.5       0        
GRID          26       0     3.6     -1.     0.5       0        
GRID          27       0     3.6     -2.     0.5       0        
GRID          28       0      3.     -1.    -0.5       0        
GRID          29       0      3.     -1.     0.5       0        
GRID          30       0      3.     -2.    -0.5       0        
GRID          31       0     2.4     -1.    -0.5       0        
GRID          32       0     3.6     -1.    -0.5       0        
GRID          33       0     2.4     -2.     0.5       0        
GRID          34       0      3.     -2.     0.5       0        
GRID          35       0     2.4     -2.    -0.5       0        
GRID          36       0     2.4      0.     0.5       0        
GRID          37       0     1.2     -1.    -0.5       0        
GRID          38       0     2.4      0.    -0.5       0        
GRID          39       0     1.8     -1.    -0.5       0        
GRID          40       0     1.2      0.    -0.5       0        
CHEXA          1       2       1       3      37       6       5       2+EL    1
+EL    1      39       7                                                        
CHEXA          2       2       5       2      39       7      33      35+EL    2
+EL    2      31       4                                                        
CHEXA          3       2       8      11      12       9      14      10+EL    3
+EL    3      13      15                                                        
CHEXA          4       2       9      12      20      40      15      13+EL    4
+EL    4      17      16                                                        
CHEXA          5       2      40      20      22      21      16      17+EL    5
+EL    5      18      23                                                        
CHEXA          6       2      21      22      36      38      23      18+EL    6
+EL    6      19      24                                                        
CHEXA          7       2      33      34      30      35       4      29+EL    7
+EL    7      28      31                                                        
CHEXA          8       2      34      27      25      30      29      26+EL    8
+EL    8      32      28                                                        
CHEXA          9       2       4      31      39       7      36      38+EL    9
+EL    9      21      22                                                        
CHEXA         10       2       7      39      37       6      22      21+EL    A
+EL    A      40      20                                                        
ENDDATA
