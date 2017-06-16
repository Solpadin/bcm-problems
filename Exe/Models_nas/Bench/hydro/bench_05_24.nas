ID D:\Job\Dima\Geo_2009\Benc,FEMAP
SOL SESTATICS
TIME 10000
CEND
  ECHO = NONE
  DISPLACEMENT(PLOT) = ALL
  SPCFORCE(PLOT) = ALL
  OLOAD(PLOT) = ALL
  FORCE(PLOT,CORNER) = ALL
  STRESS(PLOT,CORNER) = ALL
  LOAD = 2
BEGIN BULK
$ ***************************************************************************
$   Written by : FEMAP
$   Version    : 7.00
$   Translator : MSC/NASTRAN
$   From Model : D:\Job\Dima\Geo_2009\Bench\Solid\Bench_model.MOD
$   Date       : Wed Oct 13 19:15:42 2010
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
$ FEMAP Load Set 1 : Source Pressure 1
PLOAD4         1       9      1.                              28      31
PLOAD4         1      11      1.                              29      34
PLOAD4         1      13      1.                              38      26
PLOAD4         1      15      1.                              31      24
$ FEMAP Load Set 2 : End Pressure 0
PLOAD4         2       1      2.                              16       5
PLOAD4         2       2      2.                               1       9
PLOAD4         2       3      2.                              17       3
PLOAD4         2       4      2.                               5       8
$ FEMAP Property 1 : Untitled
PSOLID         1       1       0        
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           1       0     -1.    -0.5     1.8       0        
GRID           2       0     -1.    -0.5     2.4       0        
GRID           3       0     -1.     0.5     1.8       0        
GRID           4       0     -1.     0.5     1.2       0        
GRID           5       0     -1.      0.     1.8       0        
GRID           6       0    -0.5    -0.5     2.4       0        
GRID           7       0      0.      0.     2.4       0        
GRID           8       0     -1.     0.5     2.4       0        
GRID           9       0     -1.      0.     2.4       0        
GRID          10       0    -0.5      0.     2.4       0        
GRID          11       0    -0.5    -0.5     1.2       0        
GRID          12       0    -0.5    -0.5     1.8       0        
GRID          13       0    -0.5     0.5     2.4       0        
GRID          14       0    -0.5     0.5     1.2       0        
GRID          15       0    -0.5     0.5     1.8       0        
GRID          16       0     -1.    -0.5     1.2       0        
GRID          17       0     -1.      0.     1.2       0        
GRID          18       0    -0.5      0.     1.2       0        
GRID          19       0      0.     0.5     2.4       0        
GRID          20       0      0.      0.     1.8       0        
GRID          21       0      0.    -0.5     0.6       0        
GRID          22       0      0.     0.5     1.2       0        
GRID          23       0      0.      0.     0.6       0        
GRID          24       0      1.     0.5      0.       0        
GRID          25       0      1.      0.     1.2       0        
GRID          26       0      1.      0.      0.       0        
GRID          27       0      1.      0.     0.6       0        
GRID          28       0      0.    -0.5      0.       0        
GRID          29       0      0.      0.      0.       0        
GRID          30       0      0.     0.5      0.       0        
GRID          31       0     0.5      0.      0.       0        
GRID          32       0      0.     0.5     0.6       0        
GRID          33       0      1.     0.5     0.6       0        
GRID          34       0     0.5     0.5      0.       0        
GRID          35       0     0.5     0.5     0.6       0        
GRID          36       0      1.    -0.5     1.2       0        
GRID          37       0      0.    -0.5     1.2       0        
GRID          38       0     0.5    -0.5      0.       0        
GRID          39       0      1.    -0.5      0.       0        
GRID          40       0      1.    -0.5     0.6       0        
GRID          41       0     0.5    -0.5     0.6       0        
GRID          42       0     0.5      0.     1.2       0        
GRID          43       0      1.     0.5     1.8       0        
GRID          44       0      1.    -0.5     1.8       0        
GRID          45       0      1.      0.     1.8       0        
GRID          46       0      0.      0.     1.2       0        
GRID          47       0     0.5     0.5     1.2       0        
GRID          48       0      0.     0.5     1.8       0        
GRID          49       0     0.5     0.5     2.4       0        
GRID          50       0      1.     0.5     1.2       0        
GRID          51       0     0.5     0.5     1.8       0        
GRID          52       0     0.5    -0.5     2.4       0        
GRID          53       0      0.    -0.5     1.8       0        
GRID          54       0     0.5    -0.5     1.2       0        
GRID          55       0     0.5    -0.5     1.8       0        
GRID          56       0      0.    -0.5     2.4       0        
GRID          57       0      1.    -0.5     2.4       0        
GRID          58       0      1.      0.     2.4       0        
GRID          59       0      1.     0.5     2.4       0        
GRID          60       0     0.5      0.     2.4       0        
GRID          61       0    -0.5      0.     1.8       0        
GRID          62       0     0.5      0.     0.6       0        
GRID          63       0     0.5      0.     1.8       0        
CHEXA          1       1      16      17       5       1      11      18+EL    1
+EL    1      61      12                                                        
CHEXA          2       1       1       5       9       2      12      61+EL    2
+EL    2      10       6                                                        
CHEXA          3       1      17       4       3       5      18      14+EL    3
+EL    3      15      61                                                        
CHEXA          4       1       5       3       8       9      61      15+EL    4
+EL    4      13      10                                                        
CHEXA          5       1      11      18      61      12      37      46+EL    5
+EL    5      20      53                                                        
CHEXA          6       1      12      61      10       6      53      20+EL    6
+EL    6       7      56                                                        
CHEXA          7       1      18      14      15      61      46      22+EL    7
+EL    7      48      20                                                        
CHEXA          8       1      61      15      13      10      20      48+EL    8
+EL    8      19       7                                                        
CHEXA          9       1      28      29      23      21      38      31+EL    9
+EL    9      62      41                                                        
CHEXA         10       1      21      23      46      37      41      62+EL    A
+EL    A      42      54                                                        
CHEXA         11       1      29      30      32      23      31      34+EL    B
+EL    B      35      62                                                        
CHEXA         12       1      23      32      22      46      62      35+EL    C
+EL    C      47      42                                                        
CHEXA         13       1      38      31      62      41      39      26+EL    D
+EL    D      27      40                                                        
CHEXA         14       1      41      62      42      54      40      27+EL    E
+EL    E      25      36                                                        
CHEXA         15       1      31      34      35      62      26      24+EL    F
+EL    F      33      27                                                        
CHEXA         16       1      62      35      47      42      27      33+EL    G
+EL    G      50      25                                                        
CHEXA         17       1      37      46      20      53      54      42+EL    H
+EL    H      63      55                                                        
CHEXA         18       1      53      20       7      56      55      63+EL    I
+EL    I      60      52                                                        
CHEXA         19       1      46      22      48      20      42      47+EL    J
+EL    J      51      63                                                        
CHEXA         20       1      20      48      19       7      63      51+EL    K
+EL    K      49      60                                                        
CHEXA         21       1      54      42      63      55      36      25+EL    L
+EL    L      45      44                                                        
CHEXA         22       1      55      63      60      52      44      45+EL    M
+EL    M      58      57                                                        
CHEXA         23       1      42      47      51      63      25      50+EL    N
+EL    N      43      45                                                        
CHEXA         24       1      63      51      49      60      45      43+EL    O
+EL    O      59      58                                                        
ENDDATA
