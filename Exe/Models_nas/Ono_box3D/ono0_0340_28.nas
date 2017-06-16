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
  LOAD = 1
BEGIN BULK
$ ***************************************************************************
$   Written by : FEMAP
$   Version    : 7.00
$   Translator : MSC/NASTRAN
$   From Model : 
$   Date       : Mon Dec 22 15:43:21 2008
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
$ FEMAP Load Set 1 : Srcplane Vn -1.0 0.0
PLOAD4         1       8    111.                              52      22
$ FEMAP Load Set 2 : Endwall Absorb 1.0 0.0
PLOAD4         2       1    222.                               1      63
PLOAD4         2      14    222.                              63      28
PLOAD4         2      21    222.                              63      39
PLOAD4         2      22    222.                              63      59
$ FEMAP Property 1 : Untitled
PSHELL         1       1      0.       1      1.       1 0.83333      0.
$ FEMAP Property 2 : Untitled
PSHELL         2       1      0.       1      1.       1 0.83333      0.
$ FEMAP Property 3 : Untitled
PSOLID         3       1       0        
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           1       0      0.      0.     2.5       0        
GRID           2       0      0.    0.25 1.07143       0        
GRID           3       0      0.    0.25      0.       0        
GRID           4       0      0.      0.      0.       0        
GRID           5       0      0.      0. 0.35714       0        
GRID           6       0      0.      0. 1.42857       0        
GRID           7       0      0.      0. 2.14286       0        
GRID           8       0    0.35      0. 1.07143       0        
GRID           9       0    0.35      0. 2.14286       0        
GRID          10       0      0.      0. 1.78571       0        
GRID          11       0      0.      0. 1.07143       0        
GRID          12       0      0.      0. 0.71429       0        
GRID          13       0    0.35      0. 0.71429       0        
GRID          14       0      0.    0.25 1.78571       0        
GRID          15       0    0.35     0.5 2.14286       0        
GRID          16       0    0.35    0.25 1.42857       0        
GRID          17       0    0.35    0.25 1.78571       0        
GRID          18       0     0.6    0.25 0.71429       0        
GRID          19       0     0.6    0.25 1.07143       0        
GRID          20       0    0.35    0.25 1.07143       0        
GRID          21       0    0.35    0.25 0.71429       0        
GRID          22       0     0.6     0.5      0.       0        
GRID          23       0     0.6     0.5 0.71429       0        
GRID          24       0     0.6     0.5 1.07143       0        
GRID          25       0     0.6    0.25 2.14286       0        
GRID          26       0     0.6    0.25 0.35714       0        
GRID          27       0    0.35     0.5     2.5       0        
GRID          28       0     0.6     0.5     2.5       0        
GRID          29       0     0.6     0.5 2.14286       0        
GRID          30       0     0.6     0.5 1.78571       0        
GRID          31       0     0.6     0.5 1.42857       0        
GRID          32       0     0.6     0.5 0.35714       0        
GRID          33       0    0.35     0.5 0.71429       0        
GRID          34       0    0.35     0.5 1.78571       0        
GRID          35       0     0.6    0.25     2.5       0        
GRID          36       0     0.6    0.25 1.42857       0        
GRID          37       0     0.6    0.25      0.       0        
GRID          38       0     0.6    0.25 1.78571       0        
GRID          39       0     0.6      0.     2.5       0        
GRID          40       0     0.6      0. 2.14286       0        
GRID          41       0     0.6      0. 1.78571       0        
GRID          42       0     0.6      0. 1.07143       0        
GRID          43       0     0.6      0. 0.71429       0        
GRID          44       0    0.35      0.      0.       0        
GRID          45       0    0.35      0. 0.35714       0        
GRID          46       0    0.35      0. 1.42857       0        
GRID          47       0    0.35    0.25 2.14286       0        
GRID          48       0    0.35      0. 1.78571       0        
GRID          49       0     0.6      0. 0.35714       0        
GRID          50       0     0.6      0. 1.42857       0        
GRID          51       0     0.6      0.      0.       0        
GRID          52       0    0.35    0.25      0.       0        
GRID          53       0    0.35      0.     2.5       0        
GRID          54       0    0.35     0.5      0.       0        
GRID          55       0      0.     0.5      0.       0        
GRID          56       0      0.     0.5 1.07143       0        
GRID          57       0      0.     0.5 1.42857       0        
GRID          58       0      0.     0.5 2.14286       0        
GRID          59       0      0.     0.5     2.5       0        
GRID          60       0    0.35     0.5 1.42857       0        
GRID          61       0    0.35     0.5 1.07143       0        
GRID          62       0    0.35     0.5 0.35714       0        
GRID          63       0    0.35    0.25     2.5       0        
GRID          64       0      0.    0.25     2.5       0        
GRID          65       0      0.    0.25 2.14286       0        
GRID          66       0      0.    0.25 0.71429       0        
GRID          67       0    0.35    0.25 0.35714       0        
GRID          68       0      0.    0.25 0.35714       0        
GRID          69       0      0.    0.25 1.42857       0        
GRID          70       0      0.     0.5 1.78571       0        
GRID          71       0      0.     0.5 0.71429       0        
GRID          72       0      0.     0.5 0.35714       0        
CHEXA          1       3       1       7      65      64      53       9+EL    1
+EL    1      47      63                                                        
CHEXA          2       3       7      10      14      65       9      48+EL    2
+EL    2      17      47                                                        
CHEXA          3       3      10       6      69      14      48      46+EL    3
+EL    3      16      17                                                        
CHEXA          4       3       6      11       2      69      46       8+EL    4
+EL    4      20      16                                                        
CHEXA          5       3      11      12      66       2       8      13+EL    5
+EL    5      21      20                                                        
CHEXA          6       3      12       5      68      66      13      45+EL    6
+EL    6      67      21                                                        
CHEXA          7       3       5       4       3      68      45      44+EL    7
+EL    7      52      67                                                        
CHEXA          8       3      52      67      26      37      54      62+EL    8
+EL    8      32      22                                                        
CHEXA          9       3      67      21      18      26      62      33+EL    9
+EL    9      23      32                                                        
CHEXA         10       3      21      20      19      18      33      61+EL    A
+EL    A      24      23                                                        
CHEXA         11       3      20      16      36      19      61      60+EL    B
+EL    B      31      24                                                        
CHEXA         12       3      16      17      38      36      60      34+EL    C
+EL    C      30      31                                                        
CHEXA         13       3      17      47      25      38      34      15+EL    D
+EL    D      29      30                                                        
CHEXA         14       3      47      63      35      25      15      27+EL    E
+EL    E      28      29                                                        
CHEXA         15       3      44      52      67      45      51      37+EL    F
+EL    F      26      49                                                        
CHEXA         16       3      45      67      21      13      49      26+EL    G
+EL    G      18      43                                                        
CHEXA         17       3      13      21      20       8      43      18+EL    H
+EL    H      19      42                                                        
CHEXA         18       3       8      20      16      46      42      19+EL    I
+EL    I      36      50                                                        
CHEXA         19       3      46      16      17      48      50      36+EL    J
+EL    J      38      41                                                        
CHEXA         20       3      48      17      47       9      41      38+EL    K
+EL    K      25      40                                                        
CHEXA         21       3       9      47      63      53      40      25+EL    L
+EL    L      35      39                                                        
CHEXA         22       3      63      47      65      64      27      15+EL    M
+EL    M      58      59                                                        
CHEXA         23       3      47      17      14      65      15      34+EL    N
+EL    N      70      58                                                        
CHEXA         24       3      17      16      69      14      34      60+EL    O
+EL    O      57      70                                                        
CHEXA         25       3      16      20       2      69      60      61+EL    P
+EL    P      56      57                                                        
CHEXA         26       3      20      21      66       2      61      33+EL    Q
+EL    Q      71      56                                                        
CHEXA         27       3      21      67      68      66      33      62+EL    R
+EL    R      72      71                                                        
CHEXA         28       3      67      52       3      68      62      54+EL    S
+EL    S      55      72                                                        
ENDDATA
