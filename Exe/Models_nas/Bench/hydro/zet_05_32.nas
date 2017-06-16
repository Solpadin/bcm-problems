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
$   From Model : D:\Job\Dima\Geo_2009\Bench\Solid\Zet_model.MOD
$   Date       : Wed Oct 13 19:18:10 2010
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
PLOAD4         1       9      1.                              19      22
PLOAD4         1      11      1.                              20      25
PLOAD4         1      13      1.                              28      17
PLOAD4         1      15      1.                              22      16
$ FEMAP Load Set 2 : End Pressure 0
PLOAD4         2      26      2.                              63      60
PLOAD4         2      28      2.                              62      64
PLOAD4         2      30      2.                              64      66
PLOAD4         2      32      2.                              69      61
$ FEMAP Property 1 : Untitled
PSOLID         1       1       0        
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           1       0     -1.    -0.5     1.8       0        
GRID           2       0     -1.     0.5     1.8       0        
GRID           3       0     -1.     0.5     1.2       0        
GRID           4       0     -1.      0.     1.8       0        
GRID           5       0    -0.5      0.     2.4       0        
GRID           6       0    -0.5    -0.5     1.2       0        
GRID           7       0    -0.5    -0.5     1.8       0        
GRID           8       0    -0.5     0.5     1.2       0        
GRID           9       0    -0.5     0.5     1.8       0        
GRID          10       0     -1.    -0.5     1.2       0        
GRID          11       0     -1.      0.     1.2       0        
GRID          12       0    -0.5      0.     1.2       0        
GRID          13       0      0.    -0.5     0.6       0        
GRID          14       0      0.     0.5     1.2       0        
GRID          15       0      0.      0.     0.6       0        
GRID          16       0      1.     0.5      0.       0        
GRID          17       0      1.      0.      0.       0        
GRID          18       0      1.      0.     0.6       0        
GRID          19       0      0.    -0.5      0.       0        
GRID          20       0      0.      0.      0.       0        
GRID          21       0      0.     0.5      0.       0        
GRID          22       0     0.5      0.      0.       0        
GRID          23       0      0.     0.5     0.6       0        
GRID          24       0      1.     0.5     0.6       0        
GRID          25       0     0.5     0.5      0.       0        
GRID          26       0     0.5     0.5     0.6       0        
GRID          27       0      0.    -0.5     1.2       0        
GRID          28       0     0.5    -0.5      0.       0        
GRID          29       0      1.    -0.5      0.       0        
GRID          30       0      1.    -0.5     0.6       0        
GRID          31       0     0.5    -0.5     0.6       0        
GRID          32       0      1.     0.5     1.2       0        
GRID          33       0     0.5      0.     1.2       0        
GRID          34       0      0.      0.     2.4       0        
GRID          35       0      0.     0.5     2.4       0        
GRID          36       0      0.      0.     1.8       0        
GRID          37       0      1.      0.     1.2       0        
GRID          38       0      1.     0.5     1.8       0        
GRID          39       0      1.    -0.5     1.8       0        
GRID          40       0      1.      0.     1.8       0        
GRID          41       0      0.      0.     1.2       0        
GRID          42       0     0.5     0.5     1.2       0        
GRID          43       0      1.    -0.5     1.2       0        
GRID          44       0      0.     0.5     1.8       0        
GRID          45       0     0.5     0.5     2.4       0        
GRID          46       0     0.5     0.5     1.8       0        
GRID          47       0     0.5    -0.5     2.4       0        
GRID          48       0      0.    -0.5     1.8       0        
GRID          49       0     0.5    -0.5     1.2       0        
GRID          50       0     0.5    -0.5     1.8       0        
GRID          51       0      1.    -0.5     2.4       0        
GRID          52       0      1.      0.     2.4       0        
GRID          53       0      1.     0.5     2.4       0        
GRID          54       0     0.5      0.     2.4       0        
GRID          55       0     -1.    -0.5      3.       0        
GRID          56       0     -1.    -0.5     3.6       0        
GRID          57       0     -1.     0.5      3.       0        
GRID          58       0     -1.     0.5     2.4       0        
GRID          59       0     -1.      0.      3.       0        
GRID          60       0    -0.5    -0.5     3.6       0        
GRID          61       0      0.      0.     3.6       0        
GRID          62       0     -1.     0.5     3.6       0        
GRID          63       0     -1.      0.     3.6       0        
GRID          64       0    -0.5      0.     3.6       0        
GRID          65       0      0.    -0.5     2.4       0        
GRID          66       0      0.    -0.5     3.6       0        
GRID          67       0    -0.5    -0.5     2.4       0        
GRID          68       0    -0.5    -0.5      3.       0        
GRID          69       0    -0.5     0.5     3.6       0        
GRID          70       0      0.     0.5      3.       0        
GRID          71       0    -0.5     0.5     2.4       0        
GRID          72       0    -0.5     0.5      3.       0        
GRID          73       0     -1.    -0.5     2.4       0        
GRID          74       0     -1.      0.     2.4       0        
GRID          75       0      0.    -0.5      3.       0        
GRID          76       0      0.     0.5     3.6       0        
GRID          77       0      0.      0.      3.       0        
GRID          78       0    -0.5      0.     1.8       0        
GRID          79       0     0.5      0.     0.6       0        
GRID          80       0     0.5      0.     1.8       0        
GRID          81       0    -0.5      0.      3.       0        
CHEXA          1       1      10      11       4       1       6      12+EL    1
+EL    1      78       7                                                        
CHEXA          2       1       1       4      74      73       7      78+EL    2
+EL    2       5      67                                                        
CHEXA          3       1      11       3       2       4      12       8+EL    3
+EL    3       9      78                                                        
CHEXA          4       1       4       2      58      74      78       9+EL    4
+EL    4      71       5                                                        
CHEXA          5       1       6      12      78       7      27      41+EL    5
+EL    5      36      48                                                        
CHEXA          6       1       7      78       5      67      48      36+EL    6
+EL    6      34      65                                                        
CHEXA          7       1      12       8       9      78      41      14+EL    7
+EL    7      44      36                                                        
CHEXA          8       1      78       9      71       5      36      44+EL    8
+EL    8      35      34                                                        
CHEXA          9       1      19      20      15      13      28      22+EL    9
+EL    9      79      31                                                        
CHEXA         10       1      13      15      41      27      31      79+EL    A
+EL    A      33      49                                                        
CHEXA         11       1      20      21      23      15      22      25+EL    B
+EL    B      26      79                                                        
CHEXA         12       1      15      23      14      41      79      26+EL    C
+EL    C      42      33                                                        
CHEXA         13       1      28      22      79      31      29      17+EL    D
+EL    D      18      30                                                        
CHEXA         14       1      31      79      33      49      30      18+EL    E
+EL    E      37      43                                                        
CHEXA         15       1      22      25      26      79      17      16+EL    F
+EL    F      24      18                                                        
CHEXA         16       1      79      26      42      33      18      24+EL    G
+EL    G      32      37                                                        
CHEXA         17       1      27      41      36      48      49      33+EL    H
+EL    H      80      50                                                        
CHEXA         18       1      48      36      34      65      50      80+EL    I
+EL    I      54      47                                                        
CHEXA         19       1      41      14      44      36      33      42+EL    J
+EL    J      46      80                                                        
CHEXA         20       1      36      44      35      34      80      46+EL    K
+EL    K      45      54                                                        
CHEXA         21       1      49      33      80      50      43      37+EL    L
+EL    L      40      39                                                        
CHEXA         22       1      50      80      54      47      39      40+EL    M
+EL    M      52      51                                                        
CHEXA         23       1      33      42      46      80      37      32+EL    N
+EL    N      38      40                                                        
CHEXA         24       1      80      46      45      54      40      38+EL    O
+EL    O      53      52                                                        
CHEXA         25       1      73      74      59      55      67       5+EL    P
+EL    P      81      68                                                        
CHEXA         26       1      55      59      63      56      68      81+EL    Q
+EL    Q      64      60                                                        
CHEXA         27       1      74      58      57      59       5      71+EL    R
+EL    R      72      81                                                        
CHEXA         28       1      59      57      62      63      81      72+EL    S
+EL    S      69      64                                                        
CHEXA         29       1      67       5      81      68      65      34+EL    T
+EL    T      77      75                                                        
CHEXA         30       1      68      81      64      60      75      77+EL    U
+EL    U      61      66                                                        
CHEXA         31       1       5      71      72      81      34      35+EL    V
+EL    V      70      77                                                        
CHEXA         32       1      81      72      69      64      77      70+EL    W
+EL    W      76      61                                                        
ENDDATA
