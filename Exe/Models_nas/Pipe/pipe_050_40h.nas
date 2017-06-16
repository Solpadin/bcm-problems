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
$   Date       : Sat Feb 14 17:04:19 2009
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
PLOAD4         1       9      1.                              25      28
PLOAD4         1      11      1.                              26      32
PLOAD4         1      13      1.                              36      23
PLOAD4         1      15      1.                              28      21
$ FEMAP Load Set 2 : Absorbing Uptake 1.0
PLOAD4         2      27      2.                              74      61
PLOAD4         2      28      2.                              69      59
PLOAD4         2      31      2.                              75      60
PLOAD4         2      32      2.                              61      62
$ FEMAP Load Set 3 : Boundary Skews
PLOAD4         3       9      3.                              25      20
PLOAD4         3      10      3.                              17      18
PLOAD4         3      11      3.                              26      29
PLOAD4         3      12      3.                              20      19
PLOAD4         3       1      3.                              11       2
PLOAD4         3       2      3.                              10      15
PLOAD4         3       3      3.                               1       8
PLOAD4         3       4      3.                               2      14
PLOAD4         3       1      3.                              11       5
PLOAD4         3       3      3.                               1       7
PLOAD4         3       5      3.                               4      79
PLOAD4         3       7      3.                               5      13
PLOAD4         3      29      3.                              89      66
PLOAD4         3      30      3.                              64      63
PLOAD4         3      31      3.                              65      60
PLOAD4         3      32      3.                              66      62
PLOAD4         3      33      3.                              89      82
PLOAD4         3      35      3.                              64      84
PLOAD4         3      37      3.                              90      81
PLOAD4         3      39      3.                              82      55
PLOAD4         3      13      3.                              37      24
PLOAD4         3      14      3.                              38      22
PLOAD4         3      15      3.                              23      31
PLOAD4         3      16      3.                              24      41
PLOAD4         3      11      3.                              27      33
PLOAD4         3      12      3.                              29      30
PLOAD4         3      15      3.                              32      31
PLOAD4         3      16      3.                              33      41
PLOAD4         3       9      3.                              25      39
PLOAD4         3      10      3.                              17      40
PLOAD4         3      13      3.                              36      38
PLOAD4         3      14      3.                              39      34
PLOAD4         3       3      3.                               3       9
PLOAD4         3       4      3.                               8       6
PLOAD4         3       7      3.                               7      71
PLOAD4         3       8      3.                               9      85
PLOAD4         3       1      3.                              11      12
PLOAD4         3       2      3.                              10      88
PLOAD4         3       5      3.                               4      78
PLOAD4         3       6      3.                              12      89
PLOAD4         3      21      3.                              34      43
PLOAD4         3      22      3.                              53      42
PLOAD4         3      23      3.                              22      48
PLOAD4         3      24      3.                              43      56
PLOAD4         3      19      3.                              19      49
PLOAD4         3      20      3.                              46      47
PLOAD4         3      23      3.                              30      48
PLOAD4         3      24      3.                              49      56
PLOAD4         3      17      3.                              35      54
PLOAD4         3      18      3.                              52      51
PLOAD4         3      21      3.                              40      53
PLOAD4         3      22      3.                              54      50
PLOAD4         3      18      3.                              81      51
PLOAD4         3      20      3.                              91      57
PLOAD4         3      22      3.                              57      50
PLOAD4         3      24      3.                              47      42
PLOAD4         3      25      3.                              67      70
PLOAD4         3      26      3.                              79      68
PLOAD4         3      27      3.                              73      69
PLOAD4         3      28      3.                              70      58
PLOAD4         3      26      3.                              68      71
PLOAD4         3      28      3.                              58      72
PLOAD4         3      30      3.                              72      85
PLOAD4         3      32      3.                              59      63
PLOAD4         3      25      3.                              67      77
PLOAD4         3      27      3.                              73      75
PLOAD4         3      29      3.                              78      65
PLOAD4         3      31      3.                              77      76
PLOAD4         3      35      3.                              85      87
PLOAD4         3      36      3.                               6      86
PLOAD4         3      39      3.                              84      52
PLOAD4         3      40      3.                              87      35
PLOAD4         3      33      3.                              89      93
PLOAD4         3      34      3.                              88      92
PLOAD4         3      37      3.                              90      46
PLOAD4         3      38      3.                              93      19
PLOAD4         3      34      3.                              15      92
PLOAD4         3      36      3.                              14      94
PLOAD4         3      38      3.                              94      19
PLOAD4         3      40      3.                              86      18
$ FEMAP Property 1 : Acoustic space
PSOLID         1       1       0        
$ FEMAP Property 2 : Rigid body
PSOLID         2       1       0        
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           1       0     1.2     -2.      0.       0        
GRID           2       0     1.2    -1.5      0.       0        
GRID           3       0     1.2     -2.    -0.5       0        
GRID           4       0     1.8     -2.     0.5       0        
GRID           5       0     1.8     -2.      0.       0        
GRID           6       0     1.8     -1.    -0.5       0        
GRID           7       0     1.8     -2.    -0.5       0        
GRID           8       0     1.2    -1.5    -0.5       0        
GRID           9       0     1.8    -1.5    -0.5       0        
GRID          10       0     1.2    -1.5     0.5       0        
GRID          11       0     1.2     -2.     0.5       0        
GRID          12       0     1.8    -1.5     0.5       0        
GRID          13       0     2.4     -2.    -0.5       0        
GRID          14       0     1.2     -1.    -0.5       0        
GRID          15       0     1.2     -1.      0.       0        
GRID          16       0     1.8     -1.      0.       0        
GRID          17       0     0.6      0.    -0.5       0        
GRID          18       0     1.2      0.      0.       0        
GRID          19       0     1.2      0.     0.5       0        
GRID          20       0     0.6      0.      0.       0        
GRID          21       0      0.      1.     0.5       0        
GRID          22       0     1.2      1.      0.       0        
GRID          23       0      0.      1.      0.       0        
GRID          24       0     0.6      1.      0.       0        
GRID          25       0      0.      0.    -0.5       0        
GRID          26       0      0.      0.      0.       0        
GRID          27       0      0.      0.     0.5       0        
GRID          28       0      0.     0.5      0.       0        
GRID          29       0     0.6      0.     0.5       0        
GRID          30       0     1.2     0.5     0.5       0        
GRID          31       0     0.6      1.     0.5       0        
GRID          32       0      0.     0.5     0.5       0        
GRID          33       0     0.6     0.5     0.5       0        
GRID          34       0     1.2      1.    -0.5       0        
GRID          35       0     1.2      0.    -0.5       0        
GRID          36       0      0.     0.5    -0.5       0        
GRID          37       0      0.      1.    -0.5       0        
GRID          38       0     0.6      1.    -0.5       0        
GRID          39       0     0.6     0.5    -0.5       0        
GRID          40       0     1.2     0.5    -0.5       0        
GRID          41       0     1.2      1.     0.5       0        
GRID          42       0     2.4      1.      0.       0        
GRID          43       0     1.8      1.      0.       0        
GRID          44       0     1.2     0.5      0.       0        
GRID          45       0     1.8      0.      0.       0        
GRID          46       0     1.8      0.     0.5       0        
GRID          47       0     2.4     0.5     0.5       0        
GRID          48       0     1.8      1.     0.5       0        
GRID          49       0     1.8     0.5     0.5       0        
GRID          50       0     2.4      1.    -0.5       0        
GRID          51       0     2.4     0.5    -0.5       0        
GRID          52       0     1.8      0.    -0.5       0        
GRID          53       0     1.8      1.    -0.5       0        
GRID          54       0     1.8     0.5    -0.5       0        
GRID          55       0     2.4      0.    -0.5       0        
GRID          56       0     2.4      1.     0.5       0        
GRID          57       0     2.4     0.5      0.       0        
GRID          58       0     3.6     -2.    -0.5       0        
GRID          59       0     3.6    -1.5    -0.5       0        
GRID          60       0     3.6     -1.      0.       0        
GRID          61       0     3.6    -1.5      0.       0        
GRID          62       0     3.6     -1.    -0.5       0        
GRID          63       0      3.     -1.    -0.5       0        
GRID          64       0     2.4     -1.      0.       0        
GRID          65       0      3.     -1.     0.5       0        
GRID          66       0      3.     -1.      0.       0        
GRID          67       0     2.4     -2.     0.5       0        
GRID          68       0      3.     -2.    -0.5       0        
GRID          69       0     3.6     -2.      0.       0        
GRID          70       0      3.     -2.      0.       0        
GRID          71       0     2.4    -1.5    -0.5       0        
GRID          72       0      3.    -1.5    -0.5       0        
GRID          73       0      3.     -2.     0.5       0        
GRID          74       0     3.6     -2.     0.5       0        
GRID          75       0     3.6    -1.5     0.5       0        
GRID          76       0     3.6     -1.     0.5       0        
GRID          77       0      3.    -1.5     0.5       0        
GRID          78       0     2.4    -1.5     0.5       0        
GRID          79       0     2.4     -2.      0.       0        
GRID          80       0     2.4    -1.5      0.       0        
GRID          81       0     2.4      0.      0.       0        
GRID          82       0     2.4    -0.5      0.       0        
GRID          83       0     1.2     -1.     0.5       0        
GRID          84       0     2.4    -0.5    -0.5       0        
GRID          85       0     2.4     -1.    -0.5       0        
GRID          86       0     1.2    -0.5    -0.5       0        
GRID          87       0     1.8    -0.5    -0.5       0        
GRID          88       0     1.8     -1.     0.5       0        
GRID          89       0     2.4     -1.     0.5       0        
GRID          90       0     2.4    -0.5     0.5       0        
GRID          91       0     2.4      0.     0.5       0        
GRID          92       0     1.2    -0.5     0.5       0        
GRID          93       0     1.8    -0.5     0.5       0        
GRID          94       0     1.2    -0.5      0.       0        
GRID          95       0     1.8    -1.5      0.       0        
GRID          96       0     0.6     0.5      0.       0        
GRID          97       0     1.8     0.5      0.       0        
GRID          98       0      3.    -1.5      0.       0        
GRID          99       0     1.8    -0.5      0.       0        
CHEXA          1       2      11       1       2      10       4       5+EL    1
+EL    1      95      12                                                        
CHEXA          2       2      10       2      15      83      12      95+EL    2
+EL    2      16      88                                                        
CHEXA          3       2       1       3       8       2       5       7+EL    3
+EL    3       9      95                                                        
CHEXA          4       2       2       8      14      15      95       9+EL    4
+EL    4       6      16                                                        
CHEXA          5       2       4       5      95      12      67      79+EL    5
+EL    5      80      78                                                        
CHEXA          6       2      12      95      16      88      78      80+EL    6
+EL    6      64      89                                                        
CHEXA          7       2       5       7       9      95      79      13+EL    7
+EL    7      71      80                                                        
CHEXA          8       2      95       9       6      16      80      71+EL    8
+EL    8      85      64                                                        
CHEXA          9       2      25      26      20      17      36      28+EL    9
+EL    9      96      39                                                        
CHEXA         10       2      17      20      18      35      39      96+EL    A
+EL    A      44      40                                                        
CHEXA         11       2      26      27      29      20      28      32+EL    B
+EL    B      33      96                                                        
CHEXA         12       2      20      29      19      18      96      33+EL    C
+EL    C      30      44                                                        
CHEXA         13       2      36      28      96      39      37      23+EL    D
+EL    D      24      38                                                        
CHEXA         14       2      39      96      44      40      38      24+EL    E
+EL    E      22      34                                                        
CHEXA         15       2      28      32      33      96      23      21+EL    F
+EL    F      31      24                                                        
CHEXA         16       2      96      33      30      44      24      31+EL    G
+EL    G      41      22                                                        
CHEXA         17       2      35      18      45      52      40      44+EL    H
+EL    H      97      54                                                        
CHEXA         18       2      52      45      81      55      54      97+EL    I
+EL    I      57      51                                                        
CHEXA         19       2      18      19      46      45      44      30+EL    J
+EL    J      49      97                                                        
CHEXA         20       2      45      46      91      81      97      49+EL    K
+EL    K      47      57                                                        
CHEXA         21       2      40      44      97      54      34      22+EL    L
+EL    L      43      53                                                        
CHEXA         22       2      54      97      57      51      53      43+EL    M
+EL    M      42      50                                                        
CHEXA         23       2      44      30      49      97      22      41+EL    N
+EL    N      48      43                                                        
CHEXA         24       2      97      49      47      57      43      48+EL    O
+EL    O      56      42                                                        
CHEXA         25       2      67      73      70      79      78      77+EL    P
+EL    P      98      80                                                        
CHEXA         26       2      79      70      68      13      80      98+EL    Q
+EL    Q      72      71                                                        
CHEXA         27       2      73      74      69      70      77      75+EL    R
+EL    R      61      98                                                        
CHEXA         28       2      70      69      58      68      98      61+EL    S
+EL    S      59      72                                                        
CHEXA         29       2      78      77      98      80      89      65+EL    T
+EL    T      66      64                                                        
CHEXA         30       2      80      98      72      71      64      66+EL    U
+EL    U      63      85                                                        
CHEXA         31       2      77      75      61      98      65      76+EL    V
+EL    V      60      66                                                        
CHEXA         32       2      98      61      59      72      66      60+EL    W
+EL    W      62      63                                                        
CHEXA         33       2      89      64      16      88      90      82+EL    X
+EL    X      99      93                                                        
CHEXA         34       2      88      16      15      83      93      99+EL    Y
+EL    Y      94      92                                                        
CHEXA         35       2      64      85       6      16      82      84+EL    Z
+EL    Z      87      99                                                        
CHEXA         36       2      16       6      14      15      99      87+EL   10
+EL   10      86      94                                                        
CHEXA         37       2      90      82      99      93      91      81+EL   11
+EL   11      45      46                                                        
CHEXA         38       2      93      99      94      92      46      45+EL   12
+EL   12      18      19                                                        
CHEXA         39       2      82      84      87      99      81      55+EL   13
+EL   13      52      45                                                        
CHEXA         40       2      99      87      86      94      45      52+EL   14
+EL   14      35      18                                                        
ENDDATA
