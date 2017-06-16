ID C:\Users\Dima\Lame3d2\Mpl,FEMAP
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
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\Models_nas\Box3D\Solid\Box_model.MOD
$   Date       : Sun Feb 08 11:53:33 2009
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
$ FEMAP Load Set 1 : Source Vn -.001 0.0
PLOAD4         1      31      1.                              46      51
PLOAD4         1      32      1.                              66      52
PLOAD4         1      33      1.                              67      58
PLOAD4         1      34      1.                              47      73
PLOAD4         1      35      1.                              51      49
PLOAD4         1      36      1.                              52      72
$ FEMAP Load Set 2 : Boundary Skews
PLOAD4         2      39      2.                              86      64
PLOAD4         2      42      2.                              91      93
PLOAD4         2      37      2.                              27      56
PLOAD4         2      40      2.                              95      84
PLOAD4         2      37      2.                              27      65
PLOAD4         2      38      2.                              44      53
PLOAD4         2      39      2.                              96      64
PLOAD4         2      40      2.                              89      55
PLOAD4         2      41      2.                              88      71
PLOAD4         2      42      2.                              87      54
$ FEMAP Property 1 : Acoustic
PSOLID         1       1       0        
$ FEMAP Property 2 : Rigid body
PSOLID         2       1       0        
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           1       0    0.75      0.    2.05       0        
GRID           2       0    0.25      0.      3.       0        
GRID           3       0   -0.25      0.   2.525       0        
GRID           4       0   -0.25      0.    2.05       0        
GRID           5       0   -0.25      0.   1.575       0        
GRID           6       0    0.25      0.   2.525       0        
GRID           7       0    0.25      0.    2.05       0        
GRID           8       0    0.25      0.   1.575       0        
GRID           9       0    0.75      1.    2.05       0        
GRID          10       0   -0.75      1.   2.525       0        
GRID          11       0   -0.25      1.      3.       0        
GRID          12       0    0.25      1.   2.525       0        
GRID          13       0    0.25      1.    2.05       0        
GRID          14       0    0.25      1.   1.575       0        
GRID          15       0   -0.25      1.   2.525       0        
GRID          16       0   -0.25      1.    2.05       0        
GRID          17       0   -0.25      1.   1.575       0        
GRID          18       0   -0.75      1.    2.05       0        
GRID          19       0   -0.75      1.   1.575       0        
GRID          20       0   -0.75      0.   1.575       0        
GRID          21       0   -0.75      0.    2.05       0        
GRID          22       0   -0.75      0.   2.525       0        
GRID          23       0   -0.75     0.5      3.       0        
GRID          24       0   -0.75     0.5   2.525       0        
GRID          25       0   -0.75     0.5    2.05       0        
GRID          26       0   -0.75     0.5   1.575       0        
GRID          27       0    0.75      1.     1.1       0        
GRID          28       0    0.75      1.   1.575       0        
GRID          29       0    0.75      1.   2.525       0        
GRID          30       0    0.75      1.      3.       0        
GRID          31       0    0.75      0.   2.525       0        
GRID          32       0    0.75      0.   1.575       0        
GRID          33       0    0.75     0.5   1.575       0        
GRID          34       0    0.75     0.5    2.05       0        
GRID          35       0    0.75     0.5   2.525       0        
GRID          36       0    0.75      0.      3.       0        
GRID          37       0    0.75     0.5      3.       0        
GRID          38       0    0.25      1.      3.       0        
GRID          39       0   -0.75      1.      3.       0        
GRID          40       0   -0.75      0.      3.       0        
GRID          41       0   -0.25      0.      3.       0        
GRID          42       0    0.25     0.5      3.       0        
GRID          43       0   -0.25     0.5      3.       0        
GRID          44       0    0.25      1.     1.1       0        
GRID          45       0    0.25     0.5     1.1       0        
GRID          46       0    0.75      1.      0.       0        
GRID          47       0    0.75     0.5      0.       0        
GRID          48       0    0.75      0.      0.       0        
GRID          49       0   -0.25      0.      0.       0        
GRID          50       0   -0.75      1.      0.       0        
GRID          51       0    0.25     0.5      0.       0        
GRID          52       0   -0.25     0.5      0.       0        
GRID          53       0   -0.25      1.     0.8       0        
GRID          54       0   -0.75      0.     0.8       0        
GRID          55       0    0.25      0.     0.8       0        
GRID          56       0    0.75     0.5     0.8       0        
GRID          57       0   -0.25     0.5     0.8       0        
GRID          58       0   -0.75     0.5      0.       0        
GRID          59       0   -0.75      0.     0.4       0        
GRID          60       0   -0.75     0.5     0.4       0        
GRID          61       0    0.75      1.     0.4       0        
GRID          62       0    0.75      0.     0.4       0        
GRID          63       0    0.75     0.5     0.4       0        
GRID          64       0   -0.75      1.     0.8       0        
GRID          65       0    0.25      1.     0.8       0        
GRID          66       0    0.25      1.      0.       0        
GRID          67       0   -0.25      1.      0.       0        
GRID          68       0   -0.75      1.     0.4       0        
GRID          69       0   -0.25      1.     0.4       0        
GRID          70       0    0.25      1.     0.4       0        
GRID          71       0   -0.25      0.     0.8       0        
GRID          72       0   -0.75      0.      0.       0        
GRID          73       0    0.25      0.      0.       0        
GRID          74       0    0.25      0.     0.4       0        
GRID          75       0   -0.25      0.     0.4       0        
GRID          76       0   -0.25     0.5   2.525       0        
GRID          77       0   -0.25     0.5    2.05       0        
GRID          78       0   -0.25     0.5   1.575       0        
GRID          79       0    0.25     0.5   2.525       0        
GRID          80       0    0.25     0.5    2.05       0        
GRID          81       0    0.25     0.5   1.575       0        
GRID          82       0    0.25     0.5     0.4       0        
GRID          83       0   -0.25     0.5     0.4       0        
GRID          84       0    0.75      0.     0.8       0        
GRID          85       0    0.25     0.5     0.8       0        
GRID          86       0   -0.75     0.5     1.1       0        
GRID          87       0   -0.25      0.     1.1       0        
GRID          88       0    0.25      0.     1.1       0        
GRID          89       0    0.75      0.     1.1       0        
GRID          90       0   -0.25     0.5     1.1       0        
GRID          91       0   -0.75      0.     1.1       0        
GRID          92       0   -0.75      1.     1.1       0        
GRID          93       0   -0.75     0.5     0.8       0        
GRID          94       0    0.75      1.     0.8       0        
GRID          95       0    0.75     0.5     1.1       0        
GRID          96       0   -0.25      1.     1.1       0        
CHEXA          1       1      40      41       3      22      23      43+EL    1
+EL    1      76      24                                                        
CHEXA          2       1      22       3       4      21      24      76+EL    2
+EL    2      77      25                                                        
CHEXA          3       1      21       4       5      20      25      77+EL    3
+EL    3      78      26                                                        
CHEXA          4       1      20       5      87      91      26      78+EL    4
+EL    4      90      86                                                        
CHEXA          5       1      41       2       6       3      43      42+EL    5
+EL    5      79      76                                                        
CHEXA          6       1       3       6       7       4      76      79+EL    6
+EL    6      80      77                                                        
CHEXA          7       1       4       7       8       5      77      80+EL    7
+EL    7      81      78                                                        
CHEXA          8       1       5       8      88      87      78      81+EL    8
+EL    8      45      90                                                        
CHEXA          9       1       2      36      31       6      42      37+EL    9
+EL    9      35      79                                                        
CHEXA         10       1       6      31       1       7      79      35+EL    A
+EL    A      34      80                                                        
CHEXA         11       1       7       1      32       8      80      34+EL    B
+EL    B      33      81                                                        
CHEXA         12       1       8      32      89      88      81      33+EL    C
+EL    C      95      45                                                        
CHEXA         13       1      23      43      76      24      39      11+EL    D
+EL    D      15      10                                                        
CHEXA         14       1      24      76      77      25      10      15+EL    E
+EL    E      16      18                                                        
CHEXA         15       1      25      77      78      26      18      16+EL    F
+EL    F      17      19                                                        
CHEXA         16       1      26      78      90      86      19      17+EL    G
+EL    G      96      92                                                        
CHEXA         17       1      43      42      79      76      11      38+EL    H
+EL    H      12      15                                                        
CHEXA         18       1      76      79      80      77      15      12+EL    I
+EL    I      13      16                                                        
CHEXA         19       1      77      80      81      78      16      13+EL    J
+EL    J      14      17                                                        
CHEXA         20       1      78      81      45      90      17      14+EL    K
+EL    K      44      96                                                        
CHEXA         21       1      42      37      35      79      38      30+EL    L
+EL    L      29      12                                                        
CHEXA         22       1      79      35      34      80      12      29+EL    M
+EL    M       9      13                                                        
CHEXA         23       1      80      34      33      81      13       9+EL    N
+EL    N      28      14                                                        
CHEXA         24       1      81      33      95      45      14      28+EL    O
+EL    O      27      44                                                        
CHEXA         25       1      94      56      85      65      61      63+EL    P
+EL    P      82      70                                                        
CHEXA         26       1      65      85      57      53      70      82+EL    Q
+EL    Q      83      69                                                        
CHEXA         27       1      53      57      93      64      69      83+EL    R
+EL    R      60      68                                                        
CHEXA         28       1      56      84      55      85      63      62+EL    S
+EL    S      74      82                                                        
CHEXA         29       1      85      55      71      57      82      74+EL    T
+EL    T      75      83                                                        
CHEXA         30       1      57      71      54      93      83      75+EL    U
+EL    U      59      60                                                        
CHEXA         31       1      61      63      82      70      46      47+EL    V
+EL    V      51      66                                                        
CHEXA         32       1      70      82      83      69      66      51+EL    W
+EL    W      52      67                                                        
CHEXA         33       1      69      83      60      68      67      52+EL    X
+EL    X      58      50                                                        
CHEXA         34       1      63      62      74      82      47      48+EL    Y
+EL    Y      73      51                                                        
CHEXA         35       1      82      74      75      83      51      73+EL    Z
+EL    Z      49      52                                                        
CHEXA         36       1      83      75      59      60      52      49+EL   10
+EL   10      72      58                                                        
CHEXA         37       2      27      95      45      44      94      56+EL   11
+EL   11      85      65                                                        
CHEXA         38       2      44      45      90      96      65      85+EL   12
+EL   12      57      53                                                        
CHEXA         39       2      96      90      86      92      53      57+EL   13
+EL   13      93      64                                                        
CHEXA         40       2      95      89      88      45      56      84+EL   14
+EL   14      55      85                                                        
CHEXA         41       2      45      88      87      90      85      55+EL   15
+EL   15      71      57                                                        
CHEXA         42       2      90      87      91      86      57      71+EL   16
+EL   16      54      93                                                        
ENDDATA
