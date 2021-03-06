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
$   Date       : Fri Aug 04 09:07:15 2006
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
PBEAM          1       1  0.0019  1.8E-6  1.8E-6-1.07E-66.239E-8      0.+PR    1
+PR    1 -5.3E-3 -5.3E-30.094699 -5.3E-30.0046990.094699 -5.3E-30.094699+PA    1
+PA    1    YESA      1.                                                +PC    1
+PC    1  0.4313  0.4313                                                +PD    1
+PD    1                                0.0233830.0233830.0233830.023383        
$ FEMAP Property 2 : Untitled
PBEAM          2       1 3.14E-47.854E-97.854E-9      0.1.569E-8      0.+PR    2
+PR    2      0.   -0.01    0.01      0.      0.    0.01   -0.01      0.+PA    2
+PA    2    YESA      1.                                                +PC    2
+PC    2 0.88652 0.88653                                                        
$ FEMAP Material 1 : Untitled
MAT1           1 2.1E+11 1.8E+11     0.3   2700.      0.      0.        
$ FEMAP Material 2 : Untitled
MAT1           2 7.2E+10  6.E+10     0.3   2700.      0.      0.        
GRID           2       0   -0.49    0.25      0.       0        
GRID           3       0   -0.48     0.5      0.       0        
GRID           4       0   -0.47    0.75      0.       0        
GRID           5       0   -0.46      1.      0.       0        
GRID           7       0   -0.44     1.5      0.       0        
GRID           8       0   -0.43    1.75      0.       0        
GRID           9       0   -0.42      2.      0.       0        
GRID          10       0   -0.41    2.25      0.       0        
GRID          11       0    -0.4     2.5      0.       0        
GRID          12       0   -0.39    2.75      0.       0        
GRID          13       0   -0.38      3.      0.       0        
GRID          14       0   -0.37    3.25      0.       0        
GRID          15       0   -0.36     3.5      0.       0        
GRID          17       0   -0.34      4.      0.       0        
GRID          18       0   -0.33    4.25      0.       0        
GRID          19       0   -0.32     4.5      0.       0        
GRID          20       0   -0.31    4.75      0.       0        
GRID          22       0   -0.29    5.25      0.       0        
GRID          23       0   -0.28     5.5      0.       0        
GRID          24       0   -0.27    5.75      0.       0        
GRID          25       0   -0.26      6.      0.       0        
GRID          26       0   -0.25    6.25      0.       0        
GRID          27       0   -0.24     6.5      0.       0        
GRID          28       0   -0.23    6.75      0.       0        
GRID          29       0   -0.22      7.      0.       0        
GRID          30       0   -0.21    7.25      0.       0        
GRID          31       0    -0.2     7.5      0.       0        
GRID          32       0   -0.19    7.75      0.       0        
GRID          33       0   -0.18      8.      0.       0        
GRID          34       0   -0.17    8.25      0.       0        
GRID          35       0   -0.16     8.5      0.       0        
GRID          36       0   -0.15    8.75      0.       0        
GRID          37       0   -0.14      9.      0.       0        
GRID          38       0   -0.13    9.25      0.       0        
GRID          39       0   -0.12     9.5      0.       0        
GRID          40       0   -0.11    9.75      0.       0        
GRID          42       0     0.5      0.      0.       0        
GRID          43       0    0.49    0.25      0.       0        
GRID          44       0    0.48     0.5      0.       0        
GRID          45       0    0.47    0.75      0.       0        
GRID          46       0    0.46      1.      0.       0        
GRID          48       0    0.44     1.5      0.       0        
GRID          49       0    0.43    1.75      0.       0        
GRID          50       0    0.42      2.      0.       0        
GRID          51       0    0.41    2.25      0.       0        
GRID          52       0     0.4     2.5      0.       0        
GRID          53       0    0.39    2.75      0.       0        
GRID          54       0    0.38      3.      0.       0        
GRID          55       0    0.37    3.25      0.       0        
GRID          56       0    0.36     3.5      0.       0        
GRID          57       0    0.35    3.75      0.       0        
GRID          58       0    0.34      4.      0.       0        
GRID          59       0    0.33    4.25      0.       0        
GRID          60       0    0.32     4.5      0.       0        
GRID          61       0    0.31    4.75      0.       0        
GRID          63       0    0.29    5.25      0.       0        
GRID          64       0    0.28     5.5      0.       0        
GRID          65       0    0.27    5.75      0.       0        
GRID          66       0    0.26      6.      0.       0        
GRID          68       0    0.24     6.5      0.       0        
GRID          69       0    0.23    6.75      0.       0        
GRID          70       0    0.22      7.      0.       0        
GRID          71       0    0.21    7.25      0.       0        
GRID          73       0    0.19    7.75      0.       0        
GRID          74       0    0.18      8.      0.       0        
GRID          75       0    0.17    8.25      0.       0        
GRID          76       0    0.16     8.5      0.       0        
GRID          77       0    0.15    8.75      0.       0        
GRID          78       0    0.14      9.      0.       0        
GRID          79       0    0.13    9.25      0.       0        
GRID          80       0    0.12     9.5      0.       0        
GRID          81       0    0.11    9.75      0.       0        
GRID          83       0    -0.1     10.      0.       0        
GRID          84       0     0.1     10.      0.       0        
GRID          85       0    0.45    1.25      0.       0        
GRID          86       0    0.27    1.25      0.       0        
GRID          87       0    0.09    1.25      0.       0        
GRID          88       0   -0.09    1.25      0.       0        
GRID          89       0   -0.27    1.25      0.       0        
GRID          90       0   -0.45    1.25      0.       0        
GRID          92       0     0.2     2.5      0.       0        
GRID          93       0      0.     2.5      0.       0        
GRID          94       0    -0.2     2.5      0.       0        
GRID          97       0   0.175    3.75      0.       0        
GRID          98       0      0.    3.75      0.       0        
GRID          99       0  -0.175    3.75      0.       0        
GRID         100       0   -0.35    3.75      0.       0        
GRID         102       0     0.1      5.      0.       0        
GRID         103       0    -0.1      5.      0.       0        
GRID         105       0    0.25    6.25      0.       0        
GRID         106       00.083333    6.25      0.       0        
GRID         107       0 -0.0833    6.25      0.       0        
GRID         110       00.066667     7.5      0.       0        
GRID         111       0 -0.0667     7.5      0.       0        
GRID         114       0      0.    8.75      0.       0        
GRID         117       00.066667 9.83333      0.       0        
GRID         118       00.033333 9.66667      0.       0        
GRID         121       0   -0.05    9.25      0.       0        
GRID         122       0    -0.1      9.      0.       0        
GRID         125       0     0.1 8.57143      0.       0        
GRID         126       0    0.05 8.39286      0.       0        
GRID         129       0 -0.0667 7.97619      0.       0        
GRID         130       0-0.13333  7.7381      0.       0        
GRID         132       0     0.2     7.5      0.       0        
GRID         133       0 0.13333 7.31481      0.       0        
GRID         134       00.066667 7.12963      0.       0        
GRID         137       0 -0.0833 6.71296      0.       0        
GRID         138       0-0.16667 6.48148      0.       0        
GRID         141       0  0.1875 6.10795      0.       0        
GRID         142       0   0.125 5.96591      0.       0        
GRID         143       0  0.0625 5.82386      0.       0        
GRID         146       0  -0.075 5.51136      0.       0        
GRID         147       0   -0.15 5.34091      0.       0        
GRID         148       0  -0.225 5.17045      0.       0        
GRID         149       0    -0.3      5.      0.       0        
GRID         150       0     0.3      5.      0.       0        
GRID         151       0   0.225 4.85577      0.       0        
GRID         152       0    0.15 4.71154      0.       0        
GRID         153       0   0.075 4.56731      0.       0        
GRID         155       0      0. 4.42308      0.       0        
GRID         156       0 -0.0875 4.25481      0.       0        
GRID         157       0  -0.175 4.08654      0.       0        
GRID         158       0 -0.2625 3.91827      0.       0        
GRID         161       0  0.2625 3.60417      0.       0        
GRID         162       0   0.175 3.45833      0.       0        
GRID         163       0  0.0875  3.3125      0.       0        
GRID         164       0      0. 3.16667      0.       0        
GRID         166       0    -0.1      3.      0.       0        
GRID         167       0    -0.2 2.83333      0.       0        
GRID         168       0    -0.3 2.66667      0.       0        
GRID         171       0     0.3 2.35294      0.       0        
GRID         172       0     0.2 2.20588      0.       0        
GRID         173       0     0.1 2.05882      0.       0        
GRID         174       0      0. 1.91176      0.       0        
GRID         176       0 -0.1125 1.74632      0.       0        
GRID         177       0  -0.225 1.58088      0.       0        
GRID         178       0 -0.3375 1.41544      0.       0        
GRID         181       0  0.3375 1.10197      0.       0        
GRID         182       0   0.225 0.95395      0.       0        
GRID         183       0  0.1125 0.80592      0.       0        
GRID         184       0      0. 0.65789      0.       0        
GRID         186       0  -0.125 0.49342      0.       0        
GRID         187       0   -0.25 0.32895      0.       0        
GRID         188       0  -0.375 0.16447      0.       0        
GRID         189       0    -0.5      0.      0.       0        
GRID         191       0 -0.0667 9.83333      0.       0        
GRID         192       0 -0.0333 9.66667      0.       0        
GRID         194       0      0.     9.5      0.       0        
GRID         195       0    0.05    9.25      0.       0        
GRID         196       0     0.1      9.      0.       0        
GRID         199       0 -0.3375 1.10197      0.       0        
GRID         200       0  -0.225 0.95395      0.       0        
GRID         201       0 -0.1125 0.80592      0.       0        
GRID         204       0   0.125 0.49342      0.       0        
GRID         205       0    0.25 0.32895      0.       0        
GRID         206       0   0.375 0.16447      0.       0        
GRID         209       0  -0.225 4.85577      0.       0        
GRID         210       0   -0.15 4.71154      0.       0        
GRID         211       0  -0.075 4.56731      0.       0        
GRID         214       0  0.0875 4.25481      0.       0        
GRID         215       0   0.175 4.08654      0.       0        
GRID         216       0  0.2625 3.91827      0.       0        
GRID         219       0 -0.2625 3.60417      0.       0        
GRID         220       0  -0.175 3.45833      0.       0        
GRID         221       0 -0.0875  3.3125      0.       0        
GRID         224       0     0.1      3.      0.       0        
GRID         225       0     0.2 2.83333      0.       0        
GRID         226       0     0.3 2.66667      0.       0        
GRID         229       0    -0.3 2.35294      0.       0        
GRID         230       0    -0.2 2.20588      0.       0        
GRID         231       0    -0.1 2.05882      0.       0        
GRID         234       0  0.1125 1.74632      0.       0        
GRID         235       0   0.225 1.58088      0.       0        
GRID         236       0  0.3375 1.41544      0.       0        
GRID         239       0 -0.1875 6.10795      0.       0        
GRID         240       0  -0.125 5.96591      0.       0        
GRID         241       0 -0.0625 5.82386      0.       0        
GRID         243       0      0. 5.68182      0.       0        
GRID         244       0   0.075 5.51136      0.       0        
GRID         245       0    0.15 5.34091      0.       0        
GRID         246       0   0.225 5.17045      0.       0        
GRID         249       0-0.13333 7.31481      0.       0        
GRID         250       0 -0.0667 7.12963      0.       0        
GRID         252       0      0. 6.94444      0.       0        
GRID         253       00.083333 6.71296      0.       0        
GRID         254       0 0.16667 6.48148      0.       0        
GRID         257       0    -0.1 8.57143      0.       0        
GRID         258       0   -0.05 8.39286      0.       0        
GRID         260       0      0. 8.21429      0.       0        
GRID         261       00.066667 7.97619      0.       0        
GRID         262       0 0.13333  7.7381      0.       0        
CBEAM          1       1     189       2-0.23369 0.97231      0.
CBEAM          2       1       2       3-0.23369 0.97231      0.
CBEAM          3       1       3       4-0.23369 0.97231      0.
CBEAM          4       1       4       5-0.23369 0.97231      0.
CBEAM          5       1       5      90-0.23369 0.97231      0.
CBEAM          6       1      90       7-0.23369 0.97231      0.
CBEAM          7       1       7       8-0.23369 0.97231      0.
CBEAM          8       1       8       9-0.23369 0.97231      0.
CBEAM          9       1       9      10-0.23369 0.97231      0.
CBEAM         10       1      10      11-0.23369 0.97231      0.
CBEAM         11       1      11      12-0.23369 0.97231      0.
CBEAM         12       1      12      13-0.23369 0.97231      0.
CBEAM         13       1      13      14-0.23369 0.97231      0.
CBEAM         14       1      14      15-0.23369 0.97231      0.
CBEAM         15       1      15     100-0.23369 0.97231      0.
CBEAM         16       1     100      17-0.23369 0.97231      0.
CBEAM         17       1      17      18-0.23369 0.97231      0.
CBEAM         18       1      18      19-0.23369 0.97231      0.
CBEAM         19       1      19      20-0.23369 0.97231      0.
CBEAM         20       1      20     149-0.23369 0.97231      0.
CBEAM         21       1     149      22-0.23369 0.97231      0.
CBEAM         22       1      22      23-0.23369 0.97231      0.
CBEAM         23       1      23      24-0.23369 0.97231      0.
CBEAM         24       1      24      25-0.23369 0.97231      0.
CBEAM         25       1      25      26-0.23369 0.97231      0.
CBEAM         26       1      26      27-0.23369 0.97231      0.
CBEAM         27       1      27      28-0.23369 0.97231      0.
CBEAM         28       1      28      29-0.23369 0.97231      0.
CBEAM         29       1      29      30-0.23369 0.97231      0.
CBEAM         30       1      30      31-0.23369 0.97231      0.
CBEAM         31       1      31      32-0.23369 0.97231      0.
CBEAM         32       1      32      33-0.23369 0.97231      0.
CBEAM         33       1      33      34-0.23369 0.97231      0.
CBEAM         34       1      34      35-0.23369 0.97231      0.
CBEAM         35       1      35      36-0.23369 0.97231      0.
CBEAM         36       1      36      37-0.23369 0.97231      0.
CBEAM         37       1      37      38-0.23369 0.97231      0.
CBEAM         38       1      38      39-0.23369 0.97231      0.
CBEAM         39       1      39      40-0.23369 0.97231      0.
CBEAM         40       1      40      83-0.23369 0.97231      0.
CBEAM         41       1      42      43-0.23369 0.97231      0.
CBEAM         42       1      43      44-0.23369 0.97231      0.
CBEAM         43       1      44      45-0.23369 0.97231      0.
CBEAM         44       1      45      46-0.23369 0.97231      0.
CBEAM         45       1      46      85-0.23369 0.97231      0.
CBEAM         46       1      85      48-0.23369 0.97231      0.
CBEAM         47       1      48      49-0.23369 0.97231      0.
CBEAM         48       1      49      50-0.23369 0.97231      0.
CBEAM         49       1      50      51-0.23369 0.97231      0.
CBEAM         50       1      51      52-0.23369 0.97231      0.
CBEAM         51       1      52      53-0.23369 0.97231      0.
CBEAM         52       1      53      54-0.23369 0.97231      0.
CBEAM         53       1      54      55-0.23369 0.97231      0.
CBEAM         54       1      55      56-0.23369 0.97231      0.
CBEAM         55       1      56      57-0.23369 0.97231      0.
CBEAM         56       1      57      58-0.23369 0.97231      0.
CBEAM         57       1      58      59-0.23369 0.97231      0.
CBEAM         58       1      59      60-0.23369 0.97231      0.
CBEAM         59       1      60      61-0.23369 0.97231      0.
CBEAM         60       1      61     150-0.23369 0.97231      0.
CBEAM         61       1     150      63-0.23369 0.97231      0.
CBEAM         62       1      63      64-0.23369 0.97231      0.
CBEAM         63       1      64      65-0.23369 0.97231      0.
CBEAM         64       1      65      66-0.23369 0.97231      0.
CBEAM         65       1      66     105-0.23369 0.97231      0.
CBEAM         66       1     105      68-0.23369 0.97231      0.
CBEAM         67       1      68      69-0.23369 0.97231      0.
CBEAM         68       1      69      70-0.23369 0.97231      0.
CBEAM         69       1      70      71-0.23369 0.97231      0.
CBEAM         70       1      71     132-0.23369 0.97231      0.
CBEAM         71       1     132      73-0.23369 0.97231      0.
CBEAM         72       1      73      74-0.23369 0.97231      0.
CBEAM         73       1      74      75-0.23369 0.97231      0.
CBEAM         74       1      75      76-0.23369 0.97231      0.
CBEAM         75       1      76      77-0.23369 0.97231      0.
CBEAM         76       1      77      78-0.23369 0.97231      0.
CBEAM         77       1      78      79-0.23369 0.97231      0.
CBEAM         78       1      79      80-0.23369 0.97231      0.
CBEAM         79       1      80      81-0.23369 0.97231      0.
CBEAM         80       1      81      84-0.23369 0.97231      0.
CBEAM         81       1      83      84-0.23369 0.97231      0.
CBEAM         82       1      85      86-0.23369 0.97231      0.
CBEAM         83       1      86      87-0.23369 0.97231      0.
CBEAM         84       1      87      88-0.23369 0.97231      0.
CBEAM         85       1      88      89-0.23369 0.97231      0.
CBEAM         86       1      89      90-0.23369 0.97231      0.
CBEAM         87       1      52      92-0.23369 0.97231      0.
CBEAM         88       1      92      93-0.23369 0.97231      0.
CBEAM         89       1      93      94-0.23369 0.97231      0.
CBEAM         90       1      94      11-0.23369 0.97231      0.
CBEAM         91       1      57      97-0.23369 0.97231      0.
CBEAM         92       1      97      98-0.23369 0.97231      0.
CBEAM         93       1      98      99-0.23369 0.97231      0.
CBEAM         94       1      99     100-0.23369 0.97231      0.
CBEAM         95       1     150     102-0.23369 0.97231      0.
CBEAM         96       1     102     103-0.23369 0.97231      0.
CBEAM         97       1     103     149-0.23369 0.97231      0.
CBEAM         98       1     105     106-0.23369 0.97231      0.
CBEAM         99       1     106     107-0.23369 0.97231      0.
CBEAM        100       1     107      26-0.23369 0.97231      0.
CBEAM        101       1     132     110-0.23369 0.97231      0.
CBEAM        102       1     110     111-0.23369 0.97231      0.
CBEAM        103       1     111      31-0.23369 0.97231      0.
CBEAM        104       1      77     114-0.23369 0.97231      0.
CBEAM        105       1     114      36-0.23369 0.97231      0.
CBEAM        106       1      84     117-0.23369 0.97231      0.
CBEAM        107       1     117     118-0.23369 0.97231      0.
CBEAM        108       1     118     194-0.23369 0.97231      0.
CBEAM        109       1     194     121-0.23369 0.97231      0.
CBEAM        110       1     121     122-0.23369 0.97231      0.
CBEAM        111       1     122      36-0.23369 0.97231      0.
CBEAM        112       1      77     125-0.23369 0.97231      0.
CBEAM        113       1     125     126-0.23369 0.97231      0.
CBEAM        114       1     126     260-0.23369 0.97231      0.
CBEAM        115       1     260     129-0.23369 0.97231      0.
CBEAM        116       1     129     130-0.23369 0.97231      0.
CBEAM        117       1     130      31-0.23369 0.97231      0.
CBEAM        118       1     132     133-0.23369 0.97231      0.
CBEAM        119       1     133     134-0.23369 0.97231      0.
CBEAM        120       1     134     252-0.23369 0.97231      0.
CBEAM        121       1     252     137-0.23369 0.97231      0.
CBEAM        122       1     137     138-0.23369 0.97231      0.
CBEAM        123       1     138      26-0.23369 0.97231      0.
CBEAM        124       1     105     141-0.23369 0.97231      0.
CBEAM        125       1     141     142-0.23369 0.97231      0.
CBEAM        126       1     142     143-0.23369 0.97231      0.
CBEAM        127       1     143     243-0.23369 0.97231      0.
CBEAM        128       1     243     146-0.23369 0.97231      0.
CBEAM        129       1     146     147-0.23369 0.97231      0.
CBEAM        130       1     147     148-0.23369 0.97231      0.
CBEAM        131       1     148     149-0.23369 0.97231      0.
CBEAM        132       1     150     151-0.23369 0.97231      0.
CBEAM        133       1     151     152-0.23369 0.97231      0.
CBEAM        134       1     152     153-0.23369 0.97231      0.
CBEAM        135       1     153     155-0.23369 0.97231      0.
CBEAM        136       1     155     156-0.23369 0.97231      0.
CBEAM        137       1     156     157-0.23369 0.97231      0.
CBEAM        138       1     157     158-0.23369 0.97231      0.
CBEAM        139       1     158     100-0.23369 0.97231      0.
CBEAM        140       1      57     161-0.23369 0.97231      0.
CBEAM        141       1     161     162-0.23369 0.97231      0.
CBEAM        142       1     162     163-0.23369 0.97231      0.
CBEAM        143       1     163     164-0.23369 0.97231      0.
CBEAM        144       1     164     166-0.23369 0.97231      0.
CBEAM        145       1     166     167-0.23369 0.97231      0.
CBEAM        146       1     167     168-0.23369 0.97231      0.
CBEAM        147       1     168      11-0.23369 0.97231      0.
CBEAM        148       1      52     171-0.23369 0.97231      0.
CBEAM        149       1     171     172-0.23369 0.97231      0.
CBEAM        150       1     172     173-0.23369 0.97231      0.
CBEAM        151       1     173     174-0.23369 0.97231      0.
CBEAM        152       1     174     176-0.23369 0.97231      0.
CBEAM        153       1     176     177-0.23369 0.97231      0.
CBEAM        154       1     177     178-0.23369 0.97231      0.
CBEAM        155       1     178      90-0.23369 0.97231      0.
CBEAM        156       1      85     181-0.23369 0.97231      0.
CBEAM        157       1     181     182-0.23369 0.97231      0.
CBEAM        158       1     182     183-0.23369 0.97231      0.
CBEAM        159       1     183     184-0.23369 0.97231      0.
CBEAM        160       1     184     186-0.23369 0.97231      0.
CBEAM        161       1     186     187-0.23369 0.97231      0.
CBEAM        162       1     187     188-0.23369 0.97231      0.
CBEAM        163       1     188     189-0.23369 0.97231      0.
CBEAM        164       1      83     191-0.23369 0.97231      0.
CBEAM        165       1     191     192-0.23369 0.97231      0.
CBEAM        166       1     192     194-0.23369 0.97231      0.
CBEAM        167       1     194     195-0.23369 0.97231      0.
CBEAM        168       1     195     196-0.23369 0.97231      0.
CBEAM        169       1     196      77-0.23369 0.97231      0.
CBEAM        170       1      90     199-0.23369 0.97231      0.
CBEAM        171       1     199     200-0.23369 0.97231      0.
CBEAM        172       1     200     201-0.23369 0.97231      0.
CBEAM        173       1     201     184-0.23369 0.97231      0.
CBEAM        174       1     184     204-0.23369 0.97231      0.
CBEAM        175       1     204     205-0.23369 0.97231      0.
CBEAM        176       1     205     206-0.23369 0.97231      0.
CBEAM        177       1     206      42-0.23369 0.97231      0.
CBEAM        178       1     149     209-0.23369 0.97231      0.
CBEAM        179       1     209     210-0.23369 0.97231      0.
CBEAM        180       1     210     211-0.23369 0.97231      0.
CBEAM        181       1     211     155-0.23369 0.97231      0.
CBEAM        182       1     155     214-0.23369 0.97231      0.
CBEAM        183       1     214     215-0.23369 0.97231      0.
CBEAM        184       1     215     216-0.23369 0.97231      0.
CBEAM        185       1     216      57-0.23369 0.97231      0.
CBEAM        186       1     100     219-0.23369 0.97231      0.
CBEAM        187       1     219     220-0.23369 0.97231      0.
CBEAM        188       1     220     221-0.23369 0.97231      0.
CBEAM        189       1     221     164-0.23369 0.97231      0.
CBEAM        190       1     164     224-0.23369 0.97231      0.
CBEAM        191       1     224     225-0.23369 0.97231      0.
CBEAM        192       1     225     226-0.23369 0.97231      0.
CBEAM        193       1     226      52-0.23369 0.97231      0.
CBEAM        194       1      11     229-0.23369 0.97231      0.
CBEAM        195       1     229     230-0.23369 0.97231      0.
CBEAM        196       1     230     231-0.23369 0.97231      0.
CBEAM        197       1     231     174-0.23369 0.97231      0.
CBEAM        198       1     174     234-0.23369 0.97231      0.
CBEAM        199       1     234     235-0.23369 0.97231      0.
CBEAM        200       1     235     236-0.23369 0.97231      0.
CBEAM        201       1     236      85-0.23369 0.97231      0.
CBEAM        202       1      26     239-0.23369 0.97231      0.
CBEAM        203       1     239     240-0.23369 0.97231      0.
CBEAM        204       1     240     241-0.23369 0.97231      0.
CBEAM        205       1     241     243-0.23369 0.97231      0.
CBEAM        206       1     243     244-0.23369 0.97231      0.
CBEAM        207       1     244     245-0.23369 0.97231      0.
CBEAM        208       1     245     246-0.23369 0.97231      0.
CBEAM        209       1     246     150-0.23369 0.97231      0.
CBEAM        210       1      31     249-0.23369 0.97231      0.
CBEAM        211       1     249     250-0.23369 0.97231      0.
CBEAM        212       1     250     252-0.23369 0.97231      0.
CBEAM        213       1     252     253-0.23369 0.97231      0.
CBEAM        214       1     253     254-0.23369 0.97231      0.
CBEAM        215       1     254     105-0.23369 0.97231      0.
CBEAM        216       1      36     257-0.23369 0.97231      0.
CBEAM        217       1     257     258-0.23369 0.97231      0.
CBEAM        218       1     258     260-0.23369 0.97231      0.
CBEAM        219       1     260     261-0.23369 0.97231      0.
CBEAM        220       1     261     262-0.23369 0.97231      0.
CBEAM        221       1     262     132-0.23369 0.97231      0.
ENDDATA
