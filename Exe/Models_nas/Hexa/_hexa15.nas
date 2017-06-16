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
$   Date       : Mon Mar 01 17:59:06 2004
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
FORCE          1     121       0      1.   1000.      0.      0.
FORCE          1     122       0      1.   1000.      0.      0.
FORCE          1     123       0      1.   1000.      0.      0.
FORCE          1     124       0      1.   1000.      0.      0.
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
GRID          65       0      0.      8.      0.       0        
GRID          66       0      0.      8.      1.       0        
GRID          67       0      1.      8.      1.       0        
GRID          68       0      1.      8.      0.       0        
GRID          73       0      0.      9.      0.       0        
GRID          74       0      0.      9.      1.       0        
GRID          75       0      1.      9.      1.       0        
GRID          76       0      1.      9.      0.       0        
GRID          81       0      0.     10.      0.       0        
GRID          82       0      0.     10.      1.       0        
GRID          83       0      1.     10.      1.       0        
GRID          84       0      1.     10.      0.       0        
GRID          89       0      0.     11.      0.       0        
GRID          90       0      0.     11.      1.       0        
GRID          91       0      1.     11.      1.       0        
GRID          92       0      1.     11.      0.       0        
GRID          97       0      0.     12.      0.       0        
GRID          98       0      0.     12.      1.       0        
GRID          99       0      1.     12.      1.       0        
GRID         100       0      1.     12.      0.       0        
GRID         105       0      0.     13.      0.       0        
GRID         106       0      0.     13.      1.       0        
GRID         107       0      1.     13.      1.       0        
GRID         108       0      1.     13.      0.       0        
GRID         113       0      0.     14.      0.       0        
GRID         114       0      0.     14.      1.       0        
GRID         115       0      1.     14.      1.       0        
GRID         116       0      1.     14.      0.       0        
GRID         121       0      0.     15.      0.       0        
GRID         122       0      0.     15.      1.       0        
GRID         123       0      1.     15.      1.       0        
GRID         124       0      1.     15.      0.       0        
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
CHEXA          8       1      57      58      59      60      65      66+EL    8
+EL    8      67      68                                                        
CHEXA          9       1      65      66      67      68      73      74+EL    9
+EL    9      75      76                                                        
CHEXA         10       1      73      74      75      76      81      82+EL    A
+EL    A      83      84                                                        
CHEXA         11       1      81      82      83      84      89      90+EL    B
+EL    B      91      92                                                        
CHEXA         12       1      89      90      91      92      97      98+EL    C
+EL    C      99     100                                                        
CHEXA         13       1      97      98      99     100     105     106+EL    D
+EL    D     107     108                                                        
CHEXA         14       1     105     106     107     108     113     114+EL    E
+EL    E     115     116                                                        
CHEXA         15       1     113     114     115     116     121     122+EL    F
+EL    F     123     124                                                        
ENDDATA
