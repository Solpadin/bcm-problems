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
$   Date       : Fri Feb 27 18:47:58 2004
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
$ FEMAP Load Set 1 : NASTRAN 1
FORCE          1      61       0      1.   1000.      0.      0.
FORCE          1      62       0      1.   1000.      0.      0.
FORCE          1      63       0      1.   1000.      0.      0.
FORCE          1      64       0      1.   1000.      0.      0.
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
GRID          37       0      0.      6.      0.       0        
GRID          38       0      0.      6.      1.       0        
GRID          39       0      1.      6.      1.       0        
GRID          40       0      1.      6.      0.       0        
GRID          41       0      0.      7.      0.       0        
GRID          42       0      0.      7.      1.       0        
GRID          43       0      1.      7.      1.       0        
GRID          44       0      1.      7.      0.       0        
GRID          49       0      0.      8.      0.       0        
GRID          50       0      0.      8.      1.       0        
GRID          51       0      1.      8.      1.       0        
GRID          52       0      1.      8.      0.       0        
GRID          53       0      0.      9.      0.       0        
GRID          54       0      0.      9.      1.       0        
GRID          55       0      1.      9.      1.       0        
GRID          56       0      1.      9.      0.       0        
GRID          61       0      0.     10.      0.       0        
GRID          62       0      0.     10.      1.       0        
GRID          63       0      1.     10.      1.       0        
GRID          64       0      1.     10.      0.       0        
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
CTETRA        26       1      29      30      31      38                        
CTETRA        27       1      29      40      37      38                        
CTETRA        28       1      31      40      38      39                        
CTETRA        29       1      38      29      40      31                        
CTETRA        30       1      40      31      29      32                        
CTETRA        31       1      37      38      40      41                        
CTETRA        32       1      38      39      40      43                        
CTETRA        33       1      41      38      43      42                        
CTETRA        34       1      40      44      41      43                        
CTETRA        35       1      41      40      43      38                        
CTETRA        36       1      41      42      43      50                        
CTETRA        37       1      41      52      49      50                        
CTETRA        38       1      43      52      50      51                        
CTETRA        39       1      50      41      52      43                        
CTETRA        40       1      52      43      41      44                        
CTETRA        41       1      49      50      52      53                        
CTETRA        42       1      50      51      52      55                        
CTETRA        43       1      53      50      55      54                        
CTETRA        44       1      52      56      53      55                        
CTETRA        45       1      53      52      55      50                        
CTETRA        46       1      53      54      55      62                        
CTETRA        47       1      53      64      61      62                        
CTETRA        48       1      55      64      62      63                        
CTETRA        49       1      62      53      64      55                        
CTETRA        50       1      64      55      53      56                        
ENDDATA
