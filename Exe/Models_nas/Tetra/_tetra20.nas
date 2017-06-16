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
$   Date       : Fri Feb 27 18:45:13 2004
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
FORCE          1     121       0      1.   1000.      0.      0.
FORCE          1     122       0      1.   1000.      0.      0.
FORCE          1     123       0      1.   1000.      0.      0.
FORCE          1     124       0      1.   1000.      0.      0.
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
GRID          65       0      0.     11.      0.       0        
GRID          66       0      0.     11.      1.       0        
GRID          67       0      1.     11.      1.       0        
GRID          68       0      1.     11.      0.       0        
GRID          73       0      0.     12.      0.       0        
GRID          74       0      0.     12.      1.       0        
GRID          75       0      1.     12.      1.       0        
GRID          76       0      1.     12.      0.       0        
GRID          77       0      0.     13.      0.       0        
GRID          78       0      0.     13.      1.       0        
GRID          79       0      1.     13.      1.       0        
GRID          80       0      1.     13.      0.       0        
GRID          85       0      0.     14.      0.       0        
GRID          86       0      0.     14.      1.       0        
GRID          87       0      1.     14.      1.       0        
GRID          88       0      1.     14.      0.       0        
GRID          89       0      0.     15.      0.       0        
GRID          90       0      0.     15.      1.       0        
GRID          91       0      1.     15.      1.       0        
GRID          92       0      1.     15.      0.       0        
GRID          96       0      1.     16.      0.       0        
GRID          97       0      0.     16.      0.       0        
GRID          98       0      0.     16.      1.       0        
GRID          99       0      1.     16.      1.       0        
GRID         101       0      0.     17.      0.       0        
GRID         102       0      0.     17.      1.       0        
GRID         103       0      1.     17.      1.       0        
GRID         104       0      1.     17.      0.       0        
GRID         109       0      0.     18.      0.       0        
GRID         110       0      0.     18.      1.       0        
GRID         111       0      1.     18.      1.       0        
GRID         112       0      1.     18.      0.       0        
GRID         113       0      0.     19.      0.       0        
GRID         114       0      0.     19.      1.       0        
GRID         115       0      1.     19.      1.       0        
GRID         116       0      1.     19.      0.       0        
GRID         121       0      0.     20.      0.       0        
GRID         122       0      0.     20.      1.       0        
GRID         123       0      1.     20.      1.       0        
GRID         124       0      1.     20.      0.       0        
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
CTETRA        51       1      61      62      64      65                        
CTETRA        52       1      62      63      64      67                        
CTETRA        53       1      65      62      67      66                        
CTETRA        54       1      64      68      65      67                        
CTETRA        55       1      65      64      67      62                        
CTETRA        56       1      65      66      67      74                        
CTETRA        57       1      65      76      73      74                        
CTETRA        58       1      67      76      74      75                        
CTETRA        59       1      74      65      76      67                        
CTETRA        60       1      76      67      65      68                        
CTETRA        61       1      73      74      76      77                        
CTETRA        62       1      74      75      76      79                        
CTETRA        63       1      77      74      79      78                        
CTETRA        64       1      76      80      77      79                        
CTETRA        65       1      77      76      79      74                        
CTETRA        66       1      77      78      79      86                        
CTETRA        67       1      77      88      85      86                        
CTETRA        68       1      79      88      86      87                        
CTETRA        69       1      86      77      88      79                        
CTETRA        70       1      88      79      77      80                        
CTETRA        71       1      85      86      88      89                        
CTETRA        72       1      86      87      88      91                        
CTETRA        73       1      89      86      91      90                        
CTETRA        74       1      88      92      89      91                        
CTETRA        75       1      89      88      91      86                        
CTETRA        76       1      89      90      91      98                        
CTETRA        77       1      89      96      97      98                        
CTETRA        78       1      91      96      98      99                        
CTETRA        79       1      98      89      96      91                        
CTETRA        80       1      96      91      89      92                        
CTETRA        81       1      97      98      96     101                        
CTETRA        82       1      98      99      96     103                        
CTETRA        83       1     101      98     103     102                        
CTETRA        84       1      96     104     101     103                        
CTETRA        85       1     101      96     103      98                        
CTETRA        86       1     101     102     103     110                        
CTETRA        87       1     101     112     109     110                        
CTETRA        88       1     103     112     110     111                        
CTETRA        89       1     110     101     112     103                        
CTETRA        90       1     112     103     101     104                        
CTETRA        91       1     109     110     112     113                        
CTETRA        92       1     110     111     112     115                        
CTETRA        93       1     113     110     115     114                        
CTETRA        94       1     112     116     113     115                        
CTETRA        95       1     113     112     115     110                        
CTETRA        96       1     113     114     115     122                        
CTETRA        97       1     113     124     121     122                        
CTETRA        98       1     115     124     122     123                        
CTETRA        99       1     122     113     124     115                        
CTETRA       100       1     124     115     113     116                        
ENDDATA
