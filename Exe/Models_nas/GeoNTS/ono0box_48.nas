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
$   Date       : Tue Feb 05 21:50:39 2008
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
PSHELL         1       1      0.       1      1.       1 0.83333      0.
$ FEMAP Property 2 : Untitled
PSHELL         2       1      0.       1      1.       1 0.83333      0.
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           5       0    0.85   0.375      0.       0        
GRID           7       0   1.025    0.25      0.       0        
GRID           8       0    0.85    0.25      0.       0        
GRID           9       0   1.025     0.5      0.       0        
GRID          11       0     1.2   0.375      0.       0        
GRID          12       0   1.025   0.375      0.       0        
GRID          13       0     1.2      0.      0.       0        
GRID          14       0   1.025      0.      0.       0        
GRID          15       0    0.85      0.      0.       0        
GRID          16       0     1.2    0.25      0.       0        
GRID          17       0     1.2   0.125      0.       0        
GRID          18       0   1.025   0.125      0.       0        
GRID          19       0    0.85   0.125      0.       0        
GRID          20       0   0.725      0.      0.       0        
GRID          22       0     0.6   0.125      0.       0        
GRID          24       0   0.725    0.25      0.       0        
GRID          25       0   0.725   0.125      0.       0        
GRID          29       0    0.85     0.5      0.       0        
GRID          30       0    0.85   0.625      0.       0        
GRID          32       0   1.025    0.75      0.       0        
GRID          33       0    0.85    0.75      0.       0        
GRID          35       0     1.2     0.5      0.       0        
GRID          36       0     1.2   0.625      0.       0        
GRID          37       0   1.025   0.625      0.       0        
GRID          38       0     1.2      1.      0.       0        
GRID          39       0   1.025      1.      0.       0        
GRID          40       0    0.85      1.      0.       0        
GRID          41       0     1.2    0.75      0.       0        
GRID          42       0     1.2   0.875      0.       0        
GRID          43       0   1.025   0.875      0.       0        
GRID          44       0    0.85   0.875      0.       0        
GRID          45       0   0.725      1.      0.       0        
GRID          46       0     0.6      1.      0.       0        
GRID          47       0     0.6   0.875      0.       0        
GRID          48       0     0.6    0.75      0.       0        
GRID          49       0   0.725    0.75      0.       0        
GRID          50       0   0.725   0.875      0.       0        
GRID          55       0    0.35   0.625      0.       0        
GRID          57       0   0.175    0.75      0.       0        
GRID          58       0    0.35    0.75      0.       0        
GRID          59       0   0.175     0.5      0.       0        
GRID          61       0      0.   0.625      0.       0        
GRID          62       0   0.175   0.625      0.       0        
GRID          63       0      0.      1.      0.       0        
GRID          64       0   0.175      1.      0.       0        
GRID          65       0    0.35      1.      0.       0        
GRID          69       0    0.35     0.5      0.       0        
GRID          70       0    0.35   0.375      0.       0        
GRID          71       0      0.    0.75      0.       0        
GRID          72       0      0.   0.875      0.       0        
GRID          74       0   0.175   0.875      0.       0        
GRID          75       0   0.175    0.25      0.       0        
GRID          76       0    0.35    0.25      0.       0        
GRID          77       0    0.35   0.875      0.       0        
GRID          78       0   0.475      1.      0.       0        
GRID          80       0      0.     0.5      0.       0        
GRID          81       0      0.   0.375      0.       0        
GRID          82       0   0.175   0.375      0.       0        
GRID          83       0      0.      0.      0.       0        
GRID          84       0   0.175      0.      0.       0        
GRID          85       0    0.35      0.      0.       0        
GRID          89       0      0.    0.25      0.       0        
GRID          90       0      0.   0.125      0.       0        
GRID          91       0   0.175   0.125      0.       0        
GRID          92       0   0.475    0.75      0.       0        
GRID          93       0    0.35   0.125      0.       0        
GRID          94       0   0.475   0.875      0.       0        
GRID          95       0   0.475      0.      0.       0        
GRID          96       0     0.6      0.      0.       0        
GRID          98       0     0.6    0.25      0.       0        
GRID          99       0   0.475    0.25      0.       0        
GRID         100       0   0.475   0.125      0.       0        
CQUAD4         5       1       7      16      11      12                
CQUAD4         6       1       8       7      12       5                
CQUAD4         7       1      12      11      35       9                
CQUAD4         8       1       5      12       9      29                
CQUAD4         9       1      14      13      17      18                
CQUAD4        10       1      15      14      18      19                
CQUAD4        11       1      18      17      16       7                
CQUAD4        12       1      19      18       7       8                
CQUAD4        13       1      19       8      24      25                
CQUAD4        14       1      15      19      25      20                
CQUAD4        15       1      25      24      98      22                
CQUAD4        16       1      20      25      22      96                
CQUAD4        21       1      41      32      37      36                
CQUAD4        22       1      32      33      30      37                
CQUAD4        23       1      36      37       9      35                
CQUAD4        24       1      37      30      29       9                
CQUAD4        25       1      38      39      43      42                
CQUAD4        30       1      89      75      82      81                
CQUAD4        31       1      75      76      70      82                
CQUAD4        32       1      81      82      59      80                
CQUAD4        33       1      82      70      69      59                
CQUAD4        34       1      83      84      91      90                
CQUAD4        35       1      84      85      93      91                
CQUAD4        36       1      90      91      75      89                
CQUAD4        37       1      91      93      76      75                
CQUAD4        38       1      76      93     100      99                
CQUAD4        39       1      93      85      95     100                
CQUAD4        40       1      99     100      22      98                
CQUAD4        41       1     100      95      96      22                
CQUAD4        42       1      39      40      44      43                
CQUAD4        43       1      42      43      32      41                
CQUAD4        44       1      43      44      33      32                
CQUAD4        45       1      33      44      50      49                
CQUAD4        46       1      44      40      45      50                
CQUAD4        47       1      49      50      47      48                
CQUAD4        48       1      50      45      46      47                
CQUAD4        53       1      57      71      61      62                
CQUAD4        54       1      58      57      62      55                
CQUAD4        55       1      62      61      80      59                
CQUAD4        56       1      55      62      59      69                
CQUAD4        57       1      64      63      72      74                
CQUAD4        58       1      65      64      74      77                
CQUAD4        59       1      74      72      71      57                
CQUAD4        60       1      77      74      57      58                
CQUAD4        61       1      77      58      92      94                
CQUAD4        62       1      65      77      94      78                
CQUAD4        63       1      94      92      48      47                
CQUAD4        64       1      78      94      47      46                
ENDDATA
