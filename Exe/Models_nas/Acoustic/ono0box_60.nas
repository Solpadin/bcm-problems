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
$   Date       : Tue Feb 05 21:35:25 2008
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
GRID           3       0     1.2 0.33333      0.       0        
GRID           4       0   0.175     0.5      0.       0        
GRID           5       0      0.     0.5      0.       0        
GRID           6       0      0. 0.33333      0.       0        
GRID           7       0      0. 0.16667      0.       0        
GRID           8       0      0.      0.      0.       0        
GRID           9       0    0.15      0.      0.       0        
GRID          10       0     0.3      0.      0.       0        
GRID          11       0    0.45      0.      0.       0        
GRID          13       0     0.6   0.125      0.       0        
GRID          14       0     1.2 0.16667      0.       0        
GRID          15       0   0.475    0.25      0.       0        
GRID          16       0 0.31994 0.13321      0.       0        
GRID          17       0 0.16964 0.15573      0.       0        
GRID          18       0 0.20878 0.32281      0.       0        
GRID          19       0 0.46121  0.1271      0.       0        
GRID          20       0     0.6    0.25      0.       0        
GRID          22       0     0.6     0.5      0.       0        
GRID          23       0   0.475     0.5      0.       0        
GRID          24       0    0.35     0.5      0.       0        
GRID          25       0    0.35   0.375      0.       0        
GRID          26       0    0.35    0.25      0.       0        
GRID          27       0     1.2      0.      0.       0        
GRID          28       0   0.475   0.375      0.       0        
GRID          29       0    1.05      0.      0.       0        
GRID          30       0     0.9      0.      0.       0        
GRID          31       0    0.75      0.      0.       0        
GRID          32       0     0.6      0.      0.       0        
GRID          34       0   0.725    0.25      0.       0        
GRID          35       0 0.88006 0.13321      0.       0        
GRID          36       0 1.03036 0.15573      0.       0        
GRID          37       0 0.99122 0.32281      0.       0        
GRID          38       0 0.73879  0.1271      0.       0        
GRID          40       0     0.6   0.375      0.       0        
GRID          44       0    0.85   0.375      0.       0        
GRID          45       0    0.85    0.25      0.       0        
GRID          46       0   0.725   0.375      0.       0        
GRID          47       0   1.025     0.5      0.       0        
GRID          48       0     1.2     0.5      0.       0        
GRID          49       0     1.2 0.66667      0.       0        
GRID          52       0      0. 0.66667      0.       0        
GRID          53       0      0. 0.83333      0.       0        
GRID          54       0      0.      1.      0.       0        
GRID          55       0    0.15      1.      0.       0        
GRID          56       0     0.3      1.      0.       0        
GRID          57       0    0.45      1.      0.       0        
GRID          59       0     0.6   0.875      0.       0        
GRID          60       0     1.2 0.83333      0.       0        
GRID          61       0   0.475    0.75      0.       0        
GRID          62       0 0.31994 0.86679      0.       0        
GRID          63       0 0.16964 0.84427      0.       0        
GRID          64       0 0.20878 0.67719      0.       0        
GRID          65       0 0.46121  0.8729      0.       0        
GRID          66       0     0.6    0.75      0.       0        
GRID          71       0    0.35   0.625      0.       0        
GRID          72       0    0.35    0.75      0.       0        
GRID          73       0     1.2      1.      0.       0        
GRID          74       0   0.475   0.625      0.       0        
GRID          75       0    1.05      1.      0.       0        
GRID          76       0     0.9      1.      0.       0        
GRID          77       0    0.75      1.      0.       0        
GRID          78       0     0.6      1.      0.       0        
GRID          80       0   0.725    0.75      0.       0        
GRID          81       0 0.88006 0.86679      0.       0        
GRID          82       0 1.03036 0.84427      0.       0        
GRID          83       0 0.99122 0.67719      0.       0        
GRID          84       0 0.73879  0.8729      0.       0        
GRID          86       0     0.6   0.625      0.       0        
GRID          88       0   0.725     0.5      0.       0        
GRID          89       0    0.85     0.5      0.       0        
GRID          90       0    0.85   0.625      0.       0        
GRID          91       0    0.85    0.75      0.       0        
GRID          92       0   0.725   0.625      0.       0        
CQUAD4         1       1       7       8       9      17                
CQUAD4         2       1       6       7      17      18                
CQUAD4         3       1       4       5       6      18                
CQUAD4         4       1      17       9      10      16                
CQUAD4         5       1      25      24       4      18                
CTRIA3         6       1      26      25      18                        
CQUAD4         7       1      26      18      17      16                
CQUAD4         8       1      13      20      15      19                
CQUAD4         9       1      11      32      13      19                
CQUAD4        10       1      16      10      11      19                
CQUAD4        11       1      26      16      19      15                
CQUAD4        12       2      20      40      28      15                
CQUAD4        13       2      40      22      23      28                
CQUAD4        14       2      15      28      25      26                
CQUAD4        15       2      28      23      24      25                
CQUAD4        16       1      27      14      36      29                
CQUAD4        17       1      14       3      37      36                
CQUAD4        18       1      48      47      37       3                
CQUAD4        19       1      29      36      35      30                
CQUAD4        20       1      89      44      37      47                
CTRIA3        21       1      44      45      37                        
CQUAD4        22       1      37      45      35      36                
CQUAD4        23       1      20      13      38      34                
CQUAD4        24       1      32      31      38      13                
CQUAD4        25       1      30      35      38      31                
CQUAD4        26       1      35      45      34      38                
CQUAD4        27       2      40      20      34      46                
CQUAD4        28       2      22      40      46      88                
CQUAD4        29       2      46      34      45      44                
CQUAD4        30       2      88      46      44      89                
CQUAD4        31       1      54      53      63      55                
CQUAD4        32       1      53      52      64      63                
CQUAD4        33       1       5       4      64      52                
CQUAD4        34       1      55      63      62      56                
CQUAD4        35       1      24      71      64       4                
CTRIA3        36       1      71      72      64                        
CQUAD4        37       1      64      72      62      63                
CQUAD4        38       1      66      59      65      61                
CQUAD4        39       1      78      57      65      59                
CQUAD4        40       1      56      62      65      57                
CQUAD4        41       1      62      72      61      65                
CQUAD4        42       2      86      66      61      74                
CQUAD4        43       2      22      86      74      23                
CQUAD4        44       2      74      61      72      71                
CQUAD4        45       2      23      74      71      24                
CQUAD4        46       1      60      73      75      82                
CQUAD4        47       1      49      60      82      83                
CQUAD4        48       1      47      48      49      83                
CQUAD4        49       1      82      75      76      81                
CQUAD4        50       1      90      89      47      83                
CTRIA3        51       1      91      90      83                        
CQUAD4        52       1      91      83      82      81                
CQUAD4        53       1      59      66      80      84                
CQUAD4        54       1      77      78      59      84                
CQUAD4        55       1      81      76      77      84                
CQUAD4        56       1      91      81      84      80                
CQUAD4        57       2      66      86      92      80                
CQUAD4        58       2      86      22      88      92                
CQUAD4        59       2      80      92      90      91                
CQUAD4        60       2      92      88      89      90                
ENDDATA
