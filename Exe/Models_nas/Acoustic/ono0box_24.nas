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
$   Date       : Tue Feb 05 21:33:11 2008
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
GRID           3       0     1.2      0.      0.       0        
GRID           4       0   1.025      0.      0.       0        
GRID           5       0    0.85      0.      0.       0        
GRID           6       0    0.85    0.25      0.       0        
GRID           7       0   1.025    0.25      0.       0        
GRID           8       0     1.2    0.25      0.       0        
GRID          11       0     0.6     0.5      0.       0        
GRID          12       0    0.85     0.5      0.       0        
GRID          13       0   1.025     0.5      0.       0        
GRID          14       0     1.2     0.5      0.       0        
GRID          15       0     1.2      1.      0.       0        
GRID          16       0   1.025      1.      0.       0        
GRID          17       0    0.85      1.      0.       0        
GRID          18       0    0.85    0.75      0.       0        
GRID          19       0   1.025    0.75      0.       0        
GRID          20       0     1.2    0.75      0.       0        
GRID          22       0     0.6    0.75      0.       0        
GRID          25       0   0.175     0.5      0.       0        
GRID          26       0      0.     0.5      0.       0        
GRID          27       0      0.      1.      0.       0        
GRID          28       0   0.175      1.      0.       0        
GRID          29       0    0.35      1.      0.       0        
GRID          30       0    0.35    0.75      0.       0        
GRID          31       0   0.175    0.75      0.       0        
GRID          32       0      0.    0.75      0.       0        
GRID          33       0     0.6      1.      0.       0        
GRID          36       0    0.35     0.5      0.       0        
GRID          51       0      0.      0.      0.       0        
GRID          52       0   0.175      0.      0.       0        
GRID          53       0    0.35      0.      0.       0        
GRID          54       0    0.35    0.25      0.       0        
GRID          55       0   0.175    0.25      0.       0        
GRID          56       0      0.    0.25      0.       0        
GRID          59       0     0.6      0.      0.       0        
GRID          60       0     0.6    0.25      0.       0        
CQUAD4         1       1       7       8      14      13                
CQUAD4         2       1       6       7      13      12                
CQUAD4         3       1       4       3       8       7                
CQUAD4         4       1       5       4       7       6                
CQUAD4         5       1       5       6      60      59                
CQUAD4         6       2      11      60       6      12                
CQUAD4         7       1      20      19      13      14                
CQUAD4         8       1      19      18      12      13                
CQUAD4         9       1      15      16      19      20                
CQUAD4        10       1      16      17      18      19                
CQUAD4        11       1      18      17      33      22                
CQUAD4        12       2      22      11      12      18                
CQUAD4        13       1      31      32      26      25                
CQUAD4        14       1      30      31      25      36                
CQUAD4        15       1      28      27      32      31                
CQUAD4        16       1      29      28      31      30                
CQUAD4        17       1      29      30      22      33                
CQUAD4        18       2      11      22      30      36                
CQUAD4        20       1      56      55      25      26                
CQUAD4        21       1      55      54      36      25                
CQUAD4        22       1      51      52      55      56                
CQUAD4        23       1      52      53      54      55                
CQUAD4        24       1      54      53      59      60                
CQUAD4        25       2      60      11      36      54                
ENDDATA
