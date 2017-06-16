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
$   Date       : Tue Jul 15 18:24:51 2008
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
$ FEMAP Load Set 111 : Cylinder BSOURCE 0.03
$ FEMAP Property 1 : Untitled
PSHELL         1       1      0.       1      1.       1 0.83333      0.
$ FEMAP Property 2 : Untitled
PSHELL         2       1      0.       1      1.       1 0.83333      0.
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           1       0      0.      0.      0.       0        
GRID           2       0     0.1      0.      0.       0        
GRID           3       0     0.2      0.      0.       0        
GRID           4       0     0.1     0.1      0.       0        
GRID           5       0     0.3      0.      0.       0        
GRID           6       0     0.4     0.1      0.       0        
GRID           7       0     0.3     0.1      0.       0        
GRID           8       0     0.4      0.      0.       0        
GRID           9       0     0.5      0.      0.       0        
GRID          10       0     0.6      0.      0.       0        
GRID          11       0     0.7     0.1      0.       0        
GRID          12       0     0.6     0.1      0.       0        
GRID          13       0     0.7      0.      0.       0        
GRID          14       0     0.8      0.      0.       0        
GRID          15       0     0.9      0.      0.       0        
GRID          16       0      1.      0.      0.       0        
GRID          17       0     1.1      0.      0.       0        
GRID          18       0     1.2      0.      0.       0        
GRID          19       0     1.3      0.      0.       0        
GRID          20       0     1.4      0.      0.       0        
GRID          21       0     1.3     0.1      0.       0        
GRID          22       0     1.5      0.      0.       0        
GRID          23       0     1.6      0.      0.       0        
GRID          24       0     1.7     0.1      0.       0        
GRID          25       0     1.7      0.      0.       0        
GRID          26       0     1.8      0.      0.       0        
GRID          27       0     1.9      0.      0.       0        
GRID          28       0     1.9     0.1      0.       0        
GRID          29       0     1.8     0.1      0.       0        
GRID          30       0      2.      0.      0.       0        
GRID          31       0      2.     0.1      0.       0        
GRID          32       0      0.     0.1      0.       0        
GRID          35       0     0.2     0.1      0.       0        
GRID          36       0     0.5     0.1      0.       0        
GRID          39       0     0.8     0.1      0.       0        
GRID          41       0     0.9     0.1      0.       0        
GRID          42       0      1.     0.1      0.       0        
GRID          43       0     1.1     0.1      0.       0        
GRID          45       0     1.2     0.1      0.       0        
GRID          47       0     1.4     0.1      0.       0        
GRID          48       0     1.5     0.1      0.       0        
GRID          49       0     1.6     0.1      0.       0        
CQUAD4         1       1       1       2       4      32                
CQUAD4         2       1       2       3      35       4                
CQUAD4         3       1       3       5       7      35                
CQUAD4         4       1       5       8       6       7                
CQUAD4         5       1       8       9      36       6                
CQUAD4         6       1       9      10      12      36                
CQUAD4         7       1      10      13      11      12                
CQUAD4         8       1      13      14      39      11                
CQUAD4         9       1      14      15      41      39                
CQUAD4        10       1      15      16      42      41                
CQUAD4        11       1      16      17      43      42                
CQUAD4        12       1      17      18      45      43                
CQUAD4        13       1      18      19      21      45                
CQUAD4        14       1      19      20      47      21                
CQUAD4        15       1      20      22      48      47                
CQUAD4        16       1      22      23      49      48                
CQUAD4        17       1      23      25      24      49                
CQUAD4        18       1      25      26      29      24                
CQUAD4        19       1      26      27      28      29                
CQUAD4        20       1      27      30      31      28                
ENDDATA
