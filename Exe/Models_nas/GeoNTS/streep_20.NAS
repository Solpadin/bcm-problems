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
$   Date       : Mon Oct 05 00:12:25 2009
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
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           1       0      0.      1.      0.       0        
GRID           2       0      0.      0.      0.       0        
GRID           3       0     0.2      0.      0.       0        
GRID           4       0     0.2      1.      0.       0        
GRID           5       0     0.1      1.      0.       0        
GRID           6       0     0.1      0.      0.       0        
GRID           7       0     0.3      0.      0.       0        
GRID           8       0     0.4      0.      0.       0        
GRID           9       0     0.4      1.      0.       0        
GRID          10       0     0.3      1.      0.       0        
GRID          11       0     0.5      0.      0.       0        
GRID          12       0     0.6      0.      0.       0        
GRID          13       0     0.6      1.      0.       0        
GRID          14       0     0.5      1.      0.       0        
GRID          15       0     0.8      1.      0.       0        
GRID          16       0     0.7      1.      0.       0        
GRID          17       0     0.7      0.      0.       0        
GRID          18       0     0.8      0.      0.       0        
GRID          19       0      1.      0.      0.       0        
GRID          20       0     0.9      1.      0.       0        
GRID          21       0     0.9      0.      0.       0        
GRID          22       0      1.      1.      0.       0        
GRID          23       0     1.2      1.      0.       0        
GRID          24       0     1.1      1.      0.       0        
GRID          25       0     1.1      0.      0.       0        
GRID          26       0     1.2      0.      0.       0        
GRID          27       0     1.3      1.      0.       0        
GRID          28       0     1.3      0.      0.       0        
GRID          29       0     1.4      1.      0.       0        
GRID          30       0     1.4      0.      0.       0        
GRID          31       0     1.6      0.      0.       0        
GRID          32       0     1.5      1.      0.       0        
GRID          33       0     1.5      0.      0.       0        
GRID          34       0     1.7      0.      0.       0        
GRID          35       0     1.7      1.      0.       0        
GRID          36       0     1.6      1.      0.       0        
GRID          37       0     1.8      0.      0.       0        
GRID          38       0     1.8      1.      0.       0        
GRID          39       0     1.9      0.      0.       0        
GRID          40       0      2.      0.      0.       0        
GRID          41       0      2.      1.      0.       0        
GRID          42       0     1.9      1.      0.       0        
CQUAD4         1       1       6       5       1       2                
CQUAD4         2       1       3       4       5       6                
CQUAD4         3       1       7      10       4       3                
CQUAD4         4       1       8       9      10       7                
CQUAD4         5       1      11      14       9       8                
CQUAD4         6       1      12      13      14      11                
CQUAD4         7       1      17      16      13      12                
CQUAD4         8       1      18      15      16      17                
CQUAD4         9       1      21      20      15      18                
CQUAD4        10       1      19      22      20      21                
CQUAD4        11       1      25      24      22      19                
CQUAD4        12       1      26      23      24      25                
CQUAD4        13       1      28      27      23      26                
CQUAD4        14       1      30      29      27      28                
CQUAD4        15       1      33      32      29      30                
CQUAD4        16       1      31      36      32      33                
CQUAD4        17       1      34      35      36      31                
CQUAD4        18       1      37      38      35      34                
CQUAD4        19       1      39      42      38      37                
CQUAD4        20       1      40      41      42      39                
ENDDATA
