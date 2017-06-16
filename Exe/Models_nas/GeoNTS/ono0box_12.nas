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
$   Date       : Tue Feb 05 21:46:33 2008
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
GRID           2       0     1.2      0.      0.       0        
GRID           3       0     1.2    0.25      0.       0        
GRID           4       0    0.85    0.25      0.       0        
GRID           5       0    0.85      0.      0.       0        
GRID           9       0    0.85     0.5      0.       0        
GRID          10       0     1.2     0.5      0.       0        
GRID          11       0     1.2      1.      0.       0        
GRID          12       0     1.2    0.75      0.       0        
GRID          13       0    0.85    0.75      0.       0        
GRID          14       0    0.85      1.      0.       0        
GRID          15       0     0.6      1.      0.       0        
GRID          19       0      0.     0.5      0.       0        
GRID          20       0      0.      1.      0.       0        
GRID          21       0      0.    0.75      0.       0        
GRID          22       0    0.35    0.75      0.       0        
GRID          23       0    0.35      1.      0.       0        
GRID          25       0     0.6    0.75      0.       0        
GRID          27       0    0.35     0.5      0.       0        
GRID          33       0      0.      0.      0.       0        
GRID          36       0      0.    0.25      0.       0        
GRID          37       0    0.35    0.25      0.       0        
GRID          38       0    0.35      0.      0.       0        
GRID          39       0     0.6      0.      0.       0        
GRID          41       0     0.6    0.25      0.       0        
CQUAD4         1       1       4       3      10       9                
CQUAD4         2       1       5       2       3       4                
CQUAD4         3       1       5       4      41      39                
CQUAD4         5       1      12      13       9      10                
CQUAD4         6       1      11      14      13      12                
CQUAD4         7       1      13      14      15      25                
CQUAD4         9       1      22      21      19      27                
CQUAD4        10       1      23      20      21      22                
CQUAD4        11       1      23      22      25      15                
CQUAD4        16       1      36      37      27      19                
CQUAD4        17       1      33      38      37      36                
CQUAD4        18       1      37      38      39      41                
ENDDATA
