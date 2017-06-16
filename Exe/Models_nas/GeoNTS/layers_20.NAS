ID C:\Users\Dima\Lame3d2\Mpl,FEMAP
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
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\GeoNTS\Model\model_layer.MOD
$   Date       : Thu Mar 27 18:58:47 2008
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
$ FEMAP Property 1 : Phase 1
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Phase 2
PSHELL         2       1      0.       1               1              0.
$ FEMAP Material 1 : Abstract material
MAT1           1                              0.      0.      0.        
GRID           1       0      0.      0.      0.       0        
GRID           4       0      0.      1.      0.       0        
GRID           5       0     0.1      0.      0.       0        
GRID           8       0     0.1      1.      0.       0        
GRID           9       0     0.2      0.      0.       0        
GRID          10       0     0.2      1.      0.       0        
GRID          11       0     0.3      0.      0.       0        
GRID          14       0     0.3      1.      0.       0        
GRID          15       0     0.4      0.      0.       0        
GRID          16       0     0.4      1.      0.       0        
GRID          17       0     0.5      0.      0.       0        
GRID          18       0     0.6      0.      0.       0        
GRID          19       0     0.6      1.      0.       0        
GRID          20       0     0.5      1.      0.       0        
GRID          23       0     0.7      0.      0.       0        
GRID          25       0     0.8      1.      0.       0        
GRID          26       0     0.7      1.      0.       0        
GRID          27       0     0.8      0.      0.       0        
GRID          29       0     0.9      0.      0.       0        
GRID          31       0      1.      1.      0.       0        
GRID          32       0     0.9      1.      0.       0        
GRID          33       0      1.      0.      0.       0        
GRID          35       0     1.1      0.      0.       0        
GRID          38       0     1.1      1.      0.       0        
GRID          39       0     1.2      0.      0.       0        
GRID          40       0     1.2      1.      0.       0        
GRID          41       0     1.3      0.      0.       0        
GRID          42       0     1.4      0.      0.       0        
GRID          43       0     1.4      1.      0.       0        
GRID          44       0     1.3      1.      0.       0        
GRID          47       0     1.5      0.      0.       0        
GRID          49       0     1.6      1.      0.       0        
GRID          50       0     1.5      1.      0.       0        
GRID          51       0     1.6      0.      0.       0        
GRID          53       0     1.7      0.      0.       0        
GRID          54       0     1.8      0.      0.       0        
GRID          55       0     1.8      1.      0.       0        
GRID          56       0     1.7      1.      0.       0        
GRID          59       0     1.9      0.      0.       0        
GRID          60       0      2.      0.      0.       0        
GRID          61       0      2.      1.      0.       0        
GRID          62       0     1.9      1.      0.       0        
CQUAD4         1       1       1       5       8       4                
CQUAD4         2       2       5       9      10       8                
CQUAD4         3       1       9      11      14      10                
CQUAD4         4       2      11      15      16      14                
CQUAD4         5       1      15      17      20      16                
CQUAD4         6       2      17      18      19      20                
CQUAD4         7       1      18      23      26      19                
CQUAD4         8       2      23      27      25      26                
CQUAD4         9       1      27      29      32      25                
CQUAD4        10       2      29      33      31      32                
CQUAD4        11       1      33      35      38      31                
CQUAD4        12       2      35      39      40      38                
CQUAD4        13       1      39      41      44      40                
CQUAD4        14       2      41      42      43      44                
CQUAD4        15       1      42      47      50      43                
CQUAD4        16       2      47      51      49      50                
CQUAD4        17       1      51      53      56      49                
CQUAD4        18       2      53      54      55      56                
CQUAD4        19       1      54      59      62      55                
CQUAD4        20       2      59      60      61      62                
ENDDATA
