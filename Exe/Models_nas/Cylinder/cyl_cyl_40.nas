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
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\Cylinder\Solid\cylinder.MOD
$   Date       : Tue Dec 12 15:18:28 2006
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
$ FEMAP Property 1 : Cylinder
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Inclusion
PSHELL         2       1      0.       1               1              0.
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           1       0   -0.55      0.     0.0       0        
GRID           2       0-0.50813-0.21048     0.0       0        
GRID           3       0-0.38891-0.38891     0.0       0        
GRID           4       0-0.21048-0.50813     0.0       0        
GRID           5       0      0.   -0.55     0.0       0        
GRID           6       0 0.21048-0.50813     0.0       0        
GRID           7       0 0.38891-0.38891     0.0       0        
GRID           8       0 0.50813-0.21048     0.0       0        
GRID           9       0    0.55      0.     0.0       0        
GRID          10       0 0.50813 0.21048     0.0       0        
GRID          11       0 0.38891 0.38891     0.0       0        
GRID          12       0 0.21048 0.50813     0.0       0        
GRID          13       0      0.    0.55     0.0       0        
GRID          14       0-0.21048 0.50813     0.0       0        
GRID          15       0-0.38891 0.38891     0.0       0        
GRID          16       0-0.50813 0.21048     0.0       0        
GRID          19       0  -0.125-0.21651     0.0       0        
GRID          23       0-0.33875-0.22067     0.0       0        
GRID          24       0      0.-0.39986     0.0       0        
GRID          25       0 0.33875-0.22067     0.0       0        
GRID          26       0 0.33875 0.22067     0.0       0        
GRID          27       0      0. 0.39986     0.0       0        
GRID          28       0-0.33875 0.22067     0.0       0        
GRID          29       0   -0.25      0.     0.0       0        
GRID          31       0   0.125-0.21651     0.0       0        
GRID          32       0    0.25      0.     0.0       0        
GRID          33       0   0.125 0.21651     0.0       0        
GRID          34       0  -0.125 0.21651     0.0       0        
GRID          35       0      0.      0.     0.0       0        
CTRIA3         1       1       2       3      23                        
CTRIA3         2       1       1       2      23                        
CTRIA3         3       1      29       1      23                        
CTRIA3         4       1      23       3       4                        
CTRIA3         5       1      19      29      23                        
CTRIA3         6       1      19      23       4                        
CTRIA3         7       1      19       4      24                        
CTRIA3         8       1      31      19      24                        
CTRIA3         9       1      24       4       5                        
CTRIA3        10       1      24       5       6                        
CTRIA3        11       1      31      24       6                        
CTRIA3        12       1       7       8      25                        
CTRIA3        13       1       6       7      25                        
CTRIA3        14       1      31       6      25                        
CTRIA3        15       1      32      31      25                        
CTRIA3        16       1      25       8       9                        
CTRIA3        17       1      32      25       9                        
CTRIA3        18       1      10      11      26                        
CTRIA3        19       1       9      10      26                        
CTRIA3        20       1      32       9      26                        
CTRIA3        21       1      26      11      12                        
CTRIA3        22       1      33      32      26                        
CTRIA3        23       1      33      26      12                        
CTRIA3        24       1      33      12      27                        
CTRIA3        25       1      34      33      27                        
CTRIA3        26       1      27      12      13                        
CTRIA3        27       1      27      13      14                        
CTRIA3        28       1      34      27      14                        
CTRIA3        29       1      15      16      28                        
CTRIA3        30       1      14      15      28                        
CTRIA3        31       1      34      14      28                        
CTRIA3        32       1      29      34      28                        
CTRIA3        33       1       1      29      28                        
CTRIA3        34       1       1      28      16                        
CTRIA3        35       2      31      32      35                        
CTRIA3        36       2      19      31      35                        
CTRIA3        37       2      29      19      35                        
CTRIA3        38       2      35      32      33                        
CTRIA3        39       2      35      33      34                        
CTRIA3        40       2      29      35      34                        
ENDDATA
