ID c:\Users\Dima\Lame3d2\Mpl,FEMAP
SOL SESTATICS
TIME 10000
CEND
  ECHO = NONE
  DISPLACEMENT(PLOT) = ALL
  SPCFORCE(PLOT) = ALL
  OLOAD(PLOT) = ALL
  FORCE(PLOT,CORNER) = ALL
  STRESS(PLOT,CORNER) = ALL
  LOAD = 1
BEGIN BULK
$ ***************************************************************************
$   Written by : FEMAP
$   Version    : 7.00
$   Translator : MSC/NASTRAN
$   From Model : c:\Users\Dima\Lame3d2\Mpls_all\Exe\Models_nas\Box3D\Solid\Box_model.MOD
$   Date       : Wed Dec 10 10:57:42 2008
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
$ FEMAP Load Set 1 : Source Vn -.001 0.0
PLOAD4         1       1      1.                               1       3
PLOAD4         1       2      1.                               4       7
$ FEMAP Load Set 2 : Boundary Skews
PLOAD4         2       3      2.                               9      13
PLOAD4         2       4      2.                              12      15
PLOAD4         2       3      2.                              14      33
PLOAD4         2       4      2.                              10      11
PLOAD4         2       3      2.                               9      21
PLOAD4         2       4      2.                               5      15
$ FEMAP Property 1 : Acoustic space
PSOLID         1       1       0        
$ FEMAP Property 2 : Rigid body
PSOLID         2       1       0        
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           1       0      0.      1.   -0.75       0        
GRID           2       0      0.      0.   -0.75       0        
GRID           3       0      0.      0.      0.       0        
GRID           4       0      0.      1.      0.       0        
GRID           5       0     0.8      0.      0.       0        
GRID           6       0     0.8      0.    0.75       0        
GRID           7       0      0.      0.    0.75       0        
GRID           8       0      0.      1.    0.75       0        
GRID           9       0     0.8      0.   -0.75       0        
GRID          10       0     0.8      1.      0.       0        
GRID          11       0     1.1      1.    0.75       0        
GRID          12       0     0.8      1.    0.75       0        
GRID          13       0     1.1      1.   -0.75       0        
GRID          14       0     0.8      1.   -0.75       0        
GRID          15       0     1.1      0.    0.75       0        
GRID          16       0     1.1      0.   -0.75       0        
GRID          17       0      3.      0.   -0.75       0        
GRID          18       0      3.      1.   -0.75       0        
GRID          19       0      3.      1.      0.       0        
GRID          20       0      3.      0.    0.75       0        
GRID          21       0     1.1      0.      0.       0        
GRID          22       0 2.36667      0.    0.75       0        
GRID          23       0      3.      1.    0.75       0        
GRID          24       0 2.36667      1.    0.75       0        
GRID          25       0 1.73333      1.    0.75       0        
GRID          26       0 2.36667      0.   -0.75       0        
GRID          27       0 1.73333      0.   -0.75       0        
GRID          28       0 1.73333      1.   -0.75       0        
GRID          29       0      3.      0.      0.       0        
GRID          30       0 1.73333      0.    0.75       0        
GRID          31       0 1.73333      0.      0.       0        
GRID          32       0 2.36667      0.      0.       0        
GRID          33       0     1.1      1.      0.       0        
GRID          34       0 2.36667      1.   -0.75       0        
GRID          35       0 1.73333      1.      0.       0        
GRID          36       0 2.36667      1.      0.       0        
CHEXA          1       1       1       4       3       2      14      10+EL    1
+EL    1       5       9                                                        
CHEXA          2       1       4       8       7       3      10      12+EL    2
+EL    2       6       5                                                        
CHEXA          3       2       9      14      10       5      16      13+EL    3
+EL    3      33      21                                                        
CHEXA          4       2       5      10      12       6      21      33+EL    4
+EL    4      11      15                                                        
CHEXA          5       1      16      13      33      21      27      28+EL    5
+EL    5      35      31                                                        
CHEXA          6       1      21      33      11      15      31      35+EL    6
+EL    6      25      30                                                        
CHEXA          7       1      27      28      35      31      26      34+EL    7
+EL    7      36      32                                                        
CHEXA          8       1      31      35      25      30      32      36+EL    8
+EL    8      24      22                                                        
CHEXA          9       1      26      34      36      32      17      18+EL    9
+EL    9      19      29                                                        
CHEXA         10       1      32      36      24      22      29      19+EL    A
+EL    A      23      20                                                        
ENDDATA
