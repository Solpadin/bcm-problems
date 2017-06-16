ID D:\Job\Dima\Geo_2009\Benc,FEMAP
SOL SESTATICS
TIME 10000
CEND
  ECHO = NONE
  DISPLACEMENT(PLOT) = ALL
  SPCFORCE(PLOT) = ALL
  OLOAD(PLOT) = ALL
  FORCE(PLOT,CORNER) = ALL
  STRESS(PLOT,CORNER) = ALL
  LOAD = 2
BEGIN BULK
$ ***************************************************************************
$   Written by : FEMAP
$   Version    : 7.00
$   Translator : MSC/NASTRAN
$   From Model : D:\Job\Dima\Geo_2009\Bench\Solid\Zet_model.MOD
$   Date       : Wed Oct 13 19:04:46 2010
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
$ FEMAP Load Set 1 : Source Pressure 1
PLOAD4         1       2      1.                               7       9
$ FEMAP Load Set 2 : End Pressure 0
PLOAD4         2       4      2.                              18      16
$ FEMAP Property 1 : Untitled
PSOLID         1       1       0        
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           1       0     -1.    -0.5     1.2       0        
GRID           2       0     -1.     0.5     1.2       0        
GRID           3       0      0.    -0.5     2.4       0        
GRID           4       0     -1.     0.5     2.4       0        
GRID           5       0      0.    -0.5     1.2       0        
GRID           6       0      1.     0.5     1.2       0        
GRID           7       0      0.    -0.5      0.       0        
GRID           8       0      0.     0.5      0.       0        
GRID           9       0      1.     0.5      0.       0        
GRID          10       0      1.    -0.5      0.       0        
GRID          11       0      1.     0.5     2.4       0        
GRID          12       0      0.     0.5     1.2       0        
GRID          13       0      1.    -0.5     2.4       0        
GRID          14       0      1.    -0.5     1.2       0        
GRID          15       0     -1.    -0.5     2.4       0        
GRID          16       0      0.    -0.5     3.6       0        
GRID          17       0      0.     0.5     3.6       0        
GRID          18       0     -1.     0.5     3.6       0        
GRID          19       0     -1.    -0.5     3.6       0        
GRID          20       0      0.     0.5     2.4       0        
CHEXA          1       1       1       2       4      15       5      12+EL    1
+EL    1      20       3                                                        
CHEXA          2       1       7       8      12       5      10       9+EL    2
+EL    2       6      14                                                        
CHEXA          3       1       5      12      20       3      14       6+EL    3
+EL    3      11      13                                                        
CHEXA          4       1      15       4      18      19       3      20+EL    4
+EL    4      17      16                                                        
ENDDATA
