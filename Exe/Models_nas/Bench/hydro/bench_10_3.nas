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
$   From Model : D:\Job\Dima\Geo_2009\Bench\Solid\Bench_model.MOD
$   Date       : Wed Oct 13 19:09:51 2010
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
PLOAD4         1       2      1.                               9      11
$ FEMAP Load Set 2 : End Pressure 0
PLOAD4         2       1      2.                               1       4
$ FEMAP Property 1 : Untitled
PSOLID         1       1       0        
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           1       0     -1.    -0.5     1.2       0        
GRID           2       0     -1.     0.5     1.2       0        
GRID           3       0      0.    -0.5     2.4       0        
GRID           4       0     -1.     0.5     2.4       0        
GRID           5       0     -1.    -0.5     2.4       0        
GRID           6       0      0.    -0.5     1.2       0        
GRID           7       0      0.     0.5     1.2       0        
GRID           8       0      1.     0.5     1.2       0        
GRID           9       0      0.    -0.5      0.       0        
GRID          10       0      0.     0.5      0.       0        
GRID          11       0      1.     0.5      0.       0        
GRID          12       0      1.    -0.5      0.       0        
GRID          13       0      0.     0.5     2.4       0        
GRID          14       0      1.     0.5     2.4       0        
GRID          15       0      1.    -0.5     2.4       0        
GRID          16       0      1.    -0.5     1.2       0        
CHEXA          1       1       1       2       4       5       6       7+EL    1
+EL    1      13       3                                                        
CHEXA          2       1       9      10       7       6      12      11+EL    2
+EL    2       8      16                                                        
CHEXA          3       1       6       7      13       3      16       8+EL    3
+EL    3      14      15                                                        
ENDDATA
