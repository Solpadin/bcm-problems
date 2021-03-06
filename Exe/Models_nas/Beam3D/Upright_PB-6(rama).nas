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
  LOAD = 222
BEGIN BULK
$ ***************************************************************************
$   Written by : FEMAP
$   Version    : 7.00
$   Translator : MSC/NASTRAN
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\Beam3D\Solid\Upright_PB-6.MOD
$   Date       : Sun Mar 25 07:53:26 2007
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
$ FEMAP Load Set 111 : Basis Rigid
SPCD         111     335       1    111.
SPCD         111     333       1    111.
$ FEMAP Load Set 222 : Cable Marker1
SPCD         222     357       1    222.
$ FEMAP Property 3 : L50
PBEAM          3       1      0.      0.      0.      0.      0.      0.+PR    3
+PR    3      0.      0.      0.      0.      0.      0.      0.      0.+PA    3
+PA    3    YESA      1.                                                +PC    3
+PC    3                                                                        
$ FEMAP Property 5 : C10
PBEAM          5       1      0.      0.      0.      0.      0.      0.+PR    5
+PR    5      0.      0.      0.      0.      0.      0.      0.      0.+PA    5
+PA    5    YESA      1.                                                +PC    5
+PC    5                                                                        
$ FEMAP Property 6 : B18
PBEAM          6       1      0.      0.      0.      0.      0.      0.+PR    6
+PR    6      0.      0.      0.      0.      0.      0.      0.      0.+PA    6
+PA    6    YESA      1.                                                +PC    6
+PC    6                                                                        
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID         333       0  0.3125      0.      0.       0        
GRID         335       0 -0.3125      0.      0.       0        
GRID         337       0 0.26384   0.375 0.34821       0        
GRID         339       0-0.21518    0.75 0.69643       0        
GRID         341       0 0.16652   1.125 1.04464       0        
GRID         343       0 0.11786     1.5 1.39286       0        
GRID         345       0-0.11786     1.5 1.39286       0        
GRID         347       0-8.87E-2   1.725 1.60179       0        
GRID         349       00.088661   1.725 1.60179       0        
GRID         351       00.059464    1.95 1.81071       0        
GRID         353       0-5.95E-2    1.95 1.81071       0        
GRID         355       0   -0.04     2.1    1.95       0        
GRID         357       0      0.     2.1    1.95       0        
GRID         358       0    0.04     2.1    1.95       0        
CBEAM        196       5     335     333 0.57735 0.57735 0.57735
CBEAM        197       5     333     337 0.57735 0.57735 0.57735
CBEAM        198       5     335     339 0.57735 0.57735 0.57735
CBEAM        199       5     337     341 0.57735 0.57735 0.57735
CBEAM        200       5     339     345 0.57735 0.57735 0.57735
CBEAM        201       5     341     343 0.57735 0.57735 0.57735
CBEAM        202       5     343     349 0.57735 0.57735 0.57735
CBEAM        203       5     345     347 0.57735 0.57735 0.57735
CBEAM        204       5     347     353 0.57735 0.57735 0.57735
CBEAM        205       5     349     351 0.57735 0.57735 0.57735
CBEAM        206       5     351     358 0.57735 0.57735 0.57735
CBEAM        207       5     353     355 0.57735 0.57735 0.57735
CBEAM        208       5     355     357 0.57735 0.57735 0.57735
CBEAM        209       5     357     358 0.57735 0.57735 0.57735
ENDDATA
