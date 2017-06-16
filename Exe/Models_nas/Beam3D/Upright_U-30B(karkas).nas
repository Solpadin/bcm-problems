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
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\Beam3D\Solid\Upright_U-30B.MOD
$   Date       : Sat Mar 24 12:43:16 2007
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
SPCD         111     377       1    111.
SPCD         111     375       1    111.
SPCD         111     419       1    111.
SPCD         111     421       1    111.
$ FEMAP Load Set 222 : Cable Marker1
SPCD         222     418       1    222.
$ FEMAP Property 1 : C16
PBEAM          1       1      0.      0.      0.      0.      0.      0.+PR    1
+PR    1      0.      0.      0.      0.      0.      0.      0.      0.+PA    1
+PA    1    YESA      1.                                                +PC    1
+PC    1                                                                        
$ FEMAP Property 2 : 2õL75
PBEAM          2       1      0.      0.      0.      0.      0.      0.+PR    2
+PR    2      0.      0.      0.      0.      0.      0.      0.      0.+PA    2
+PA    2    YESA      1.                                                +PC    2
+PC    2                                                                        
$ FEMAP Property 3 : L50
PBEAM          3       1      0.      0.      0.      0.      0.      0.+PR    3
+PR    3      0.      0.      0.      0.      0.      0.      0.      0.+PA    3
+PA    3    YESA      1.                                                +PC    3
+PC    3                                                                        
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID         375       0    0.75      0.      0.       0        
GRID         377       0   -0.75      0.      0.       0        
GRID         379       0 0.743260.049544      0.       0        
GRID         380       0 0.59125 1.16726      0.       0        
GRID         381       0-0.59664 1.12763      0.       0        
GRID         383       0-0.59125 1.16726      0.       0        
GRID         387       0 0.43803  2.2939      0.       0        
GRID         388       0 0.43264 2.33353      0.       0        
GRID         391       0-0.58586  1.2069      0.       0        
GRID         393       0-0.43264 2.33353      0.       0        
GRID         395       0 0.42725 2.37317      0.       0        
GRID         401       0-0.27807 3.47008      0.       0        
GRID         403       0-0.26998 3.52953      0.       0        
GRID         404       0-0.18239  4.1736      0.       0        
GRID         407       0 0.18239  4.1736      0.       0        
GRID         408       0 0.17835 4.20333      0.       0        
GRID         409       0 0.18643 4.14388      0.       0        
GRID         413       0-9.48E-2 4.81768      0.       0        
GRID         414       0   -0.07      5.      0.       0        
GRID         417       0    0.07      5.      0.       0        
GRID         418       0      0.      5.      0.       0        
GRID         419       0   -0.75      0.    1.22       0        
GRID         420       0-0.27403     3.5      0.       0        
GRID         421       0    0.75      0.    1.22       0        
GRID         422       0 0.27403     3.5      0.       0        
CBEAM        244       1     375     379 0.57735 0.57735 0.57735
CBEAM        245       1     377     381 0.57735 0.57735 0.57735
CBEAM        246       1     379     380 0.57735 0.57735 0.57735
CBEAM        247       1     381     383 0.57735 0.57735 0.57735
CBEAM        248       1     383     391 0.57735 0.57735 0.57735
CBEAM        249       1     380     387 0.57735 0.57735 0.57735
CBEAM        250       1     387     388 0.57735 0.57735 0.57735
CBEAM        251       1     388     395 0.57735 0.57735 0.57735
CBEAM        252       1     391     393 0.57735 0.57735 0.57735
CBEAM        253       1     393     401 0.57735 0.57735 0.57735
CBEAM        254       1     395     422 0.57735 0.57735 0.57735
CBEAM        255       1     422     409 0.57735 0.57735 0.57735
CBEAM        256       1     420     403 0.57735 0.57735 0.57735
CBEAM        257       1     401     420 0.57735 0.57735 0.57735
CBEAM        258       1     403     404 0.57735 0.57735 0.57735
CBEAM        259       1     404     413 0.57735 0.57735 0.57735
CBEAM        260       1     407     408 0.57735 0.57735 0.57735
CBEAM        261       1     409     407 0.57735 0.57735 0.57735
CBEAM        262       1     408     417 0.57735 0.57735 0.57735
CBEAM        263       1     413     414 0.57735 0.57735 0.57735
CBEAM        264       1     414     418 0.57735 0.57735 0.57735
CBEAM        265       1     417     418 0.57735 0.57735 0.57735
CBEAM        266       2     419     420 0.57735 0.57735 0.57735
CBEAM        267       2     421     422 0.57735 0.57735 0.57735
ENDDATA
