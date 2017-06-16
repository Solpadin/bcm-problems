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
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\Beam3D\Solid\Upright_U-60B.MOD
$   Date       : Sun Mar 25 07:27:52 2007
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
SPCD         111     457       1    111.
SPCD         111     453       1    111.
$ FEMAP Load Set 222 : Cable Marker1
SPCD         222     496       1    222.
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
$ FEMAP Property 4 : C20
PBEAM          4       1      0.      0.      0.      0.      0.      0.+PR    4
+PR    4      0.      0.      0.      0.      0.      0.      0.      0.+PA    4
+PA    4    YESA      1.                                                +PC    4
+PC    4                                                                        
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID         453       0    0.75      0.      0.       0        
GRID         455       0 0.743260.049544      0.       0        
GRID         457       0   -0.75      0.      0.       0        
GRID         459       0-0.59664 1.12763      0.       0        
GRID         461       0-0.59125 1.16726      0.       0        
GRID         463       0-0.58586  1.2069      0.       0        
GRID         465       0 0.59125 1.16726      0.       0        
GRID         467       0 0.43803  2.2939      0.       0        
GRID         469       0 0.43264 2.33353      0.       0        
GRID         471       0 0.42725 2.37317      0.       0        
GRID         473       0-0.43264 2.33353      0.       0        
GRID         475       0-0.27403     3.5      0.       0        
GRID         477       0-0.26998 3.52953      0.       0        
GRID         479       0 0.27403     3.5      0.       0        
GRID         481       0-0.27807 3.47008      0.       0        
GRID         485       0 0.18643 4.14388      0.       0        
GRID         486       0 0.18239  4.1736      0.       0        
GRID         487       0 0.17835 4.20333      0.       0        
GRID         489       0-0.18239  4.1736      0.       0        
GRID         491       0-9.48E-2 4.81768      0.       0        
GRID         492       0   -0.07      5.      0.       0        
GRID         493       0    0.07      5.      0.       0        
GRID         496       0      0.      5.      0.       0        
CBEAM        297       4     453     455 0.57735 0.57735 0.57735
CBEAM        298       4     455     465 0.57735 0.57735 0.57735
CBEAM        299       4     457     459 0.57735 0.57735 0.57735
CBEAM        300       4     459     461 0.57735 0.57735 0.57735
CBEAM        301       4     461     463 0.57735 0.57735 0.57735
CBEAM        302       4     463     473 0.57735 0.57735 0.57735
CBEAM        303       4     465     467 0.57735 0.57735 0.57735
CBEAM        304       4     467     469 0.57735 0.57735 0.57735
CBEAM        305       4     469     471 0.57735 0.57735 0.57735
CBEAM        306       4     471     479 0.57735 0.57735 0.57735
CBEAM        307       4     473     481 0.57735 0.57735 0.57735
CBEAM        308       4     475     477 0.57735 0.57735 0.57735
CBEAM        309       4     477     489 0.57735 0.57735 0.57735
CBEAM        310       4     479     485 0.57735 0.57735 0.57735
CBEAM        311       4     481     475 0.57735 0.57735 0.57735
CBEAM        312       4     486     487 0.57735 0.57735 0.57735
CBEAM        313       4     485     486 0.57735 0.57735 0.57735
CBEAM        314       4     487     493 0.57735 0.57735 0.57735
CBEAM        315       4     489     491 0.57735 0.57735 0.57735
CBEAM        316       4     491     492 0.57735 0.57735 0.57735
CBEAM        317       4     493     496 0.57735 0.57735 0.57735
CBEAM        318       4     492     496 0.57735 0.57735 0.57735
ENDDATA
