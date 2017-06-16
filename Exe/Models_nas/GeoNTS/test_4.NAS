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
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\GeoNTS\Model\model_test.MOD
$   Date       : Mon Jun 09 13:10:00 2008
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
$ FEMAP Load Set 111 : Cylinder BSOURCE 0.3
$ FEMAP Property 1 : Untitled
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Untitled
PSHELL         2       1      0.       1               1              0.
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID         235       0      0.      0.      0.       0        
GRID         236       0     0.2      0.      0.       0        
GRID         237       0     0.2     0.1      0.       0        
GRID         238       0     0.1     0.1      0.       0        
GRID         239       0     0.1      0.      0.       0        
GRID         240       0      0.     0.2      0.       0        
GRID         241       0      0.     0.1      0.       0        
GRID         243       0     0.2     0.2      0.       0        
GRID         244       0     0.1     0.2      0.       0        
CQUAD4       201       1     239     238     241     235                
CQUAD4       202       1     236     237     238     239                
CQUAD4       203       1     238     244     240     241                
CQUAD4       204       1     237     243     244     238                
ENDDATA
