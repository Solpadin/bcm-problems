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
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\GeoNTS\Model\model_cylinder.MOD
$   Date       : Wed Jul 16 17:11:16 2008
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
$ FEMAP Load Set 111 : Cylinder BSOURCE 0.12
$ FEMAP Property 1 : Untitled
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Untitled
PSHELL         2       1      0.       1               1              0.
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID        2153       0      0.     0.4      0.       0        
GRID        2154       0      0.      0.      0.       0        
GRID        2157       0     0.4     0.4      0.       0        
GRID        2158       0     0.4      0.      0.       0        
GRID        2159       0     1.2      0.      0.       0        
GRID        2161       0     0.8     0.4      0.       0        
GRID        2162       0     0.8      0.      0.       0        
GRID        2165       0     1.2     0.4      0.       0        
GRID        2167       0      2.      0.      0.       0        
GRID        2168       0      2.     0.4      0.       0        
GRID        2169       0     1.6     0.4      0.       0        
GRID        2170       0     1.6      0.      0.       0        
CQUAD4      1372       1    2158    2157    2153    2154                
CQUAD4      1373       1    2162    2161    2157    2158                
CQUAD4      1374       1    2159    2165    2161    2162                
CQUAD4      1375       1    2170    2169    2165    2159                
CQUAD4      1376       1    2167    2168    2169    2170                
ENDDATA
