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
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\Models_nas\GeoNTS\Model\model_streep2.MOD
$   Date       : Mon Oct 05 00:44:11 2009
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
$ FEMAP Property 1 : Untitled
PSHELL         1       1      0.       1               1              0.
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           1       0      0.      1.      0.       0        
GRID           2       0      0.      0.      0.       0        
GRID           3       0      2.      0.      0.       0        
GRID           4       0      2.      1.      0.       0        
GRID           5       0      1.      1.      0.       0        
GRID           6       0      1.      0.      0.       0        
CQUAD4         1       1       6       5       1       2                
CQUAD4         2       1       3       4       5       6                
ENDDATA
