ID FEMAP,FEMAP
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
$   From Model : 
$   Date       : Mon Oct 09 07:35:41 2006
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
PBEAM          1       1      0.      0.      0.      0.      0.      0.+PR    1
+PR    1      0.      0.      0.      0.      0.      0.      0.      0.+PA    1
+PA    1    YESA      1.                                                +PC    1
+PC    1                                                                        
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID         134       0    -0.5     15.    -0.5       0        
GRID         136       0      0.     16.      0.       0        
GRID         139       0     0.5     15.     0.5       0        
GRID         144       0     0.5     15.    -0.5       0        
GRID         147       0    -0.5     15.     0.5       0        
CBEAM         81       1     136     134      1.      0.      0.
CBEAM         83       1     136     139      1.      0.      0.
CBEAM         85       1     136     144      1.      0.      0.
CBEAM         87       1     136     147      1.      0.      0.
ENDDATA
