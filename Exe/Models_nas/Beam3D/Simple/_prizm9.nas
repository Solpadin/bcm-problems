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
$   Date       : Mon Oct 09 07:43:55 2006
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
GRID         128       0      1.     14.     -1.       0        
GRID         134       0    -0.5     15.    -0.5       0        
GRID         139       0     0.5     15.     0.5       0        
GRID         144       0     0.5     15.    -0.5       0        
GRID         147       0    -0.5     15.     0.5       0        
GRID         331       0     -1.     14.      1.       0        
GRID         333       0      1.     14.      1.       0        
GRID         364       0     -1.     14.     -1.       0        
GRID         366       0      0.     16.      0.       0        
GRID         367       0      0.     17.      0.       0        
CBEAM         81       1     366     134      1.      0.      0.
CBEAM         82       1     134     364      1.      0.      0.
CBEAM         83       1     366     139      1.      0.      0.
CBEAM         84       1     139     333      1.      0.      0.
CBEAM         85       1     366     144      1.      0.      0.
CBEAM         86       1     144     128      1.      0.      0.
CBEAM         87       1     366     147      1.      0.      0.
CBEAM         88       1     147     331      1.      0.      0.
CBEAM        221       1     366     367      1.      0.      0.
ENDDATA
