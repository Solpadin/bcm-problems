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
$   Date       : Mon Oct 09 07:33:15 2006
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
GRID          54       0      1.     14.      0.       0        
GRID          60       0     -1.     14.      0.       0        
GRID          63       0      0.     14.      1.       0        
GRID         116       0      1.     13.      1.       0        
GRID         121       0     -1.     13.      1.       0        
GRID         126       0     -1.     13.     -1.       0        
GRID         128       0      1.     14.     -1.       0        
GRID         129       0      1.     13.     -1.       0        
GRID         134       0    -0.5     15.    -0.5       0        
GRID         136       0      0.     16.      0.       0        
GRID         139       0     0.5     15.     0.5       0        
GRID         144       0     0.5     15.    -0.5       0        
GRID         147       0    -0.5     15.     0.5       0        
GRID         331       0     -1.     14.      1.       0        
GRID         333       0      1.     14.      1.       0        
GRID         363       0      0.     14.     -1.       0        
GRID         364       0     -1.     14.     -1.       0        
CBEAM         21       1     333      54      1.      0.      0.
CBEAM         22       1      54     128      1.      0.      0.
CBEAM         23       1     364      60      1.      0.      0.
CBEAM         24       1      60     331      1.      0.      0.
CBEAM         68       1     116     333      1.      0.      0.
CBEAM         72       1     121     331      1.      0.      0.
CBEAM         76       1     126     364      1.      0.      0.
CBEAM         77       1     128     129      1.      0.      0.
CBEAM         81       1     136     134      1.      0.      0.
CBEAM         82       1     134     364      1.      0.      0.
CBEAM         83       1     136     139      1.      0.      0.
CBEAM         84       1     139     333      1.      0.      0.
CBEAM         85       1     136     144      1.      0.      0.
CBEAM         86       1     144     128      1.      0.      0.
CBEAM         87       1     136     147      1.      0.      0.
CBEAM         88       1     147     331      1.      0.      0.
CBEAM        217       1     331      63      0.      0.      1.
CBEAM        218       1      63     333      0.      0.      1.
CBEAM        219       1     128     363      0.      0.      1.
CBEAM        220       1     363     364      0.      0.      1.
ENDDATA
