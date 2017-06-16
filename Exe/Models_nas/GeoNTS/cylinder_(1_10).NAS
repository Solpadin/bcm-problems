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
$   Date       : Wed Jul 16 15:11:56 2008
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
$ FEMAP Load Set 111 : Cylinder BSOURCE 0.06
$ FEMAP Property 1 : Untitled
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Untitled
PSHELL         2       1      0.       1               1              0.
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID         235       0      0.      0.      0.       0        
GRID         236       0     0.2     0.2      0.       0        
GRID         239       0      0.     0.2      0.       0        
GRID         252       0     0.4      0.      0.       0        
GRID         253       0     0.2      0.      0.       0        
GRID         254       0     0.4     0.2      0.       0        
GRID         266       0     0.6     0.2      0.       0        
GRID         277       0     0.6      0.      0.       0        
GRID         278       0     0.8     0.2      0.       0        
GRID         288       0      1.      0.      0.       0        
GRID         289       0     0.8      0.      0.       0        
GRID         290       0      1.     0.2      0.       0        
GRID         302       0     1.2     0.2      0.       0        
GRID         313       0     1.2      0.      0.       0        
GRID         325       0     1.4      0.      0.       0        
GRID         327       0     1.4     0.2      0.       0        
GRID         337       0     1.6      0.      0.       0        
GRID         339       0     1.6     0.2      0.       0        
GRID         348       0      2.      0.      0.       0        
GRID         349       0     1.8      0.      0.       0        
GRID         350       0      2.     0.2      0.       0        
GRID         351       0     1.8     0.2      0.       0        
CQUAD4       201       1     253     236     239     235                
CQUAD4       206       1     252     254     236     253                
CQUAD4       211       1     277     266     254     252                
CQUAD4       216       1     289     278     266     277                
CQUAD4       221       1     288     290     278     289                
CQUAD4       226       1     313     302     290     288                
CQUAD4       231       1     325     327     302     313                
CQUAD4       236       1     337     339     327     325                
CQUAD4       241       1     349     351     339     337                
CQUAD4       246       1     348     350     351     349                
ENDDATA
