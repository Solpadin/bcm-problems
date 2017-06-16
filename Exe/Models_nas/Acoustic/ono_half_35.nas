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
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\Acoustic\Solid\ono_box.MOD
$   Date       : Mon Mar 19 00:00:47 2007
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
$ FEMAP Property 1 : Ono box
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Inclusion
PSHELL         2       1      0.       1               1              0.
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID         254       0  0.2625     0.5     0.0       0        
GRID         255       0   0.175     0.5     0.0       0        
GRID         256       0  0.0875     0.5     0.0       0        
GRID         257       0      0.     0.5     0.0       0        
GRID         258       0      0.     0.4     0.0       0        
GRID         259       0      0.     0.3     0.0       0        
GRID         260       0      0.     0.2     0.0       0        
GRID         261       0      0.     0.1     0.0       0        
GRID         262       0      0.      0.     0.0       0        
GRID         263       0     0.1      0.     0.0       0        
GRID         264       0     0.2      0.     0.0       0        
GRID         265       0     0.3      0.     0.0       0        
GRID         266       0     0.4      0.     0.0       0        
GRID         267       0     0.5      0.     0.0       0        
GRID         268       0     0.6      0.     0.0       0        
GRID         269       0     0.60.083333     0.0       0        
GRID         270       0     0.6 0.16667     0.0       0        
GRID         271       0     0.6    0.25     0.0       0        
GRID         273       0    0.43    0.25     0.0       0        
GRID         275       0    0.35   0.415     0.0       0        
GRID         276       0 0.31034 0.10101     0.0       0        
GRID         277       0 0.33031 0.21012     0.0       0        
GRID         278       0 0.28509  0.3137     0.0       0        
GRID         279       0 0.26989 0.40819     0.0       0        
GRID         280       0  0.1923 0.30637     0.0       0        
GRID         281       00.095582 0.30259     0.0       0        
GRID         282       0  0.1008 0.20218     0.0       0        
GRID         283       0 0.10097 0.10101     0.0       0        
GRID         284       0 0.20618 0.20516     0.0       0        
GRID         285       0 0.20449 0.10178     0.0       0        
GRID         286       00.091261 0.40166     0.0       0        
GRID         287       0 0.18168 0.40409     0.0       0        
GRID         288       0 0.41925 0.18081     0.0       0        
GRID         289       0 0.408580.092089     0.0       0        
GRID         290       0  0.5097 0.17101     0.0       0        
GRID         291       0 0.504370.086607     0.0       0        
GRID         293       0     0.6 0.33333     0.0       0        
GRID         294       0     0.6 0.41667     0.0       0        
GRID         295       0     0.6     0.5     0.0       0        
GRID         296       0 0.51667     0.5     0.0       0        
GRID         297       0 0.43333     0.5     0.0       0        
GRID         298       0    0.35     0.5     0.0       0        
GRID         300       0    0.35    0.33     0.0       0        
GRID         302       0   0.515    0.25     0.0       0        
GRID         303       0 0.45148 0.38439     0.0       0        
GRID         304       0 0.52244 0.41161     0.0       0        
GRID         305       0 0.52211 0.34496     0.0       0        
CQUAD4       186       1     277     278     280     284                
CQUAD4       187       1     276     277     284     285                
CQUAD4       188       1     264     265     276     285                
CQUAD4       189       1     284     280     281     282                
CQUAD4       190       1     285     284     282     283                
CQUAD4       191       1     263     264     285     283                
CQUAD4       192       1     261     262     263     283                
CQUAD4       193       1     260     261     283     282                
CQUAD4       194       1     259     260     282     281                
CQUAD4       195       1     258     259     281     286                
CQUAD4       196       1     256     257     258     286                
CQUAD4       197       1     286     281     280     287                
CQUAD4       198       1     255     256     286     287                
CQUAD4       199       1     287     280     278     279                
CQUAD4       200       1     254     255     287     279                
CQUAD4       201       1     270     271     302     290                
CQUAD4       202       1     269     270     290     291                
CQUAD4       203       1     267     268     269     291                
CQUAD4       204       1     290     302     273     288                
CQUAD4       205       1     291     290     288     289                
CQUAD4       206       1     266     267     291     289                
CQUAD4       207       1     276     265     266     289                
CQUAD4       208       1     277     276     289     288                
CTRIA3       209       1     277     288     273                        
CQUAD4       210       1     278     277     273     300                
CQUAD4       211       1     279     278     300     275                
CQUAD4       212       1     298     254     279     275                
CQUAD4       213       2     297     298     275     303                
CQUAD4       214       2     296     297     303     304                
CQUAD4       215       2     294     295     296     304                
CQUAD4       216       2     303     275     300     273                
CQUAD4       217       2     303     273     302     305                
CTRIA3       218       2     304     303     305                        
CQUAD4       219       2     293     294     304     305                
CQUAD4       220       2     271     293     305     302                
ENDDATA
