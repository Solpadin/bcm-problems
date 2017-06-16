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
$   Date       : Wed Jul 16 15:08:41 2008
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
GRID         240       0     0.2     0.4      0.       0        
GRID         243       0      0.     0.4      0.       0        
GRID         244       0     0.2     0.6      0.       0        
GRID         247       0      0.     0.6      0.       0        
GRID         250       0      0.      1.      0.       0        
GRID         251       0      0.     0.8      0.       0        
GRID         252       0     0.4      0.      0.       0        
GRID         253       0     0.2      0.      0.       0        
GRID         254       0     0.4     0.2      0.       0        
GRID         262       0     0.2      1.      0.       0        
GRID         263       0     0.2     0.8      0.       0        
GRID         266       0     0.6     0.2      0.       0        
GRID         268       0     0.6     0.4      0.       0        
GRID         269       0     0.4     0.4      0.       0        
GRID         271       0     0.4     0.6      0.       0        
GRID         274       0     0.4      1.      0.       0        
GRID         275       0     0.4     0.8      0.       0        
GRID         277       0     0.6      0.      0.       0        
GRID         278       0     0.8     0.2      0.       0        
GRID         280       0     0.8     0.4      0.       0        
GRID         283       0     0.6     0.6      0.       0        
GRID         286       0     0.6      1.      0.       0        
GRID         287       0     0.6     0.8      0.       0        
GRID         288       0      1.      0.      0.       0        
GRID         289       0     0.8      0.      0.       0        
GRID         290       0      1.     0.2      0.       0        
GRID         292       0      1.     0.4      0.       0        
GRID         295       0     0.8     0.6      0.       0        
GRID         298       0     0.8      1.      0.       0        
GRID         299       0     0.8     0.8      0.       0        
GRID         302       0     1.2     0.2      0.       0        
GRID         304       0     1.2     0.4      0.       0        
GRID         307       0      1.     0.6      0.       0        
GRID         310       0      1.      1.      0.       0        
GRID         311       0      1.     0.8      0.       0        
GRID         313       0     1.2      0.      0.       0        
GRID         316       0     1.4     0.4      0.       0        
GRID         318       0     1.4     0.6      0.       0        
GRID         319       0     1.2     0.6      0.       0        
GRID         320       0     1.4     0.8      0.       0        
GRID         321       0     1.4      1.      0.       0        
GRID         322       0     1.2      1.      0.       0        
GRID         323       0     1.2     0.8      0.       0        
GRID         325       0     1.4      0.      0.       0        
GRID         327       0     1.4     0.2      0.       0        
GRID         328       0     1.6     0.4      0.       0        
GRID         337       0     1.6      0.      0.       0        
GRID         339       0     1.6     0.2      0.       0        
GRID         340       0     1.8     0.4      0.       0        
GRID         343       0     1.6     0.6      0.       0        
GRID         346       0     1.6      1.      0.       0        
GRID         347       0     1.6     0.8      0.       0        
GRID         348       0      2.      0.      0.       0        
GRID         349       0     1.8      0.      0.       0        
GRID         350       0      2.     0.2      0.       0        
GRID         351       0     1.8     0.2      0.       0        
GRID         352       0      2.     0.4      0.       0        
GRID         354       0      2.     0.6      0.       0        
GRID         355       0     1.8     0.6      0.       0        
GRID         356       0      2.     0.8      0.       0        
GRID         357       0      2.      1.      0.       0        
GRID         358       0     1.8      1.      0.       0        
GRID         359       0     1.8     0.8      0.       0        
CQUAD4       201       1     253     236     239     235                
CQUAD4       202       1     236     240     243     239                
CQUAD4       203       1     240     244     247     243                
CQUAD4       204       1     244     263     251     247                
CQUAD4       205       1     263     262     250     251                
CQUAD4       206       1     252     254     236     253                
CQUAD4       207       1     254     269     240     236                
CQUAD4       208       1     269     271     244     240                
CQUAD4       209       1     271     275     263     244                
CQUAD4       210       1     275     274     262     263                
CQUAD4       211       1     277     266     254     252                
CQUAD4       212       1     266     268     269     254                
CQUAD4       213       1     268     283     271     269                
CQUAD4       214       1     283     287     275     271                
CQUAD4       215       1     287     286     274     275                
CQUAD4       216       1     289     278     266     277                
CQUAD4       217       1     278     280     268     266                
CQUAD4       218       1     280     295     283     268                
CQUAD4       219       1     295     299     287     283                
CQUAD4       220       1     299     298     286     287                
CQUAD4       221       1     288     290     278     289                
CQUAD4       222       1     290     292     280     278                
CQUAD4       223       1     292     307     295     280                
CQUAD4       224       1     307     311     299     295                
CQUAD4       225       1     311     310     298     299                
CQUAD4       226       1     313     302     290     288                
CQUAD4       227       1     302     304     292     290                
CQUAD4       228       1     304     319     307     292                
CQUAD4       229       1     319     323     311     307                
CQUAD4       230       1     323     322     310     311                
CQUAD4       231       1     325     327     302     313                
CQUAD4       232       1     327     316     304     302                
CQUAD4       233       1     316     318     319     304                
CQUAD4       234       1     318     320     323     319                
CQUAD4       235       1     320     321     322     323                
CQUAD4       236       1     337     339     327     325                
CQUAD4       237       1     339     328     316     327                
CQUAD4       238       1     328     343     318     316                
CQUAD4       239       1     343     347     320     318                
CQUAD4       240       1     347     346     321     320                
CQUAD4       241       1     349     351     339     337                
CQUAD4       242       1     351     340     328     339                
CQUAD4       243       1     340     355     343     328                
CQUAD4       244       1     355     359     347     343                
CQUAD4       245       1     359     358     346     347                
CQUAD4       246       1     348     350     351     349                
CQUAD4       247       1     350     352     340     351                
CQUAD4       248       1     352     354     355     340                
CQUAD4       249       1     354     356     359     355                
CQUAD4       250       1     356     357     358     359                
ENDDATA
