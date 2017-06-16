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
$   Date       : Mon Jun 16 12:24:30 2008
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
GRID         238       0     0.1     0.1      0.       0        
GRID         239       0     0.1      0.      0.       0        
GRID         241       0      0.     0.1      0.       0        
GRID         246       0     0.2      0.      0.       0        
GRID         248       0     0.4     0.1      0.       0        
GRID         249       0     0.3     0.1      0.       0        
GRID         250       0     0.3      0.      0.       0        
GRID         251       0     0.2     0.2      0.       0        
GRID         252       0     0.2     0.1      0.       0        
GRID         255       0      0.     0.2      0.       0        
GRID         258       0     0.1     0.3      0.       0        
GRID         259       0     0.1     0.2      0.       0        
GRID         261       0      0.     0.3      0.       0        
GRID         263       0     0.1     0.4      0.       0        
GRID         266       0     0.4     0.3      0.       0        
GRID         267       0     0.3     0.3      0.       0        
GRID         268       0     0.3     0.2      0.       0        
GRID         270       0     0.2     0.3      0.       0        
GRID         273       0     0.4      0.      0.       0        
GRID         274       0     0.5     0.1      0.       0        
GRID         275       0     0.5      0.      0.       0        
GRID         277       0     0.6      0.      0.       0        
GRID         278       0     0.8      0.      0.       0        
GRID         279       0     0.8     0.1      0.       0        
GRID         280       0     0.7     0.1      0.       0        
GRID         281       0     0.7      0.      0.       0        
GRID         282       0     0.6     0.2      0.       0        
GRID         283       0     0.6     0.1      0.       0        
GRID         284       0     0.8     0.2      0.       0        
GRID         285       0     0.4     0.2      0.       0        
GRID         286       0     0.5     0.3      0.       0        
GRID         287       0     0.5     0.2      0.       0        
GRID         291       0     0.8     0.3      0.       0        
GRID         292       0     0.7     0.3      0.       0        
GRID         293       0     0.7     0.2      0.       0        
GRID         295       0     0.6     0.3      0.       0        
GRID         296       0     0.8     0.4      0.       0        
GRID         298       0      0.     0.4      0.       0        
GRID         299       0     0.1     0.5      0.       0        
GRID         301       0      0.     0.5      0.       0        
GRID         302       0     0.2     0.4      0.       0        
GRID         304       0     0.4     0.5      0.       0        
GRID         305       0     0.3     0.5      0.       0        
GRID         306       0     0.3     0.4      0.       0        
GRID         307       0     0.2     0.6      0.       0        
GRID         308       0     0.2     0.5      0.       0        
GRID         310       0      0.     0.6      0.       0        
GRID         311       0     0.1     0.7      0.       0        
GRID         312       0     0.1     0.6      0.       0        
GRID         313       0      0.     0.8      0.       0        
GRID         314       0      0.     0.7      0.       0        
GRID         315       0     0.1     0.8      0.       0        
GRID         317       0     0.3     0.7      0.       0        
GRID         318       0     0.3     0.6      0.       0        
GRID         319       0     0.2     0.8      0.       0        
GRID         320       0     0.2     0.7      0.       0        
GRID         321       0     0.4     0.8      0.       0        
GRID         322       0     0.3     0.8      0.       0        
GRID         323       0     0.4     0.4      0.       0        
GRID         324       0     0.5     0.5      0.       0        
GRID         325       0     0.5     0.4      0.       0        
GRID         327       0     0.6     0.4      0.       0        
GRID         329       0     0.8     0.5      0.       0        
GRID         330       0     0.7     0.5      0.       0        
GRID         331       0     0.7     0.4      0.       0        
GRID         332       0     0.6     0.6      0.       0        
GRID         333       0     0.6     0.5      0.       0        
GRID         334       0     0.8     0.6      0.       0        
GRID         335       0     0.4     0.6      0.       0        
GRID         336       0     0.5     0.7      0.       0        
GRID         337       0     0.5     0.6      0.       0        
GRID         339       0     0.4     0.7      0.       0        
GRID         340       0     0.5     0.8      0.       0        
GRID         341       0     0.8     0.7      0.       0        
GRID         342       0     0.7     0.7      0.       0        
GRID         343       0     0.7     0.6      0.       0        
GRID         344       0     0.6     0.8      0.       0        
GRID         345       0     0.6     0.7      0.       0        
GRID         346       0     0.8     0.8      0.       0        
GRID         347       0     0.7     0.8      0.       0        
CQUAD4       201       1     239     238     241     235                
CQUAD4       202       1     246     252     238     239                
CQUAD4       203       1     238     259     255     241                
CQUAD4       204       1     252     251     259     238                
CQUAD4       205       1     250     249     252     246                
CQUAD4       206       1     273     248     249     250                
CQUAD4       207       1     249     268     251     252                
CQUAD4       208       1     248     285     268     249                
CQUAD4       209       1     259     258     261     255                
CQUAD4       210       1     251     270     258     259                
CQUAD4       211       1     258     263     298     261                
CQUAD4       212       1     270     302     263     258                
CQUAD4       213       1     268     267     270     251                
CQUAD4       214       1     285     266     267     268                
CQUAD4       215       1     267     306     302     270                
CQUAD4       216       1     266     323     306     267                
CQUAD4       217       1     275     274     248     273                
CQUAD4       218       1     277     283     274     275                
CQUAD4       219       1     274     287     285     248                
CQUAD4       220       1     283     282     287     274                
CQUAD4       221       1     281     280     283     277                
CQUAD4       222       1     278     279     280     281                
CQUAD4       223       1     280     293     282     283                
CQUAD4       224       1     279     284     293     280                
CQUAD4       225       1     287     286     266     285                
CQUAD4       226       1     282     295     286     287                
CQUAD4       227       1     286     325     323     266                
CQUAD4       228       1     295     327     325     286                
CQUAD4       229       1     293     292     295     282                
CQUAD4       230       1     284     291     292     293                
CQUAD4       231       1     292     331     327     295                
CQUAD4       232       1     291     296     331     292                
CQUAD4       233       1     263     299     301     298                
CQUAD4       234       1     302     308     299     263                
CQUAD4       235       1     299     312     310     301                
CQUAD4       236       1     308     307     312     299                
CQUAD4       237       1     306     305     308     302                
CQUAD4       238       1     323     304     305     306                
CQUAD4       239       1     305     318     307     308                
CQUAD4       240       1     304     335     318     305                
CQUAD4       241       1     312     311     314     310                
CQUAD4       242       1     307     320     311     312                
CQUAD4       243       1     311     315     313     314                
CQUAD4       244       1     320     319     315     311                
CQUAD4       245       1     318     317     320     307                
CQUAD4       246       1     335     339     317     318                
CQUAD4       247       1     317     322     319     320                
CQUAD4       248       1     339     321     322     317                
CQUAD4       249       1     325     324     304     323                
CQUAD4       250       1     327     333     324     325                
CQUAD4       251       1     324     337     335     304                
CQUAD4       252       1     333     332     337     324                
CQUAD4       253       1     331     330     333     327                
CQUAD4       254       1     296     329     330     331                
CQUAD4       255       1     330     343     332     333                
CQUAD4       256       1     329     334     343     330                
CQUAD4       257       1     337     336     339     335                
CQUAD4       258       1     332     345     336     337                
CQUAD4       259       1     336     340     321     339                
CQUAD4       260       1     345     344     340     336                
CQUAD4       261       1     343     342     345     332                
CQUAD4       262       1     334     341     342     343                
CQUAD4       263       1     342     347     344     345                
CQUAD4       264       1     341     346     347     342                
ENDDATA
