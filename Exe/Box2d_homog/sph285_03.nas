ID C:\Users\Dima\Lame3d2\IDE,FEMAP
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
$   From Model : C:\Users\Dima\Lame3d2\IDENT_12septmber2006\Ident\Exe\box2d_fulleren\Solid\model.MOD
$   Date       : Tue Sep 26 21:38:15 2006
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
$ FEMAP Load Set 111 : Inclusion BSOURCE 0.34355 0.34355 0.285 0.285 0.
$ FEMAP Property 1 : Matrix
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Fulleren
PSHELL         2       1      0.       1               1              0.
$ FEMAP Material 1 : Untitled
MAT1           1     7.8     2.3     0.3      0.      0.      0.        
GRID         224       0      0.  0.6871      0.       0        
GRID         225       0      0. 0.58894      0.       0        
GRID         226       0      0. 0.49079      0.       0        
GRID         227       0      0. 0.39263      0.       0        
GRID         228       0      0. 0.29447      0.       0        
GRID         229       0      0. 0.19631      0.       0        
GRID         230       0      0.0.098157      0.       0        
GRID         231       0      0.      0.      0.       0        
GRID         232       00.098157      0.      0.       0        
GRID         233       0 0.19631      0.      0.       0        
GRID         234       0 0.29447      0.      0.       0        
GRID         235       0 0.39263      0.      0.       0        
GRID         236       0 0.49079      0.      0.       0        
GRID         237       0 0.58894      0.      0.       0        
GRID         238       0  0.6871      0.      0.       0        
GRID         239       0  0.68710.098157      0.       0        
GRID         240       0  0.6871 0.19631      0.       0        
GRID         241       0  0.6871 0.29447      0.       0        
GRID         242       0  0.6871 0.39263      0.       0        
GRID         243       0  0.6871 0.49079      0.       0        
GRID         244       0  0.6871 0.58894      0.       0        
GRID         245       0  0.6871  0.6871      0.       0        
GRID         246       0 0.58894  0.6871      0.       0        
GRID         247       0 0.49079  0.6871      0.       0        
GRID         248       0 0.39263  0.6871      0.       0        
GRID         249       0 0.29447  0.6871      0.       0        
GRID         250       0 0.19631  0.6871      0.       0        
GRID         251       00.098157  0.6871      0.       0        
GRID         254       0 0.12523 0.52674      0.       0        
GRID         256       0 0.29406 0.62422      0.       0        
GRID         257       0 0.39304 0.62422      0.       0        
GRID         259       0 0.56187 0.52674      0.       0        
GRID         260       0 0.61136 0.44103      0.       0        
GRID         264       0 0.486050.096733      0.       0        
GRID         268       0 0.12523 0.16036      0.       0        
GRID         270       0 0.583640.059208      0.       0        
GRID         271       0 0.58364 0.62789      0.       0        
GRID         272       0 0.103460.059208      0.       0        
GRID         273       0 0.10346 0.62789      0.       0        
GRID         274       0 0.62855 0.34355      0.       0        
GRID         277       0 0.48605 0.59037      0.       0        
GRID         280       0 0.20105 0.59037      0.       0        
GRID         282       00.075738 0.44103      0.       0        
GRID         283       0 0.05855 0.34355      0.       0        
GRID         284       00.075738 0.24607      0.       0        
GRID         286       0 0.201050.096733      0.       0        
GRID         287       0 0.29406 0.06288      0.       0        
GRID         288       0 0.39304 0.06288      0.       0        
GRID         290       0 0.56187 0.16036      0.       0        
GRID         291       0 0.61136 0.24607      0.       0        
GRID         292       0 0.44686 0.16461      0.       0        
GRID         293       0 0.39565 0.25331      0.       0        
GRID         294       0 0.34355 0.34355      0.       0        
GRID         295       0 0.29145 0.43379      0.       0        
GRID         296       0 0.24024 0.52249      0.       0        
GRID         297       0 0.29145 0.25331      0.       0        
GRID         298       0 0.24024 0.16461      0.       0        
GRID         299       0 0.34355 0.16063      0.       0        
GRID         300       0 0.23935 0.34355      0.       0        
GRID         301       0 0.13693 0.34355      0.       0        
GRID         302       0 0.18514 0.25209      0.       0        
GRID         303       0 0.18514 0.43501      0.       0        
GRID         304       0 0.39565 0.43379      0.       0        
GRID         305       0 0.44686 0.52249      0.       0        
GRID         306       0 0.34355 0.52647      0.       0        
GRID         307       0 0.44775 0.34355      0.       0        
GRID         308       0 0.55017 0.34355      0.       0        
GRID         309       0 0.50196 0.43501      0.       0        
GRID         310       0 0.50196 0.25209      0.       0        
CTRIA3       370       1     288     235     236                        
CTRIA3       371       1     264     288     236                        
CTRIA3       372       1     236     237     270                        
CTRIA3       373       1     264     236     270                        
CTRIA3       374       1     290     264     270                        
CTRIA3       375       1     270     237     238                        
CTRIA3       376       1     270     238     239                        
CTRIA3       377       1     290     270     239                        
CTRIA3       378       1     290     239     240                        
CTRIA3       379       1     291     290     240                        
CTRIA3       380       1     291     240     241                        
CTRIA3       381       1     274     291     241                        
CTRIA3       382       1     274     241     242                        
CTRIA3       383       1     260     274     242                        
CTRIA3       384       1     260     242     243                        
CTRIA3       385       1     259     260     243                        
CTRIA3       386       1     259     243     244                        
CTRIA3       387       1     245     246     271                        
CTRIA3       388       1     244     245     271                        
CTRIA3       389       1     259     244     271                        
CTRIA3       390       1     277     259     271                        
CTRIA3       391       1     271     246     247                        
CTRIA3       392       1     277     271     247                        
CTRIA3       393       1     257     277     247                        
CTRIA3       394       1     257     247     248                        
CTRIA3       395       1     256     257     248                        
CTRIA3       396       1     256     248     249                        
CTRIA3       397       1     256     249     250                        
CTRIA3       398       1     280     256     250                        
CTRIA3       399       1     234     235     288                        
CTRIA3       400       1     234     288     287                        
CTRIA3       401       1     233     234     287                        
CTRIA3       402       1     233     287     286                        
CTRIA3       403       1     286     268     272                        
CTRIA3       404       1     233     286     272                        
CTRIA3       405       1     232     233     272                        
CTRIA3       406       1     231     232     272                        
CTRIA3       407       1     230     231     272                        
CTRIA3       408       1     230     272     268                        
CTRIA3       409       1     229     230     268                        
CTRIA3       410       1     229     268     284                        
CTRIA3       411       1     228     229     284                        
CTRIA3       412       1     228     284     283                        
CTRIA3       413       1     227     228     283                        
CTRIA3       414       1     227     283     282                        
CTRIA3       415       1     226     227     282                        
CTRIA3       416       1     226     282     254                        
CTRIA3       417       1     225     226     254                        
CTRIA3       418       1     250     251     273                        
CTRIA3       419       1     280     250     273                        
CTRIA3       420       1     254     280     273                        
CTRIA3       421       1     225     254     273                        
CTRIA3       422       1     224     225     273                        
CTRIA3       423       1     224     273     251                        
CTRIA3       424       2     297     298     299                        
CTRIA3       425       2     293     294     297                        
CTRIA3       426       2     293     297     299                        
CTRIA3       427       2     292     293     299                        
CTRIA3       428       2     288     264     292                        
CTRIA3       429       2     288     292     299                        
CTRIA3       430       2     287     288     299                        
CTRIA3       431       2     287     299     298                        
CTRIA3       432       2     286     287     298                        
CTRIA3       433       2     297     294     300                        
CTRIA3       434       2     300     301     302                        
CTRIA3       435       2     297     300     302                        
CTRIA3       436       2     298     297     302                        
CTRIA3       437       2     268     286     298                        
CTRIA3       438       2     268     298     302                        
CTRIA3       439       2     284     268     302                        
CTRIA3       440       2     284     302     301                        
CTRIA3       441       2     283     284     301                        
CTRIA3       442       2     295     296     303                        
CTRIA3       443       2     300     294     295                        
CTRIA3       444       2     300     295     303                        
CTRIA3       445       2     301     300     303                        
CTRIA3       446       2     282     283     301                        
CTRIA3       447       2     282     301     303                        
CTRIA3       448       2     254     282     303                        
CTRIA3       449       2     254     303     296                        
CTRIA3       450       2     280     254     296                        
CTRIA3       451       2     295     294     304                        
CTRIA3       452       2     304     305     306                        
CTRIA3       453       2     295     304     306                        
CTRIA3       454       2     296     295     306                        
CTRIA3       455       2     256     280     296                        
CTRIA3       456       2     256     296     306                        
CTRIA3       457       2     257     256     306                        
CTRIA3       458       2     257     306     305                        
CTRIA3       459       2     277     257     305                        
CTRIA3       460       2     307     308     309                        
CTRIA3       461       2     304     294     307                        
CTRIA3       462       2     304     307     309                        
CTRIA3       463       2     305     304     309                        
CTRIA3       464       2     259     277     305                        
CTRIA3       465       2     259     305     309                        
CTRIA3       466       2     260     259     309                        
CTRIA3       467       2     260     309     308                        
CTRIA3       468       2     274     260     308                        
CTRIA3       469       2     307     294     293                        
CTRIA3       470       2     293     292     310                        
CTRIA3       471       2     307     293     310                        
CTRIA3       472       2     308     307     310                        
CTRIA3       473       2     292     264     290                        
CTRIA3       474       2     310     292     290                        
CTRIA3       475       2     310     290     291                        
CTRIA3       476       2     308     310     291                        
CTRIA3       477       2     274     308     291                        
ENDDATA
