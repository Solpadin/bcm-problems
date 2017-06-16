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
$   Date       : Wed Jul 16 15:52:03 2008
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
$ FEMAP Load Set 111 : Cylinder BSOURCE 0.015
$ FEMAP Property 1 : Untitled
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Untitled
PSHELL         2       1      0.       1               1              0.
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID         369       0      0.    0.05      0.       0        
GRID         370       0      0.      0.      0.       0        
GRID         371       0    0.05      0.      0.       0        
GRID         372       0    0.05    0.05      0.       0        
GRID         454       0     0.2      0.      0.       0        
GRID         455       0     0.2    0.05      0.       0        
GRID         457       0     0.1    0.05      0.       0        
GRID         458       0     0.1      0.      0.       0        
GRID         459       0    0.15      0.      0.       0        
GRID         460       0    0.15    0.05      0.       0        
GRID         518       0     0.3    0.05      0.       0        
GRID         522       0    0.25      0.      0.       0        
GRID         523       0    0.25    0.05      0.       0        
GRID         580       0     0.4      0.      0.       0        
GRID         581       0     0.4    0.05      0.       0        
GRID         584       0     0.3      0.      0.       0        
GRID         585       0    0.35      0.      0.       0        
GRID         586       0    0.35    0.05      0.       0        
GRID         644       0     0.5    0.05      0.       0        
GRID         648       0    0.45      0.      0.       0        
GRID         649       0    0.45    0.05      0.       0        
GRID         710       0     0.5      0.      0.       0        
GRID         711       0    0.55      0.      0.       0        
GRID         712       0    0.55    0.05      0.       0        
GRID         772       0     0.6    0.05      0.       0        
GRID         773       0     0.6      0.      0.       0        
GRID         774       0    0.65      0.      0.       0        
GRID         775       0    0.65    0.05      0.       0        
GRID         832       0     0.8      0.      0.       0        
GRID         835       0     0.7    0.05      0.       0        
GRID         836       0     0.7      0.      0.       0        
GRID         837       0    0.75      0.      0.       0        
GRID         838       0    0.75    0.05      0.       0        
GRID         898       0     0.8    0.05      0.       0        
GRID         900       0    0.85      0.      0.       0        
GRID         901       0    0.85    0.05      0.       0        
GRID         959       0      1.    0.05      0.       0        
GRID         961       0     0.9    0.05      0.       0        
GRID         962       0     0.9      0.      0.       0        
GRID         963       0    0.95      0.      0.       0        
GRID         964       0    0.95    0.05      0.       0        
GRID        1021       0     1.1      0.      0.       0        
GRID        1025       0      1.      0.      0.       0        
GRID        1026       0    1.05      0.      0.       0        
GRID        1027       0    1.05    0.05      0.       0        
GRID        1084       0     1.2      0.      0.       0        
GRID        1085       0     1.2    0.05      0.       0        
GRID        1087       0     1.1    0.05      0.       0        
GRID        1089       0    1.15      0.      0.       0        
GRID        1090       0    1.15    0.05      0.       0        
GRID        1147       0     1.3      0.      0.       0        
GRID        1148       0     1.3    0.05      0.       0        
GRID        1152       0    1.25      0.      0.       0        
GRID        1153       0    1.25    0.05      0.       0        
GRID        1211       0     1.4    0.05      0.       0        
GRID        1215       0    1.35      0.      0.       0        
GRID        1216       0    1.35    0.05      0.       0        
GRID        1277       0     1.4      0.      0.       0        
GRID        1278       0    1.45      0.      0.       0        
GRID        1279       0    1.45    0.05      0.       0        
GRID        1336       0     1.6      0.      0.       0        
GRID        1337       0     1.6    0.05      0.       0        
GRID        1339       0     1.5    0.05      0.       0        
GRID        1340       0     1.5      0.      0.       0        
GRID        1341       0    1.55      0.      0.       0        
GRID        1342       0    1.55    0.05      0.       0        
GRID        1399       0     1.7      0.      0.       0        
GRID        1400       0     1.7    0.05      0.       0        
GRID        1404       0    1.65      0.      0.       0        
GRID        1405       0    1.65    0.05      0.       0        
GRID        1462       0     1.8      0.      0.       0        
GRID        1467       0    1.75      0.      0.       0        
GRID        1468       0    1.75    0.05      0.       0        
GRID        1526       0     1.9    0.05      0.       0        
GRID        1528       0     1.8    0.05      0.       0        
GRID        1530       0    1.85      0.      0.       0        
GRID        1531       0    1.85    0.05      0.       0        
GRID        1588       0      2.      0.      0.       0        
GRID        1589       0      2.    0.05      0.       0        
GRID        1592       0     1.9      0.      0.       0        
GRID        1593       0    1.95      0.      0.       0        
GRID        1594       0    1.95    0.05      0.       0        
CQUAD4       252       1     458     457     372     371                
CQUAD4       254       1     371     372     369     370                
CQUAD4       292       1     454     455     460     459                
CQUAD4       294       1     459     460     457     458                
CQUAD4       332       1     584     518     523     522                
CQUAD4       334       1     522     523     455     454                
CQUAD4       372       1     580     581     586     585                
CQUAD4       374       1     585     586     518     584                
CQUAD4       412       1     710     644     649     648                
CQUAD4       414       1     648     649     581     580                
CQUAD4       452       1     773     772     712     711                
CQUAD4       454       1     711     712     644     710                
CQUAD4       492       1     836     835     775     774                
CQUAD4       494       1     774     775     772     773                
CQUAD4       532       1     832     898     838     837                
CQUAD4       534       1     837     838     835     836                
CQUAD4       572       1     962     961     901     900                
CQUAD4       574       1     900     901     898     832                
CQUAD4       612       1    1025     959     964     963                
CQUAD4       614       1     963     964     961     962                
CQUAD4       652       1    1021    1087    1027    1026                
CQUAD4       654       1    1026    1027     959    1025                
CQUAD4       692       1    1084    1085    1090    1089                
CQUAD4       694       1    1089    1090    1087    1021                
CQUAD4       732       1    1147    1148    1153    1152                
CQUAD4       734       1    1152    1153    1085    1084                
CQUAD4       772       1    1277    1211    1216    1215                
CQUAD4       774       1    1215    1216    1148    1147                
CQUAD4       812       1    1340    1339    1279    1278                
CQUAD4       814       1    1278    1279    1211    1277                
CQUAD4       852       1    1336    1337    1342    1341                
CQUAD4       854       1    1341    1342    1339    1340                
CQUAD4       892       1    1399    1400    1405    1404                
CQUAD4       894       1    1404    1405    1337    1336                
CQUAD4       932       1    1462    1528    1468    1467                
CQUAD4       934       1    1467    1468    1400    1399                
CQUAD4       972       1    1592    1526    1531    1530                
CQUAD4       974       1    1530    1531    1528    1462                
CQUAD4      1012       1    1588    1589    1594    1593                
CQUAD4      1014       1    1593    1594    1526    1592                
ENDDATA
