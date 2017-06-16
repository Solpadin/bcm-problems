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
$   Date       : Mon Jun 16 12:17:31 2008
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
GRID         247       0     0.4      0.      0.       0        
GRID         248       0     0.4     0.1      0.       0        
GRID         249       0     0.3     0.1      0.       0        
GRID         250       0     0.3      0.      0.       0        
GRID         251       0     0.2     0.2      0.       0        
GRID         252       0     0.2     0.1      0.       0        
GRID         253       0     0.4     0.2      0.       0        
GRID         255       0      0.     0.2      0.       0        
GRID         258       0     0.1     0.3      0.       0        
GRID         259       0     0.1     0.2      0.       0        
GRID         260       0      0.     0.4      0.       0        
GRID         261       0      0.     0.3      0.       0        
GRID         263       0     0.1     0.4      0.       0        
GRID         266       0     0.4     0.3      0.       0        
GRID         267       0     0.3     0.3      0.       0        
GRID         268       0     0.3     0.2      0.       0        
GRID         269       0     0.2     0.4      0.       0        
GRID         270       0     0.2     0.3      0.       0        
GRID         271       0     0.4     0.4      0.       0        
GRID         272       0     0.3     0.4      0.       0        
CQUAD4       201       1     239     238     241     235                
CQUAD4       202       1     246     252     238     239                
CQUAD4       203       1     238     259     255     241                
CQUAD4       204       1     252     251     259     238                
CQUAD4       205       1     250     249     252     246                
CQUAD4       206       1     247     248     249     250                
CQUAD4       207       1     249     268     251     252                
CQUAD4       208       1     248     253     268     249                
CQUAD4       209       1     259     258     261     255                
CQUAD4       210       1     251     270     258     259                
CQUAD4       211       1     258     263     260     261                
CQUAD4       212       1     270     269     263     258                
CQUAD4       213       1     268     267     270     251                
CQUAD4       214       1     253     266     267     268                
CQUAD4       215       1     267     272     269     270                
CQUAD4       216       1     266     271     272     267                
ENDDATA
