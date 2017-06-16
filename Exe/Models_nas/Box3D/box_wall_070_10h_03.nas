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
  LOAD = 2
BEGIN BULK
$ ***************************************************************************
$   Written by : FEMAP
$   Version    : 7.00
$   Translator : MSC/NASTRAN
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\Models_nas\Box3D\Solid\Box_model.MOD
$   Date       : Wed Feb 11 07:01:52 2009
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
$ FEMAP Load Set 1 : Source Vn -.001 0.0
PLOAD4         1     137      1.                             214     218
PLOAD4         1     138      1.                             215     234
$ FEMAP Load Set 2 : Boundary Skews
PLOAD4         2     139      2.                             236     299
PLOAD4         2     140      2.                             226     268
PLOAD4         2     139      2.                             254     261
PLOAD4         2     140      2.                             264     291
PLOAD4         2     139      2.                             236     312
PLOAD4         2     140      2.                             221     268
$ FEMAP Property 1 : Acoustic
PSOLID         1       1       0        
$ FEMAP Property 2 : Rigid body
PSOLID         2       1       0        
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID         213       0      0.      1.   -0.75       0        
GRID         214       0      0.      0.   -0.75       0        
GRID         215       0      0.      0.      0.       0        
GRID         218       0      0.      1.      0.       0        
GRID         221       0     0.8      0.      0.       0        
GRID         222       0     0.8      0.    0.75       0        
GRID         223       0      0.      0.    0.75       0        
GRID         226       0     0.8      1.    0.75       0        
GRID         234       0      0.      1.    0.75       0        
GRID         236       0     0.8      0.   -0.75       0        
GRID         254       0     0.8      1.   -0.75       0        
GRID         261       0    0.83      1.      0.       0        
GRID         264       0     0.8      1.      0.       0        
GRID         266       0    0.83      0.   -0.75       0        
GRID         268       0    0.83      0.    0.75       0        
GRID         278       0      3.      1.   -0.75       0        
GRID         279       0      3.      1.      0.       0        
GRID         280       0      3.      1.    0.75       0        
GRID         282       0      3.      0.      0.       0        
GRID         289       0 2.27667      1.   -0.75       0        
GRID         291       0    0.83      1.    0.75       0        
GRID         293       0 1.55333      0.    0.75       0        
GRID         299       0    0.83      1.   -0.75       0        
GRID         302       0 1.55333      1.    0.75       0        
GRID         303       0 2.27667      1.    0.75       0        
GRID         308       0 1.55333      1.   -0.75       0        
GRID         309       0 1.55333      1.      0.       0        
GRID         310       0 2.27667      1.      0.       0        
GRID         312       0    0.83      0.      0.       0        
GRID         314       0 1.55333      0.   -0.75       0        
GRID         315       0 2.27667      0.   -0.75       0        
GRID         316       0      3.      0.   -0.75       0        
GRID         318       0      3.      0.    0.75       0        
GRID         319       0 2.27667      0.    0.75       0        
GRID         321       0 1.55333      0.      0.       0        
GRID         322       0 2.27667      0.      0.       0        
CHEXA        137       1     214     215     221     236     213     218+EL   3T
+EL   3T     264     254                                                        
CHEXA        138       1     215     223     222     221     218     234+EL   3U
+EL   3U     226     264                                                        
CHEXA        139       2     236     254     264     221     266     299+EL   3V
+EL   3V     261     312                                                        
CHEXA        140       2     221     264     226     222     312     261+EL   3W
+EL   3W     291     268                                                        
CHEXA        141       1     316     278     289     315     282     279+EL   3X
+EL   3X     310     322                                                        
CHEXA        142       1     315     289     308     314     322     310+EL   3Y
+EL   3Y     309     321                                                        
CHEXA        143       1     314     308     299     266     321     309+EL   3Z
+EL   3Z     261     312                                                        
CHEXA        144       1     282     279     310     322     318     280+EL   40
+EL   40     303     319                                                        
CHEXA        145       1     322     310     309     321     319     303+EL   41
+EL   41     302     293                                                        
CHEXA        146       1     321     309     261     312     293     302+EL   42
+EL   42     291     268                                                        
ENDDATA
