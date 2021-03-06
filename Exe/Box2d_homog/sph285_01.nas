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
$   Date       : Tue Sep 26 21:12:37 2006
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
$ FEMAP Load Set 111 : Inclusion BSOURCE 0.4806 0.4806 0.285 0.285 0.
$ FEMAP Property 1 : Matrix
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Fulleren
PSHELL         2       1      0.       1               1              0.
$ FEMAP Material 1 : Untitled
MAT1           1     7.8     2.3     0.3      0.      0.      0.        
GRID         363       0      0.  0.9612      0.       0        
GRID         364       0      0. 0.82389      0.       0        
GRID         365       0      0. 0.68657      0.       0        
GRID         366       0      0. 0.54926      0.       0        
GRID         367       0      0. 0.41194      0.       0        
GRID         368       0      0. 0.27463      0.       0        
GRID         369       0      0. 0.13731      0.       0        
GRID         370       0      0.      0.      0.       0        
GRID         371       0 0.13731      0.      0.       0        
GRID         372       0 0.27463      0.      0.       0        
GRID         373       0 0.41194      0.      0.       0        
GRID         374       0 0.54926      0.      0.       0        
GRID         375       0 0.68657      0.      0.       0        
GRID         376       0 0.82389      0.      0.       0        
GRID         377       0  0.9612      0.      0.       0        
GRID         378       0  0.9612 0.13731      0.       0        
GRID         379       0  0.9612 0.27463      0.       0        
GRID         380       0  0.9612 0.41194      0.       0        
GRID         381       0  0.9612 0.54926      0.       0        
GRID         382       0  0.9612 0.68657      0.       0        
GRID         383       0  0.9612 0.82389      0.       0        
GRID         384       0  0.9612  0.9612      0.       0        
GRID         385       0 0.82389  0.9612      0.       0        
GRID         386       0 0.68657  0.9612      0.       0        
GRID         387       0 0.54926  0.9612      0.       0        
GRID         388       0 0.41194  0.9612      0.       0        
GRID         389       0 0.27463  0.9612      0.       0        
GRID         390       0 0.13731  0.9612      0.       0        
GRID         392       0 0.23378  0.6231      0.       0        
GRID         395       0  0.6231 0.72742      0.       0        
GRID         396       0 0.72742  0.6231      0.       0        
GRID         398       0 0.72742  0.3381      0.       0        
GRID         401       0  0.3381 0.23378      0.       0        
GRID         402       0 0.23378  0.3381      0.       0        
GRID         403       0 0.12436 0.63706      0.       0        
GRID         404       0 0.85064 0.23382      0.       0        
GRID         405       0 0.23382 0.11056      0.       0        
GRID         406       00.091792 0.12313      0.       0        
GRID         407       0 0.17952 0.21607      0.       0        
GRID         408       0  0.3814    0.09      0.       0        
GRID         409       0 0.838070.091792      0.       0        
GRID         410       0 0.74513 0.17952      0.       0        
GRID         411       0  0.7219 0.86686      0.       0        
GRID         412       0  0.8712  0.3814      0.       0        
GRID         413       0 0.86828 0.84141      0.       0        
GRID         414       0 0.78036 0.74902      0.       0        
GRID         415       0 0.31559 0.83266      0.       0        
GRID         416       0 0.21463 0.73539      0.       0        
GRID         417       0 0.17404 0.84373      0.       0        
GRID         418       0  0.1025 0.74571      0.       0        
GRID         419       0  0.7656  0.4806      0.       0        
GRID         422       0  0.4806  0.7656      0.       0        
GRID         423       0  0.3381 0.72742      0.       0        
GRID         425       0  0.1956  0.4806      0.       0        
GRID         428       0  0.4806  0.1956      0.       0        
GRID         429       0  0.6231 0.23378      0.       0        
GRID         431       0 0.32493  0.4806      0.       0        
GRID         432       0  0.4806  0.4806      0.       0        
GRID         433       0 0.63627  0.4806      0.       0        
GRID         434       0 0.55843 0.61541      0.       0        
GRID         435       0 0.40277 0.61541      0.       0        
GRID         436       0 0.40277 0.34579      0.       0        
GRID         437       0 0.55843 0.34579      0.       0        
CTRIA3       552       1     403     365     366                        
CTRIA3       553       1     425     392     403                        
CTRIA3       554       1     425     403     366                        
CTRIA3       555       1     425     366     367                        
CTRIA3       556       1     371     405     406                        
CTRIA3       557       1     370     371     406                        
CTRIA3       558       1     369     370     406                        
CTRIA3       559       1     368     369     406                        
CTRIA3       560       1     368     406     407                        
CTRIA3       561       1     402     425     367                        
CTRIA3       562       1     402     367     368                        
CTRIA3       563       1     402     368     407                        
CTRIA3       564       1     407     406     405                        
CTRIA3       565       1     401     402     407                        
CTRIA3       566       1     401     407     405                        
CTRIA3       567       1     405     371     372                        
CTRIA3       568       1     372     373     408                        
CTRIA3       569       1     405     372     408                        
CTRIA3       570       1     401     405     408                        
CTRIA3       571       1     428     401     408                        
CTRIA3       572       1     408     373     374                        
CTRIA3       573       1     428     408     374                        
CTRIA3       574       1     378     404     409                        
CTRIA3       575       1     377     378     409                        
CTRIA3       576       1     376     377     409                        
CTRIA3       577       1     375     376     409                        
CTRIA3       578       1     375     409     410                        
CTRIA3       579       1     429     428     374                        
CTRIA3       580       1     429     374     375                        
CTRIA3       581       1     429     375     410                        
CTRIA3       582       1     410     409     404                        
CTRIA3       583       1     398     429     410                        
CTRIA3       584       1     398     410     404                        
CTRIA3       585       1     404     378     379                        
CTRIA3       586       1     379     380     412                        
CTRIA3       587       1     404     379     412                        
CTRIA3       588       1     398     404     412                        
CTRIA3       589       1     419     398     412                        
CTRIA3       590       1     412     380     381                        
CTRIA3       591       1     419     412     381                        
CTRIA3       592       1     385     411     413                        
CTRIA3       593       1     384     385     413                        
CTRIA3       594       1     383     384     413                        
CTRIA3       595       1     382     383     413                        
CTRIA3       596       1     382     413     414                        
CTRIA3       597       1     396     419     381                        
CTRIA3       598       1     396     381     382                        
CTRIA3       599       1     396     382     414                        
CTRIA3       600       1     414     413     411                        
CTRIA3       601       1     395     396     414                        
CTRIA3       602       1     395     414     411                        
CTRIA3       603       1     411     385     386                        
CTRIA3       604       1     411     386     387                        
CTRIA3       605       1     395     411     387                        
CTRIA3       606       1     422     395     387                        
CTRIA3       607       1     422     387     388                        
CTRIA3       608       1     422     388     415                        
CTRIA3       609       1     423     422     415                        
CTRIA3       610       1     423     415     416                        
CTRIA3       611       1     392     423     416                        
CTRIA3       612       1     403     392     416                        
CTRIA3       613       1     415     388     389                        
CTRIA3       614       1     389     390     417                        
CTRIA3       615       1     415     389     417                        
CTRIA3       616       1     416     415     417                        
CTRIA3       617       1     403     416     418                        
CTRIA3       618       1     365     403     418                        
CTRIA3       619       1     364     365     418                        
CTRIA3       620       1     418     416     417                        
CTRIA3       621       1     364     418     417                        
CTRIA3       622       1     364     417     390                        
CTRIA3       623       1     363     364     390                        
CTRIA3       624       2     431     432     435                        
CTRIA3       625       2     392     425     431                        
CTRIA3       626       2     392     431     435                        
CTRIA3       627       2     423     392     435                        
CTRIA3       628       2     435     432     434                        
CTRIA3       629       2     422     423     435                        
CTRIA3       630       2     422     435     434                        
CTRIA3       631       2     395     422     434                        
CTRIA3       632       2     434     432     433                        
CTRIA3       633       2     396     395     434                        
CTRIA3       634       2     396     434     433                        
CTRIA3       635       2     419     396     433                        
CTRIA3       636       2     402     401     436                        
CTRIA3       637       2     431     425     402                        
CTRIA3       638       2     431     402     436                        
CTRIA3       639       2     432     431     436                        
CTRIA3       640       2     436     401     428                        
CTRIA3       641       2     428     429     437                        
CTRIA3       642       2     436     428     437                        
CTRIA3       643       2     432     436     437                        
CTRIA3       644       2     433     432     437                        
CTRIA3       645       2     437     429     398                        
CTRIA3       646       2     433     437     398                        
CTRIA3       647       2     419     433     398                        
ENDDATA
