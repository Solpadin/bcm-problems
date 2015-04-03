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
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\Box2d_homog\Solid\model.mod
$   Date       : Sat Nov 05 12:54:46 2005
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
$ FEMAP Load Set 111 : Inclusion BSOURCE 0.5 0.5 0.25 0.25 0.
$ FEMAP Property 1 : Matrix
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Fulleren
PSHELL         2       1      0.       1               1              0.
$ FEMAP Material 1 : Untitled
MAT1           1     7.8     2.3     0.3      0.      0.      0.        
GRID        2402       0   0.375     0.5      0.       0        
GRID        2404       0     0.5   0.625      0.       0        
GRID        2406       0   0.375 0.71651      0.       0        
GRID        2408       0 0.42225 0.61942      0.       0        
GRID        2409       0      0.      1.      0.       0        
GRID        2410       0      0. 0.83333      0.       0        
GRID        2411       0      0. 0.66667      0.       0        
GRID        2413       0   0.125     0.5      0.       0        
GRID        2415       0 0.28349   0.625      0.       0        
GRID        2419       0     0.5      1.      0.       0        
GRID        2420       0 0.33333      1.      0.       0        
GRID        2421       0 0.16667      1.      0.       0        
GRID        2422       0 0.34955 0.84548      0.       0        
GRID        2423       0 0.18989 0.78984      0.       0        
GRID        2424       0 0.14959 0.64552      0.       0        
GRID        2425       0    0.25     0.5      0.       0        
GRID        2427       0     0.5     0.5      0.       0        
GRID        2430       0   0.375 0.28349      0.       0        
GRID        2432       0 0.42225 0.38058      0.       0        
GRID        2433       0      0.      0.      0.       0        
GRID        2434       0      0. 0.16667      0.       0        
GRID        2435       0      0. 0.33333      0.       0        
GRID        2436       0      0.     0.5      0.       0        
GRID        2439       0 0.28349   0.375      0.       0        
GRID        2442       0     0.5   0.125      0.       0        
GRID        2444       0 0.33333      0.      0.       0        
GRID        2445       0 0.16667      0.      0.       0        
GRID        2446       0 0.34955 0.15452      0.       0        
GRID        2447       0 0.18989 0.21016      0.       0        
GRID        2448       0 0.14959 0.35448      0.       0        
GRID        2450       0   0.625     0.5      0.       0        
GRID        2453       0     0.5    0.75      0.       0        
GRID        2454       0   0.625 0.71651      0.       0        
GRID        2456       0 0.57775 0.61942      0.       0        
GRID        2457       0      1.      1.      0.       0        
GRID        2458       0      1. 0.83333      0.       0        
GRID        2459       0      1. 0.66667      0.       0        
GRID        2460       0      1.     0.5      0.       0        
GRID        2463       0 0.71651   0.625      0.       0        
GRID        2466       0     0.5   0.875      0.       0        
GRID        2468       0 0.66667      1.      0.       0        
GRID        2469       0 0.83333      1.      0.       0        
GRID        2470       0 0.65045 0.84548      0.       0        
GRID        2471       0 0.81011 0.78984      0.       0        
GRID        2472       0 0.85041 0.64552      0.       0        
GRID        2476       0     0.5   0.375      0.       0        
GRID        2477       0     0.5    0.25      0.       0        
GRID        2478       0   0.625 0.28349      0.       0        
GRID        2480       0 0.57775 0.38058      0.       0        
GRID        2481       0      1.      0.      0.       0        
GRID        2482       0      1. 0.16667      0.       0        
GRID        2483       0      1. 0.33333      0.       0        
GRID        2485       0   0.875     0.5      0.       0        
GRID        2486       0    0.75     0.5      0.       0        
GRID        2487       0 0.71651   0.375      0.       0        
GRID        2491       0     0.5      0.      0.       0        
GRID        2492       0 0.66667      0.      0.       0        
GRID        2493       0 0.83333      0.      0.       0        
GRID        2494       0 0.65045 0.15452      0.       0        
GRID        2495       0 0.81011 0.21016      0.       0        
GRID        2496       0 0.85041 0.35448      0.       0        
CTRIA3      3439       2    2406    2415    2408                        
CTRIA3      3440       2    2453    2406    2408                        
CTRIA3      3441       2    2404    2453    2408                        
CTRIA3      3442       2    2427    2404    2408                        
CTRIA3      3443       2    2402    2427    2408                        
CTRIA3      3444       2    2402    2408    2415                        
CTRIA3      3445       2    2425    2402    2415                        
CQUAD4      3446       1    2413    2425    2415    2424                
CQUAD4      3447       1    2411    2436    2413    2424                
CQUAD4      3448       1    2406    2453    2466    2422                
CQUAD4      3449       1    2415    2406    2422    2423                
CTRIA3      3450       1    2424    2415    2423                        
CQUAD4      3451       1    2410    2411    2424    2423                
CQUAD4      3452       1    2422    2466    2419    2420                
CQUAD4      3453       1    2423    2422    2420    2421                
CQUAD4      3454       1    2409    2410    2423    2421                
CTRIA3      3455       2    2439    2430    2432                        
CTRIA3      3456       2    2430    2477    2432                        
CTRIA3      3457       2    2477    2476    2432                        
CTRIA3      3458       2    2476    2427    2432                        
CTRIA3      3459       2    2427    2402    2432                        
CTRIA3      3460       2    2432    2402    2439                        
CTRIA3      3461       2    2402    2425    2439                        
CQUAD4      3462       1    2425    2413    2448    2439                
CQUAD4      3463       1    2436    2435    2448    2413                
CQUAD4      3464       1    2477    2430    2446    2442                
CQUAD4      3465       1    2430    2439    2447    2446                
CTRIA3      3466       1    2439    2448    2447                        
CQUAD4      3467       1    2435    2434    2447    2448                
CQUAD4      3468       1    2442    2446    2444    2491                
CQUAD4      3469       1    2446    2447    2445    2444                
CQUAD4      3470       1    2434    2433    2445    2447                
CTRIA3      3471       2    2463    2454    2456                        
CTRIA3      3472       2    2454    2453    2456                        
CTRIA3      3473       2    2453    2404    2456                        
CTRIA3      3474       2    2404    2427    2456                        
CTRIA3      3475       2    2427    2450    2456                        
CTRIA3      3476       2    2456    2450    2463                        
CTRIA3      3477       2    2450    2486    2463                        
CQUAD4      3478       1    2486    2485    2472    2463                
CQUAD4      3479       1    2460    2459    2472    2485                
CQUAD4      3480       1    2453    2454    2470    2466                
CQUAD4      3481       1    2454    2463    2471    2470                
CTRIA3      3482       1    2463    2472    2471                        
CQUAD4      3483       1    2459    2458    2471    2472                
CQUAD4      3484       1    2466    2470    2468    2419                
CQUAD4      3485       1    2470    2471    2469    2468                
CQUAD4      3486       1    2458    2457    2469    2471                
CTRIA3      3487       2    2478    2487    2480                        
CTRIA3      3488       2    2477    2478    2480                        
CTRIA3      3489       2    2476    2477    2480                        
CTRIA3      3490       2    2427    2476    2480                        
CTRIA3      3491       2    2450    2427    2480                        
CTRIA3      3492       2    2450    2480    2487                        
CTRIA3      3493       2    2486    2450    2487                        
CQUAD4      3494       1    2485    2486    2487    2496                
CQUAD4      3495       1    2483    2460    2485    2496                
CQUAD4      3496       1    2478    2477    2442    2494                
CQUAD4      3497       1    2487    2478    2494    2495                
CTRIA3      3498       1    2496    2487    2495                        
CQUAD4      3499       1    2482    2483    2496    2495                
CQUAD4      3500       1    2494    2442    2491    2492                
CQUAD4      3501       1    2495    2494    2492    2493                
CQUAD4      3502       1    2481    2482    2495    2493                
ENDDATA
