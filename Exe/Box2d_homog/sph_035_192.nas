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
$   Date       : Sat Nov 05 11:47:40 2005
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
$ FEMAP Load Set 111 : Inclusion BSOURCE 0.5 0.5 0.35 0.35 0.
$ FEMAP Property 1 : Matrix
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Fulleren
PSHELL         2       1      0.       1               1              0.
$ FEMAP Material 1 : Untitled
MAT1           1     7.8     2.3     0.3      0.      0.      0.        
GRID        2507       0      0.      1.      0.       0        
GRID        2508       0      0.     0.9      0.       0        
GRID        2509       0      0.     0.8      0.       0        
GRID        2510       0      0.     0.7      0.       0        
GRID        2511       0      0.     0.6      0.       0        
GRID        2516       0 0.29428 0.78316      0.       0        
GRID        2519       0     0.5      1.      0.       0        
GRID        2520       0     0.4      1.      0.       0        
GRID        2521       0     0.3      1.      0.       0        
GRID        2522       0     0.2      1.      0.       0        
GRID        2523       0     0.1      1.      0.       0        
GRID        2524       0 0.13338 0.90826      0.       0        
GRID        2525       0 0.17722 0.81387      0.       0        
GRID        2526       0 0.24307 0.90015      0.       0        
GRID        2527       0 0.35483  0.8944      0.       0        
GRID        2528       00.096532 0.74983      0.       0        
GRID        2529       0 0.10498 0.64403      0.       0        
GRID        2530       00.081066  0.8347      0.       0        
GRID        2536       0     0.5 0.73333      0.       0        
GRID        2538       0 0.39184 0.83287      0.       0        
GRID        2540       0 0.21684 0.70572      0.       0        
GRID        2541       0 0.16713 0.60816      0.       0        
GRID        2542       0 0.32229 0.68724      0.       0        
GRID        2543       0 0.35582 0.59684      0.       0        
GRID        2544       0 0.41644 0.67867      0.       0        
GRID        2545       0 0.40411 0.76073      0.       0        
GRID        2546       0 0.24645 0.59947      0.       0        
GRID        2547       0      0.      0.      0.       0        
GRID        2548       0      0.     0.1      0.       0        
GRID        2549       0      0.     0.2      0.       0        
GRID        2550       0      0.     0.3      0.       0        
GRID        2551       0      0.     0.4      0.       0        
GRID        2552       0      0.     0.5      0.       0        
GRID        2553       0    0.15     0.5      0.       0        
GRID        2555       0 0.21684 0.29428      0.       0        
GRID        2556       0 0.29428 0.21684      0.       0        
GRID        2557       0 0.39184 0.16713      0.       0        
GRID        2560       0     0.4      0.      0.       0        
GRID        2561       0     0.3      0.      0.       0        
GRID        2562       0     0.2      0.      0.       0        
GRID        2563       0     0.1      0.      0.       0        
GRID        2564       0 0.13338 0.09174      0.       0        
GRID        2565       0 0.17722 0.18613      0.       0        
GRID        2566       0 0.243070.099855      0.       0        
GRID        2567       0 0.35483  0.1056      0.       0        
GRID        2568       00.096532 0.25017      0.       0        
GRID        2569       0 0.10498 0.35597      0.       0        
GRID        2570       00.081066  0.1653      0.       0        
GRID        2572       0 0.26667     0.5      0.       0        
GRID        2573       0 0.38333     0.5      0.       0        
GRID        2574       0     0.5     0.5      0.       0        
GRID        2581       0 0.16713 0.39184      0.       0        
GRID        2582       0 0.32229 0.31276      0.       0        
GRID        2583       0 0.35582 0.40316      0.       0        
GRID        2584       0 0.41644 0.32133      0.       0        
GRID        2585       0 0.40411 0.23927      0.       0        
GRID        2586       0 0.24645 0.40053      0.       0        
GRID        2587       0      1.      1.      0.       0        
GRID        2588       0      1.     0.9      0.       0        
GRID        2589       0      1.     0.8      0.       0        
GRID        2590       0      1.     0.7      0.       0        
GRID        2591       0      1.     0.6      0.       0        
GRID        2595       0 0.78316 0.70572      0.       0        
GRID        2598       0     0.5    0.85      0.       0        
GRID        2600       0     0.6      1.      0.       0        
GRID        2601       0     0.7      1.      0.       0        
GRID        2602       0     0.8      1.      0.       0        
GRID        2603       0     0.9      1.      0.       0        
GRID        2604       0 0.86662 0.90826      0.       0        
GRID        2605       0 0.82278 0.81387      0.       0        
GRID        2606       0 0.75693 0.90015      0.       0        
GRID        2607       0 0.64517  0.8944      0.       0        
GRID        2608       0 0.90347 0.74983      0.       0        
GRID        2609       0 0.89502 0.64403      0.       0        
GRID        2610       0 0.91893  0.8347      0.       0        
GRID        2615       0     0.5 0.61667      0.       0        
GRID        2618       0 0.60816 0.83287      0.       0        
GRID        2619       0 0.70572 0.78316      0.       0        
GRID        2621       0 0.83287 0.60816      0.       0        
GRID        2622       0 0.67771 0.68724      0.       0        
GRID        2623       0 0.64418 0.59684      0.       0        
GRID        2624       0 0.58356 0.67867      0.       0        
GRID        2625       0 0.59589 0.76073      0.       0        
GRID        2626       0 0.75355 0.59947      0.       0        
GRID        2627       0      1.      0.      0.       0        
GRID        2628       0      1.     0.1      0.       0        
GRID        2629       0      1.     0.2      0.       0        
GRID        2630       0      1.     0.3      0.       0        
GRID        2631       0      1.     0.4      0.       0        
GRID        2632       0      1.     0.5      0.       0        
GRID        2633       0    0.85     0.5      0.       0        
GRID        2635       0 0.78316 0.29428      0.       0        
GRID        2636       0 0.70572 0.21684      0.       0        
GRID        2637       0 0.60816 0.16713      0.       0        
GRID        2638       0     0.5    0.15      0.       0        
GRID        2639       0     0.5      0.      0.       0        
GRID        2640       0     0.6      0.      0.       0        
GRID        2641       0     0.7      0.      0.       0        
GRID        2642       0     0.8      0.      0.       0        
GRID        2643       0     0.9      0.      0.       0        
GRID        2644       0 0.86662 0.09174      0.       0        
GRID        2645       0 0.82278 0.18613      0.       0        
GRID        2646       0 0.756930.099855      0.       0        
GRID        2647       0 0.64517  0.1056      0.       0        
GRID        2648       0 0.90347 0.25017      0.       0        
GRID        2649       0 0.89502 0.35597      0.       0        
GRID        2650       0 0.91893  0.1653      0.       0        
GRID        2652       0 0.73333     0.5      0.       0        
GRID        2653       0 0.61667     0.5      0.       0        
GRID        2655       0     0.5 0.38333      0.       0        
GRID        2656       0     0.5 0.26667      0.       0        
GRID        2661       0 0.83287 0.39184      0.       0        
GRID        2662       0 0.67771 0.31276      0.       0        
GRID        2663       0 0.64418 0.40316      0.       0        
GRID        2664       0 0.58356 0.32133      0.       0        
GRID        2665       0 0.59589 0.23927      0.       0        
GRID        2666       0 0.75355 0.40053      0.       0        
CTRIA3      3545       1    2522    2523    2524                        
CTRIA3      3546       1    2524    2525    2526                        
CTRIA3      3547       1    2522    2524    2526                        
CTRIA3      3548       1    2521    2522    2526                        
CTRIA3      3549       1    2520    2521    2527                        
CTRIA3      3550       1    2598    2519    2520                        
CTRIA3      3551       1    2598    2520    2527                        
CTRIA3      3552       1    2538    2598    2527                        
CTRIA3      3553       1    2527    2521    2526                        
CTRIA3      3554       1    2516    2538    2527                        
CTRIA3      3555       1    2516    2527    2526                        
CTRIA3      3556       1    2516    2526    2525                        
CTRIA3      3557       1    2540    2516    2525                        
CTRIA3      3558       1    2511    2552    2553                        
CTRIA3      3559       1    2553    2541    2529                        
CTRIA3      3560       1    2511    2553    2529                        
CTRIA3      3561       1    2510    2511    2529                        
CTRIA3      3562       1    2540    2525    2528                        
CTRIA3      3563       1    2529    2541    2540                        
CTRIA3      3564       1    2529    2540    2528                        
CTRIA3      3565       1    2510    2529    2528                        
CTRIA3      3566       1    2525    2524    2530                        
CTRIA3      3567       1    2528    2525    2530                        
CTRIA3      3568       1    2509    2510    2528                        
CTRIA3      3569       1    2509    2528    2530                        
CTRIA3      3570       1    2508    2509    2530                        
CTRIA3      3571       1    2508    2530    2524                        
CTRIA3      3572       1    2508    2524    2523                        
CTRIA3      3573       1    2507    2508    2523                        
CTRIA3      3574       2    2542    2544    2545                        
CTRIA3      3575       2    2516    2542    2545                        
CTRIA3      3576       2    2538    2516    2545                        
CTRIA3      3577       2    2598    2538    2545                        
CTRIA3      3578       2    2536    2598    2545                        
CTRIA3      3579       2    2536    2545    2544                        
CTRIA3      3580       2    2615    2536    2544                        
CTRIA3      3581       2    2544    2542    2543                        
CTRIA3      3582       2    2573    2574    2615                        
CTRIA3      3583       2    2615    2544    2543                        
CTRIA3      3584       2    2573    2615    2543                        
CTRIA3      3585       2    2542    2516    2540                        
CTRIA3      3586       2    2540    2541    2546                        
CTRIA3      3587       2    2542    2540    2546                        
CTRIA3      3588       2    2543    2542    2546                        
CTRIA3      3589       2    2572    2573    2543                        
CTRIA3      3590       2    2572    2543    2546                        
CTRIA3      3591       2    2553    2572    2546                        
CTRIA3      3592       2    2553    2546    2541                        
CTRIA3      3593       1    2563    2562    2564                        
CTRIA3      3594       1    2565    2564    2566                        
CTRIA3      3595       1    2564    2562    2566                        
CTRIA3      3596       1    2562    2561    2566                        
CTRIA3      3597       1    2561    2560    2567                        
CTRIA3      3598       1    2639    2638    2560                        
CTRIA3      3599       1    2560    2638    2567                        
CTRIA3      3600       1    2638    2557    2567                        
CTRIA3      3601       1    2561    2567    2566                        
CTRIA3      3602       1    2557    2556    2567                        
CTRIA3      3603       1    2567    2556    2566                        
CTRIA3      3604       1    2566    2556    2565                        
CTRIA3      3605       1    2556    2555    2565                        
CTRIA3      3606       1    2552    2551    2553                        
CTRIA3      3607       1    2581    2553    2569                        
CTRIA3      3608       1    2553    2551    2569                        
CTRIA3      3609       1    2551    2550    2569                        
CTRIA3      3610       1    2565    2555    2568                        
CTRIA3      3611       1    2581    2569    2555                        
CTRIA3      3612       1    2555    2569    2568                        
CTRIA3      3613       1    2569    2550    2568                        
CTRIA3      3614       1    2564    2565    2570                        
CTRIA3      3615       1    2565    2568    2570                        
CTRIA3      3616       1    2550    2549    2568                        
CTRIA3      3617       1    2568    2549    2570                        
CTRIA3      3618       1    2549    2548    2570                        
CTRIA3      3619       1    2570    2548    2564                        
CTRIA3      3620       1    2564    2548    2563                        
CTRIA3      3621       1    2548    2547    2563                        
CTRIA3      3622       2    2584    2582    2585                        
CTRIA3      3623       2    2582    2556    2585                        
CTRIA3      3624       2    2556    2557    2585                        
CTRIA3      3625       2    2557    2638    2585                        
CTRIA3      3626       2    2638    2656    2585                        
CTRIA3      3627       2    2585    2656    2584                        
CTRIA3      3628       2    2656    2655    2584                        
CTRIA3      3629       2    2582    2584    2583                        
CTRIA3      3630       2    2574    2573    2655                        
CTRIA3      3631       2    2584    2655    2583                        
CTRIA3      3632       2    2655    2573    2583                        
CTRIA3      3633       2    2556    2582    2555                        
CTRIA3      3634       2    2581    2555    2586                        
CTRIA3      3635       2    2555    2582    2586                        
CTRIA3      3636       2    2582    2583    2586                        
CTRIA3      3637       2    2573    2572    2583                        
CTRIA3      3638       2    2583    2572    2586                        
CTRIA3      3639       2    2572    2553    2586                        
CTRIA3      3640       2    2586    2553    2581                        
CTRIA3      3641       1    2603    2602    2604                        
CTRIA3      3642       1    2605    2604    2606                        
CTRIA3      3643       1    2604    2602    2606                        
CTRIA3      3644       1    2602    2601    2606                        
CTRIA3      3645       1    2601    2600    2607                        
CTRIA3      3646       1    2519    2598    2600                        
CTRIA3      3647       1    2600    2598    2607                        
CTRIA3      3648       1    2598    2618    2607                        
CTRIA3      3649       1    2601    2607    2606                        
CTRIA3      3650       1    2618    2619    2607                        
CTRIA3      3651       1    2607    2619    2606                        
CTRIA3      3652       1    2606    2619    2605                        
CTRIA3      3653       1    2619    2595    2605                        
CTRIA3      3654       1    2632    2591    2633                        
CTRIA3      3655       1    2621    2633    2609                        
CTRIA3      3656       1    2633    2591    2609                        
CTRIA3      3657       1    2591    2590    2609                        
CTRIA3      3658       1    2605    2595    2608                        
CTRIA3      3659       1    2621    2609    2595                        
CTRIA3      3660       1    2595    2609    2608                        
CTRIA3      3661       1    2609    2590    2608                        
CTRIA3      3662       1    2604    2605    2610                        
CTRIA3      3663       1    2605    2608    2610                        
CTRIA3      3664       1    2590    2589    2608                        
CTRIA3      3665       1    2608    2589    2610                        
CTRIA3      3666       1    2589    2588    2610                        
CTRIA3      3667       1    2610    2588    2604                        
CTRIA3      3668       1    2604    2588    2603                        
CTRIA3      3669       1    2588    2587    2603                        
CTRIA3      3670       2    2624    2622    2625                        
CTRIA3      3671       2    2622    2619    2625                        
CTRIA3      3672       2    2619    2618    2625                        
CTRIA3      3673       2    2618    2598    2625                        
CTRIA3      3674       2    2598    2536    2625                        
CTRIA3      3675       2    2625    2536    2624                        
CTRIA3      3676       2    2536    2615    2624                        
CTRIA3      3677       2    2622    2624    2623                        
CTRIA3      3678       2    2574    2653    2615                        
CTRIA3      3679       2    2624    2615    2623                        
CTRIA3      3680       2    2615    2653    2623                        
CTRIA3      3681       2    2619    2622    2595                        
CTRIA3      3682       2    2621    2595    2626                        
CTRIA3      3683       2    2595    2622    2626                        
CTRIA3      3684       2    2622    2623    2626                        
CTRIA3      3685       2    2653    2652    2623                        
CTRIA3      3686       2    2623    2652    2626                        
CTRIA3      3687       2    2652    2633    2626                        
CTRIA3      3688       2    2626    2633    2621                        
CTRIA3      3689       1    2642    2643    2644                        
CTRIA3      3690       1    2644    2645    2646                        
CTRIA3      3691       1    2642    2644    2646                        
CTRIA3      3692       1    2641    2642    2646                        
CTRIA3      3693       1    2640    2641    2647                        
CTRIA3      3694       1    2638    2639    2640                        
CTRIA3      3695       1    2638    2640    2647                        
CTRIA3      3696       1    2637    2638    2647                        
CTRIA3      3697       1    2647    2641    2646                        
CTRIA3      3698       1    2636    2637    2647                        
CTRIA3      3699       1    2636    2647    2646                        
CTRIA3      3700       1    2636    2646    2645                        
CTRIA3      3701       1    2635    2636    2645                        
CTRIA3      3702       1    2631    2632    2633                        
CTRIA3      3703       1    2633    2661    2649                        
CTRIA3      3704       1    2631    2633    2649                        
CTRIA3      3705       1    2630    2631    2649                        
CTRIA3      3706       1    2635    2645    2648                        
CTRIA3      3707       1    2649    2661    2635                        
CTRIA3      3708       1    2649    2635    2648                        
CTRIA3      3709       1    2630    2649    2648                        
CTRIA3      3710       1    2645    2644    2650                        
CTRIA3      3711       1    2648    2645    2650                        
CTRIA3      3712       1    2629    2630    2648                        
CTRIA3      3713       1    2629    2648    2650                        
CTRIA3      3714       1    2628    2629    2650                        
CTRIA3      3715       1    2628    2650    2644                        
CTRIA3      3716       1    2628    2644    2643                        
CTRIA3      3717       1    2627    2628    2643                        
CTRIA3      3718       2    2662    2664    2665                        
CTRIA3      3719       2    2636    2662    2665                        
CTRIA3      3720       2    2637    2636    2665                        
CTRIA3      3721       2    2638    2637    2665                        
CTRIA3      3722       2    2656    2638    2665                        
CTRIA3      3723       2    2656    2665    2664                        
CTRIA3      3724       2    2655    2656    2664                        
CTRIA3      3725       2    2664    2662    2663                        
CTRIA3      3726       2    2653    2574    2655                        
CTRIA3      3727       2    2655    2664    2663                        
CTRIA3      3728       2    2653    2655    2663                        
CTRIA3      3729       2    2662    2636    2635                        
CTRIA3      3730       2    2635    2661    2666                        
CTRIA3      3731       2    2662    2635    2666                        
CTRIA3      3732       2    2663    2662    2666                        
CTRIA3      3733       2    2652    2653    2663                        
CTRIA3      3734       2    2652    2663    2666                        
CTRIA3      3735       2    2633    2652    2666                        
CTRIA3      3736       2    2633    2666    2661                        
ENDDATA
