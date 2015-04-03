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
$   Date       : Wed Nov 09 07:08:03 2005
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
$ FEMAP Property 1 : Untitled
PSHELL         1       1      0.       1      1.       1 0.83333      0.
$ FEMAP Property 2 : Untitled
PSHELL         2       1      0.       1      1.       1 0.83333      0.
$ FEMAP Material 1 : Untitled
MAT1           1     7.8     2.3     0.3      0.      0.      0.        
GRID        2529       0      0.      1.      0.       0        
GRID        2530       0      0. 0.91667      0.       0        
GRID        2531       0      0. 0.83333      0.       0        
GRID        2532       0      0.    0.75      0.       0        
GRID        2533       0      0. 0.66667      0.       0        
GRID        2534       0      0. 0.58333      0.       0        
GRID        2536       0     0.1     0.5      0.       0        
GRID        2540       0 0.32366 0.74271      0.       0        
GRID        2544       0     0.5      1.      0.       0        
GRID        2545       0 0.41667      1.      0.       0        
GRID        2546       0 0.33333      1.      0.       0        
GRID        2547       0    0.25      1.      0.       0        
GRID        2548       0 0.16667      1.      0.       0        
GRID        2549       00.083333      1.      0.       0        
GRID        2550       0  0.2978 0.92517      0.       0        
GRID        2551       0 0.24523 0.84718      0.       0        
GRID        2552       0 0.21515 0.74942      0.       0        
GRID        2553       0 0.17466 0.67205      0.       0        
GRID        2554       0 0.13763 0.58924      0.       0        
GRID        2555       0 0.34822 0.84542      0.       0        
GRID        2556       0 0.38745 0.92038      0.       0        
GRID        2557       0 0.42859 0.85046      0.       0        
GRID        2558       0 0.15155 0.81831      0.       0        
GRID        2559       0 0.07778 0.87722      0.       0        
GRID        2560       00.082952 0.68411      0.       0        
GRID        2561       00.053517  0.5872      0.       0        
GRID        2562       00.075208  0.7838      0.       0        
GRID        2563       0 0.13972  0.7411      0.       0        
GRID        2564       0 0.22333 0.93703      0.       0        
GRID        2565       0 0.15771  0.9131      0.       0        
GRID        2573       0 0.40729 0.78532      0.       0        
GRID        2575       0 0.25729 0.67634      0.       0        
GRID        2576       0 0.21468 0.59271      0.       0        
GRID        2582       0      0.      0.      0.       0        
GRID        2583       0      0.0.083333      0.       0        
GRID        2584       0      0. 0.16667      0.       0        
GRID        2585       0      0.    0.25      0.       0        
GRID        2586       0      0. 0.33333      0.       0        
GRID        2587       0      0. 0.41667      0.       0        
GRID        2588       0      0.     0.5      0.       0        
GRID        2590       0     0.2     0.5      0.       0        
GRID        2597       0     0.5      0.      0.       0        
GRID        2598       0 0.41667      0.      0.       0        
GRID        2599       0 0.33333      0.      0.       0        
GRID        2600       0    0.25      0.      0.       0        
GRID        2601       0 0.16667      0.      0.       0        
GRID        2602       00.083333      0.      0.       0        
GRID        2603       0  0.29780.074827      0.       0        
GRID        2604       0 0.24523 0.15282      0.       0        
GRID        2605       0 0.21515 0.25058      0.       0        
GRID        2606       0 0.17466 0.32795      0.       0        
GRID        2607       0 0.13763 0.41076      0.       0        
GRID        2608       0 0.34822 0.15458      0.       0        
GRID        2609       0 0.387450.079618      0.       0        
GRID        2610       0 0.42859 0.14954      0.       0        
GRID        2611       0 0.15155 0.18169      0.       0        
GRID        2612       0 0.07778 0.12278      0.       0        
GRID        2613       00.082952 0.31589      0.       0        
GRID        2614       00.053517  0.4128      0.       0        
GRID        2615       00.075208  0.2162      0.       0        
GRID        2616       0 0.13972  0.2589      0.       0        
GRID        2617       0 0.223330.062967      0.       0        
GRID        2618       0 0.157710.086897      0.       0        
GRID        2626       0 0.40729 0.21468      0.       0        
GRID        2627       0 0.32366 0.25729      0.       0        
GRID        2628       0 0.25729 0.32366      0.       0        
GRID        2629       0 0.21468 0.40729      0.       0        
GRID        2635       0      1.      1.      0.       0        
GRID        2636       0      1. 0.91667      0.       0        
GRID        2637       0      1. 0.83333      0.       0        
GRID        2638       0      1.    0.75      0.       0        
GRID        2639       0      1. 0.66667      0.       0        
GRID        2640       0      1. 0.58333      0.       0        
GRID        2642       0     0.9     0.5      0.       0        
GRID        2644       0 0.78532 0.59271      0.       0        
GRID        2646       0 0.67634 0.74271      0.       0        
GRID        2647       0 0.59271 0.78532      0.       0        
GRID        2648       0     0.5     0.8      0.       0        
GRID        2649       0     0.5     0.9      0.       0        
GRID        2651       0 0.58333      1.      0.       0        
GRID        2652       0 0.66667      1.      0.       0        
GRID        2653       0    0.75      1.      0.       0        
GRID        2654       0 0.83333      1.      0.       0        
GRID        2655       0 0.91667      1.      0.       0        
GRID        2656       0  0.7022 0.92517      0.       0        
GRID        2657       0 0.75477 0.84718      0.       0        
GRID        2658       0 0.78485 0.74942      0.       0        
GRID        2659       0 0.82534 0.67205      0.       0        
GRID        2660       0 0.86237 0.58924      0.       0        
GRID        2661       0 0.65178 0.84542      0.       0        
GRID        2662       0 0.61255 0.92038      0.       0        
GRID        2663       0 0.57141 0.85046      0.       0        
GRID        2664       0 0.84845 0.81831      0.       0        
GRID        2665       0 0.92222 0.87722      0.       0        
GRID        2666       0 0.91705 0.68411      0.       0        
GRID        2667       0 0.94648  0.5872      0.       0        
GRID        2668       0 0.92479  0.7838      0.       0        
GRID        2669       0 0.86028  0.7411      0.       0        
GRID        2670       0 0.77667 0.93703      0.       0        
GRID        2671       0 0.84229  0.9131      0.       0        
GRID        2681       0 0.74271 0.67634      0.       0        
GRID        2688       0      1.      0.      0.       0        
GRID        2689       0      1.0.083333      0.       0        
GRID        2690       0      1. 0.16667      0.       0        
GRID        2691       0      1.    0.25      0.       0        
GRID        2692       0      1. 0.33333      0.       0        
GRID        2693       0      1. 0.41667      0.       0        
GRID        2694       0      1.     0.5      0.       0        
GRID        2696       0     0.8     0.5      0.       0        
GRID        2698       0 0.74271 0.32366      0.       0        
GRID        2699       0 0.67634 0.25729      0.       0        
GRID        2700       0 0.59271 0.21468      0.       0        
GRID        2701       0     0.5     0.2      0.       0        
GRID        2702       0     0.5     0.1      0.       0        
GRID        2704       0 0.58333      0.      0.       0        
GRID        2705       0 0.66667      0.      0.       0        
GRID        2706       0    0.75      0.      0.       0        
GRID        2707       0 0.83333      0.      0.       0        
GRID        2708       0 0.91667      0.      0.       0        
GRID        2709       0  0.70220.074827      0.       0        
GRID        2710       0 0.75477 0.15282      0.       0        
GRID        2711       0 0.78485 0.25058      0.       0        
GRID        2712       0 0.82534 0.32795      0.       0        
GRID        2713       0 0.86237 0.41076      0.       0        
GRID        2714       0 0.65178 0.15458      0.       0        
GRID        2715       0 0.612550.079618      0.       0        
GRID        2716       0 0.57141 0.14954      0.       0        
GRID        2717       0 0.84845 0.18169      0.       0        
GRID        2718       0 0.92222 0.12278      0.       0        
GRID        2719       0 0.91705 0.31589      0.       0        
GRID        2720       0 0.94648  0.4128      0.       0        
GRID        2721       0 0.92479  0.2162      0.       0        
GRID        2722       0 0.86028  0.2589      0.       0        
GRID        2723       0 0.776670.062967      0.       0        
GRID        2724       0 0.842290.086897      0.       0        
GRID        2735       0 0.78532 0.40729      0.       0        
CTRIA3      3597       1    2550    2555    2556                        
CTRIA3      3598       1    2546    2550    2556                        
CTRIA3      3599       1    2545    2546    2556                        
CTRIA3      3600       1    2556    2555    2557                        
CTRIA3      3601       1    2649    2544    2545                        
CTRIA3      3602       1    2649    2545    2556                        
CTRIA3      3603       1    2649    2556    2557                        
CTRIA3      3604       1    2648    2649    2557                        
CTRIA3      3605       1    2573    2648    2557                        
CTRIA3      3606       1    2573    2557    2555                        
CTRIA3      3607       1    2555    2550    2551                        
CTRIA3      3608       1    2540    2573    2555                        
CTRIA3      3609       1    2540    2555    2551                        
CTRIA3      3610       1    2540    2551    2552                        
CTRIA3      3611       1    2575    2540    2552                        
CTRIA3      3612       1    2575    2552    2553                        
CTRIA3      3613       1    2576    2575    2553                        
CTRIA3      3614       1    2576    2553    2554                        
CTRIA3      3615       1    2590    2576    2554                        
CTRIA3      3616       1    2536    2590    2554                        
CTRIA3      3617       1    2536    2554    2561                        
CTRIA3      3618       1    2588    2536    2561                        
CTRIA3      3619       1    2534    2588    2561                        
CTRIA3      3620       1    2554    2553    2560                        
CTRIA3      3621       1    2561    2554    2560                        
CTRIA3      3622       1    2533    2534    2561                        
CTRIA3      3623       1    2533    2561    2560                        
CTRIA3      3624       1    2558    2562    2563                        
CTRIA3      3625       1    2552    2558    2563                        
CTRIA3      3626       1    2553    2552    2563                        
CTRIA3      3627       1    2560    2553    2563                        
CTRIA3      3628       1    2532    2533    2560                        
CTRIA3      3629       1    2560    2563    2562                        
CTRIA3      3630       1    2532    2560    2562                        
CTRIA3      3631       1    2531    2532    2562                        
CTRIA3      3632       1    2562    2558    2559                        
CTRIA3      3633       1    2531    2562    2559                        
CTRIA3      3634       1    2530    2531    2559                        
CTRIA3      3635       1    2547    2548    2564                        
CTRIA3      3636       1    2550    2546    2547                        
CTRIA3      3637       1    2550    2547    2564                        
CTRIA3      3638       1    2551    2550    2564                        
CTRIA3      3639       1    2558    2552    2551                        
CTRIA3      3640       1    2551    2564    2565                        
CTRIA3      3641       1    2558    2551    2565                        
CTRIA3      3642       1    2559    2558    2565                        
CTRIA3      3643       1    2565    2564    2548                        
CTRIA3      3644       1    2565    2548    2549                        
CTRIA3      3645       1    2559    2565    2549                        
CTRIA3      3646       1    2530    2559    2549                        
CTRIA3      3647       1    2529    2530    2549                        
CTRIA3      3667       1    2608    2603    2609                        
CTRIA3      3668       1    2603    2599    2609                        
CTRIA3      3669       1    2599    2598    2609                        
CTRIA3      3670       1    2608    2609    2610                        
CTRIA3      3671       1    2597    2702    2598                        
CTRIA3      3672       1    2598    2702    2609                        
CTRIA3      3673       1    2609    2702    2610                        
CTRIA3      3674       1    2702    2701    2610                        
CTRIA3      3675       1    2701    2626    2610                        
CTRIA3      3676       1    2610    2626    2608                        
CTRIA3      3677       1    2603    2608    2604                        
CTRIA3      3678       1    2626    2627    2608                        
CTRIA3      3679       1    2608    2627    2604                        
CTRIA3      3680       1    2604    2627    2605                        
CTRIA3      3681       1    2627    2628    2605                        
CTRIA3      3682       1    2605    2628    2606                        
CTRIA3      3683       1    2628    2629    2606                        
CTRIA3      3684       1    2606    2629    2607                        
CTRIA3      3685       1    2629    2590    2607                        
CTRIA3      3686       1    2590    2536    2607                        
CTRIA3      3687       1    2607    2536    2614                        
CTRIA3      3688       1    2536    2588    2614                        
CTRIA3      3689       1    2588    2587    2614                        
CTRIA3      3690       1    2606    2607    2613                        
CTRIA3      3691       1    2607    2614    2613                        
CTRIA3      3692       1    2587    2586    2614                        
CTRIA3      3693       1    2614    2586    2613                        
CTRIA3      3694       1    2615    2611    2616                        
CTRIA3      3695       1    2611    2605    2616                        
CTRIA3      3696       1    2605    2606    2616                        
CTRIA3      3697       1    2606    2613    2616                        
CTRIA3      3698       1    2586    2585    2613                        
CTRIA3      3699       1    2616    2613    2615                        
CTRIA3      3700       1    2613    2585    2615                        
CTRIA3      3701       1    2585    2584    2615                        
CTRIA3      3702       1    2611    2615    2612                        
CTRIA3      3703       1    2615    2584    2612                        
CTRIA3      3704       1    2584    2583    2612                        
CTRIA3      3705       1    2601    2600    2617                        
CTRIA3      3706       1    2599    2603    2600                        
CTRIA3      3707       1    2600    2603    2617                        
CTRIA3      3708       1    2603    2604    2617                        
CTRIA3      3709       1    2605    2611    2604                        
CTRIA3      3710       1    2617    2604    2618                        
CTRIA3      3711       1    2604    2611    2618                        
CTRIA3      3712       1    2611    2612    2618                        
CTRIA3      3713       1    2617    2618    2601                        
CTRIA3      3714       1    2601    2618    2602                        
CTRIA3      3715       1    2618    2612    2602                        
CTRIA3      3716       1    2612    2583    2602                        
CTRIA3      3717       1    2583    2582    2602                        
CTRIA3      3737       1    2661    2656    2662                        
CTRIA3      3738       1    2656    2652    2662                        
CTRIA3      3739       1    2652    2651    2662                        
CTRIA3      3740       1    2661    2662    2663                        
CTRIA3      3741       1    2544    2649    2651                        
CTRIA3      3742       1    2651    2649    2662                        
CTRIA3      3743       1    2662    2649    2663                        
CTRIA3      3744       1    2649    2648    2663                        
CTRIA3      3745       1    2648    2647    2663                        
CTRIA3      3746       1    2663    2647    2661                        
CTRIA3      3747       1    2656    2661    2657                        
CTRIA3      3748       1    2647    2646    2661                        
CTRIA3      3749       1    2661    2646    2657                        
CTRIA3      3750       1    2657    2646    2658                        
CTRIA3      3751       1    2646    2681    2658                        
CTRIA3      3752       1    2658    2681    2659                        
CTRIA3      3753       1    2681    2644    2659                        
CTRIA3      3754       1    2659    2644    2660                        
CTRIA3      3755       1    2644    2696    2660                        
CTRIA3      3756       1    2696    2642    2660                        
CTRIA3      3757       1    2660    2642    2667                        
CTRIA3      3758       1    2642    2694    2667                        
CTRIA3      3759       1    2694    2640    2667                        
CTRIA3      3760       1    2659    2660    2666                        
CTRIA3      3761       1    2660    2667    2666                        
CTRIA3      3762       1    2640    2639    2667                        
CTRIA3      3763       1    2667    2639    2666                        
CTRIA3      3764       1    2668    2664    2669                        
CTRIA3      3765       1    2664    2658    2669                        
CTRIA3      3766       1    2658    2659    2669                        
CTRIA3      3767       1    2659    2666    2669                        
CTRIA3      3768       1    2639    2638    2666                        
CTRIA3      3769       1    2669    2666    2668                        
CTRIA3      3770       1    2666    2638    2668                        
CTRIA3      3771       1    2638    2637    2668                        
CTRIA3      3772       1    2664    2668    2665                        
CTRIA3      3773       1    2668    2637    2665                        
CTRIA3      3774       1    2637    2636    2665                        
CTRIA3      3775       1    2654    2653    2670                        
CTRIA3      3776       1    2652    2656    2653                        
CTRIA3      3777       1    2653    2656    2670                        
CTRIA3      3778       1    2656    2657    2670                        
CTRIA3      3779       1    2658    2664    2657                        
CTRIA3      3780       1    2670    2657    2671                        
CTRIA3      3781       1    2657    2664    2671                        
CTRIA3      3782       1    2664    2665    2671                        
CTRIA3      3783       1    2670    2671    2654                        
CTRIA3      3784       1    2654    2671    2655                        
CTRIA3      3785       1    2671    2665    2655                        
CTRIA3      3786       1    2665    2636    2655                        
CTRIA3      3787       1    2636    2635    2655                        
CTRIA3      3807       1    2709    2714    2715                        
CTRIA3      3808       1    2705    2709    2715                        
CTRIA3      3809       1    2704    2705    2715                        
CTRIA3      3810       1    2715    2714    2716                        
CTRIA3      3811       1    2702    2597    2704                        
CTRIA3      3812       1    2702    2704    2715                        
CTRIA3      3813       1    2702    2715    2716                        
CTRIA3      3814       1    2701    2702    2716                        
CTRIA3      3815       1    2700    2701    2716                        
CTRIA3      3816       1    2700    2716    2714                        
CTRIA3      3817       1    2714    2709    2710                        
CTRIA3      3818       1    2699    2700    2714                        
CTRIA3      3819       1    2699    2714    2710                        
CTRIA3      3820       1    2699    2710    2711                        
CTRIA3      3821       1    2698    2699    2711                        
CTRIA3      3822       1    2698    2711    2712                        
CTRIA3      3823       1    2735    2698    2712                        
CTRIA3      3824       1    2735    2712    2713                        
CTRIA3      3825       1    2696    2735    2713                        
CTRIA3      3826       1    2642    2696    2713                        
CTRIA3      3827       1    2642    2713    2720                        
CTRIA3      3828       1    2694    2642    2720                        
CTRIA3      3829       1    2693    2694    2720                        
CTRIA3      3830       1    2713    2712    2719                        
CTRIA3      3831       1    2720    2713    2719                        
CTRIA3      3832       1    2692    2693    2720                        
CTRIA3      3833       1    2692    2720    2719                        
CTRIA3      3834       1    2717    2721    2722                        
CTRIA3      3835       1    2711    2717    2722                        
CTRIA3      3836       1    2712    2711    2722                        
CTRIA3      3837       1    2719    2712    2722                        
CTRIA3      3838       1    2691    2692    2719                        
CTRIA3      3839       1    2719    2722    2721                        
CTRIA3      3840       1    2691    2719    2721                        
CTRIA3      3841       1    2690    2691    2721                        
CTRIA3      3842       1    2721    2717    2718                        
CTRIA3      3843       1    2690    2721    2718                        
CTRIA3      3844       1    2689    2690    2718                        
CTRIA3      3845       1    2706    2707    2723                        
CTRIA3      3846       1    2709    2705    2706                        
CTRIA3      3847       1    2709    2706    2723                        
CTRIA3      3848       1    2710    2709    2723                        
CTRIA3      3849       1    2717    2711    2710                        
CTRIA3      3850       1    2710    2723    2724                        
CTRIA3      3851       1    2717    2710    2724                        
CTRIA3      3852       1    2718    2717    2724                        
CTRIA3      3853       1    2724    2723    2707                        
CTRIA3      3854       1    2724    2707    2708                        
CTRIA3      3855       1    2718    2724    2708                        
CTRIA3      3856       1    2689    2718    2708                        
CTRIA3      3857       1    2688    2689    2708                        
ENDDATA