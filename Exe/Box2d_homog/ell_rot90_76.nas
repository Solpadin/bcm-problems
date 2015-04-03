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
$   Date       : Sat Nov 05 17:20:03 2005
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
$ FEMAP Load Set 111 : Inclusion BSOURCE 0.5 0.5 0.35 0.15 90.
$ FEMAP Property 1 : Matrix
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Fulleren
PSHELL         2       1      0.       1               1              0.
$ FEMAP Material 1 : Untitled
MAT1           1     7.8     2.3     0.3      0.      0.      0.        
GRID        2813       0     0.5 0.61667      0.       0        
GRID        2816       0   0.425 0.80311      0.       0        
GRID        2818       0      0.      1.      0.       0        
GRID        2819       0      0.   0.875      0.       0        
GRID        2820       0      0.    0.75      0.       0        
GRID        2821       0      0.   0.625      0.       0        
GRID        2822       0      0.     0.5      0.       0        
GRID        2823       0 0.11667     0.5      0.       0        
GRID        2826       0  0.3701   0.675      0.       0        
GRID        2829       0     0.5      1.      0.       0        
GRID        2830       0   0.375      1.      0.       0        
GRID        2831       0    0.25      1.      0.       0        
GRID        2832       0   0.125      1.      0.       0        
GRID        2833       0 0.11575 0.87063      0.       0        
GRID        2834       0 0.10948 0.74491      0.       0        
GRID        2835       0 0.10944 0.62168      0.       0        
GRID        2836       0  0.2277 0.86343      0.       0        
GRID        2837       0 0.21258 0.73543      0.       0        
GRID        2838       0 0.21074 0.61766      0.       0        
GRID        2839       0 0.33231 0.84578      0.       0        
GRID        2840       0 0.30065 0.71743      0.       0        
GRID        2841       0 0.28744 0.61139      0.       0        
GRID        2845       0     0.5 0.73333      0.       0        
GRID        2846       0     0.5    0.85      0.       0        
GRID        2848       0  0.6299   0.675      0.       0        
GRID        2849       0      1.      1.      0.       0        
GRID        2850       0      1.   0.875      0.       0        
GRID        2851       0      1.    0.75      0.       0        
GRID        2852       0      1.   0.625      0.       0        
GRID        2858       0   0.575 0.80311      0.       0        
GRID        2861       0   0.625      1.      0.       0        
GRID        2862       0    0.75      1.      0.       0        
GRID        2863       0   0.875      1.      0.       0        
GRID        2864       0 0.88425 0.87063      0.       0        
GRID        2865       0 0.89052 0.74491      0.       0        
GRID        2866       0 0.89056 0.62168      0.       0        
GRID        2867       0  0.7723 0.86343      0.       0        
GRID        2868       0 0.78742 0.73543      0.       0        
GRID        2869       0 0.78926 0.61766      0.       0        
GRID        2870       0 0.66769 0.84578      0.       0        
GRID        2871       0 0.69935 0.71743      0.       0        
GRID        2872       0 0.71256 0.61139      0.       0        
GRID        2875       0     0.5 0.38333      0.       0        
GRID        2876       0     0.5 0.26667      0.       0        
GRID        2879       0  0.3701   0.325      0.       0        
GRID        2880       0      0.      0.      0.       0        
GRID        2881       0      0.   0.125      0.       0        
GRID        2882       0      0.    0.25      0.       0        
GRID        2883       0      0.   0.375      0.       0        
GRID        2886       0 0.23333     0.5      0.       0        
GRID        2887       0    0.35     0.5      0.       0        
GRID        2889       0   0.425 0.19689      0.       0        
GRID        2890       0     0.5    0.15      0.       0        
GRID        2891       0     0.5      0.      0.       0        
GRID        2892       0   0.375      0.      0.       0        
GRID        2893       0    0.25      0.      0.       0        
GRID        2894       0   0.125      0.      0.       0        
GRID        2895       0 0.11575 0.12937      0.       0        
GRID        2896       0 0.10948 0.25509      0.       0        
GRID        2897       0 0.10944 0.37832      0.       0        
GRID        2898       0  0.2277 0.13657      0.       0        
GRID        2899       0 0.21258 0.26457      0.       0        
GRID        2900       0 0.21074 0.38234      0.       0        
GRID        2901       0 0.33231 0.15422      0.       0        
GRID        2902       0 0.30065 0.28257      0.       0        
GRID        2903       0 0.28744 0.38861      0.       0        
GRID        2905       0     0.5     0.5      0.       0        
GRID        2909       0   0.575 0.19689      0.       0        
GRID        2911       0      1.      0.      0.       0        
GRID        2912       0      1.   0.125      0.       0        
GRID        2913       0      1.    0.25      0.       0        
GRID        2914       0      1.   0.375      0.       0        
GRID        2915       0      1.     0.5      0.       0        
GRID        2916       0 0.88333     0.5      0.       0        
GRID        2917       0 0.76667     0.5      0.       0        
GRID        2918       0    0.65     0.5      0.       0        
GRID        2919       0  0.6299   0.325      0.       0        
GRID        2923       0   0.625      0.      0.       0        
GRID        2924       0    0.75      0.      0.       0        
GRID        2925       0   0.875      0.      0.       0        
GRID        2926       0 0.88425 0.12937      0.       0        
GRID        2927       0 0.89052 0.25509      0.       0        
GRID        2928       0 0.89056 0.37832      0.       0        
GRID        2929       0  0.7723 0.13657      0.       0        
GRID        2930       0 0.78742 0.26457      0.       0        
GRID        2931       0 0.78926 0.38234      0.       0        
GRID        2932       0 0.66769 0.15422      0.       0        
GRID        2933       0 0.69935 0.28257      0.       0        
GRID        2934       0 0.71256 0.38861      0.       0        
CTRIA3      3865       2    2845    2846    2816                        
CQUAD4      3866       2    2813    2845    2816    2826                
CQUAD4      3867       2    2887    2905    2813    2826                
CQUAD4      3868       1    2830    2831    2836    2839                
CTRIA3      3869       1    2846    2829    2830                        
CQUAD4      3870       1    2816    2846    2830    2839                
CQUAD4      3871       1    2839    2836    2837    2840                
CQUAD4      3872       1    2826    2816    2839    2840                
CQUAD4      3873       1    2840    2837    2838    2841                
CQUAD4      3874       1    2887    2826    2840    2841                
CQUAD4      3875       1    2886    2887    2841    2838                
CQUAD4      3876       1    2836    2831    2832    2833                
CQUAD4      3877       1    2837    2836    2833    2834                
CQUAD4      3878       1    2838    2837    2834    2835                
CQUAD4      3879       1    2823    2886    2838    2835                
CQUAD4      3880       1    2821    2822    2823    2835                
CQUAD4      3881       1    2820    2821    2835    2834                
CQUAD4      3882       1    2819    2820    2834    2833                
CQUAD4      3883       1    2818    2819    2833    2832                
CTRIA3      3884       2    2846    2845    2858                        
CQUAD4      3885       2    2845    2813    2848    2858                
CQUAD4      3886       2    2905    2918    2848    2813                
CQUAD4      3887       1    2862    2861    2870    2867                
CTRIA3      3888       1    2829    2846    2861                        
CQUAD4      3889       1    2846    2858    2870    2861                
CQUAD4      3890       1    2867    2870    2871    2868                
CQUAD4      3891       1    2858    2848    2871    2870                
CQUAD4      3892       1    2868    2871    2872    2869                
CQUAD4      3893       1    2848    2918    2872    2871                
CQUAD4      3894       1    2918    2917    2869    2872                
CQUAD4      3895       1    2862    2867    2864    2863                
CQUAD4      3896       1    2867    2868    2865    2864                
CQUAD4      3897       1    2868    2869    2866    2865                
CQUAD4      3898       1    2917    2916    2866    2869                
CQUAD4      3899       1    2915    2852    2866    2916                
CQUAD4      3900       1    2852    2851    2865    2866                
CQUAD4      3901       1    2851    2850    2864    2865                
CQUAD4      3902       1    2850    2849    2863    2864                
CTRIA3      3903       2    2890    2876    2889                        
CQUAD4      3904       2    2876    2875    2879    2889                
CQUAD4      3905       2    2905    2887    2879    2875                
CQUAD4      3906       1    2893    2892    2901    2898                
CTRIA3      3907       1    2891    2890    2892                        
CQUAD4      3908       1    2890    2889    2901    2892                
CQUAD4      3909       1    2898    2901    2902    2899                
CQUAD4      3910       1    2889    2879    2902    2901                
CQUAD4      3911       1    2899    2902    2903    2900                
CQUAD4      3912       1    2879    2887    2903    2902                
CQUAD4      3913       1    2887    2886    2900    2903                
CQUAD4      3914       1    2893    2898    2895    2894                
CQUAD4      3915       1    2898    2899    2896    2895                
CQUAD4      3916       1    2899    2900    2897    2896                
CQUAD4      3917       1    2886    2823    2897    2900                
CQUAD4      3918       1    2822    2883    2897    2823                
CQUAD4      3919       1    2883    2882    2896    2897                
CQUAD4      3920       1    2882    2881    2895    2896                
CQUAD4      3921       1    2881    2880    2894    2895                
CTRIA3      3922       2    2876    2890    2909                        
CQUAD4      3923       2    2875    2876    2909    2919                
CQUAD4      3924       2    2918    2905    2875    2919                
CQUAD4      3925       1    2923    2924    2929    2932                
CTRIA3      3926       1    2890    2891    2923                        
CQUAD4      3927       1    2909    2890    2923    2932                
CQUAD4      3928       1    2932    2929    2930    2933                
CQUAD4      3929       1    2919    2909    2932    2933                
CQUAD4      3930       1    2933    2930    2931    2934                
CQUAD4      3931       1    2918    2919    2933    2934                
CQUAD4      3932       1    2917    2918    2934    2931                
CQUAD4      3933       1    2929    2924    2925    2926                
CQUAD4      3934       1    2930    2929    2926    2927                
CQUAD4      3935       1    2931    2930    2927    2928                
CQUAD4      3936       1    2916    2917    2931    2928                
CQUAD4      3937       1    2914    2915    2916    2928                
CQUAD4      3938       1    2913    2914    2928    2927                
CQUAD4      3939       1    2912    2913    2927    2926                
CQUAD4      3940       1    2911    2912    2926    2925                
ENDDATA
