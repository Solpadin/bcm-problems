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
$   Date       : Wed Jun 13 20:06:20 2007
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
$ FEMAP Load Set 111 : Inclusion BSOURCE 1.0473 1.0473 0.894 0.894 0.
$ FEMAP Property 1 : Matrix
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Fulleren
PSHELL         2       1      0.       1               1              0.
$ FEMAP Material 1 : mat
MAT1           1                              0.      0.      0.        
GRID           1       0  2.0946      0.      0.       0        
GRID           2       0  2.0946  0.3491      0.       0        
GRID           3       0  2.0946  0.6982      0.       0        
GRID           4       0  2.0946  1.0473      0.       0        
GRID           8       0 1.60217 0.28359      0.       0        
GRID          10       0  1.0473  0.1033      0.       0        
GRID          12       0  1.3964      0.      0.       0        
GRID          13       0  1.7455      0.      0.       0        
GRID          18       0 1.57278 0.32404      0.       0        
GRID          22       0 1.33901  0.1495      0.       0        
GRID          24       0 1.81101 0.49243      0.       0        
GRID          25       0  1.9451 0.75559      0.       0        
GRID          27       0  1.6433  1.0473      0.       0        
GRID          33       0 1.32356 0.19706      0.       0        
GRID          35       0 1.77056 0.52182      0.       0        
GRID          36       0 1.89754 0.77104      0.       0        
GRID          37       0 1.37497 0.48352      0.       0        
GRID          38       0 1.35846 0.76395      0.       0        
GRID          39       0 1.66747 0.77602      0.       0        
GRID          40       0      0.      0.      0.       0        
GRID          41       0      0.  0.3491      0.       0        
GRID          42       0      0.  0.6982      0.       0        
GRID          43       0      0.  1.0473      0.       0        
GRID          44       0  0.1033  1.0473      0.       0        
GRID          45       0  0.1495 0.75559      0.       0        
GRID          46       0 0.28359 0.49243      0.       0        
GRID          47       0 0.49243 0.28359      0.       0        
GRID          48       0 0.75559  0.1495      0.       0        
GRID          50       0  1.0473      0.      0.       0        
GRID          51       0  0.6982      0.      0.       0        
GRID          52       0  0.3491      0.      0.       0        
GRID          58       0 0.77104 0.19706      0.       0        
GRID          66       0  0.4513  1.0473      0.       0        
GRID          67       0  0.7493  1.0473      0.       0        
GRID          68       0  1.0473  1.0473      0.       0        
GRID          69       0  1.0473  0.7493      0.       0        
GRID          70       0  1.0473  0.4513      0.       0        
GRID          71       0  1.0473  0.1533      0.       0        
GRID          73       0 0.52182 0.32404      0.       0        
GRID          74       0 0.32404 0.52182      0.       0        
GRID          75       0 0.19706 0.77104      0.       0        
GRID          76       0 0.71963 0.48352      0.       0        
GRID          77       0 0.73614 0.76395      0.       0        
GRID          78       0 0.42713 0.77602      0.       0        
GRID          79       0  2.0946  2.0946      0.       0        
GRID          80       0  2.0946  1.7455      0.       0        
GRID          81       0  2.0946  1.3964      0.       0        
GRID          83       0  1.9913  1.0473      0.       0        
GRID          84       0  1.9451 1.33901      0.       0        
GRID          87       0 1.33901  1.9451      0.       0        
GRID          89       0  1.0473  2.0946      0.       0        
GRID          90       0  1.3964  2.0946      0.       0        
GRID          91       0  1.7455  2.0946      0.       0        
GRID          93       0  1.9413  1.0473      0.       0        
GRID          94       0 1.89754 1.32356      0.       0        
GRID          95       0 1.77056 1.57278      0.       0        
GRID          96       0 1.57278 1.77056      0.       0        
GRID          97       0 1.32356 1.89754      0.       0        
GRID         101       0 1.60217 1.81101      0.       0        
GRID         102       0 1.81101 1.60217      0.       0        
GRID         106       0  1.3453  1.0473      0.       0        
GRID         108       0  1.0473  1.3453      0.       0        
GRID         109       0  1.0473  1.6433      0.       0        
GRID         115       0 1.37497 1.61108      0.       0        
GRID         116       0 1.35846 1.33065      0.       0        
GRID         117       0 1.66747 1.31858      0.       0        
GRID         118       0      0.  2.0946      0.       0        
GRID         119       0      0.  1.7455      0.       0        
GRID         120       0      0.  1.3964      0.       0        
GRID         126       0 0.75559  1.9451      0.       0        
GRID         129       0  0.6982  2.0946      0.       0        
GRID         130       0  0.3491  2.0946      0.       0        
GRID         132       0  0.1533  1.0473      0.       0        
GRID         133       0 0.19706 1.32356      0.       0        
GRID         134       0 0.32404 1.57278      0.       0        
GRID         135       0 0.52182 1.77056      0.       0        
GRID         136       0 0.77104 1.89754      0.       0        
GRID         137       0  1.0473  1.9413      0.       0        
GRID         138       0  1.0473  1.9913      0.       0        
GRID         140       0 0.49243 1.81101      0.       0        
GRID         141       0 0.28359 1.60217      0.       0        
GRID         142       0  0.1495 1.33901      0.       0        
GRID         154       0 0.71963 1.61108      0.       0        
GRID         155       0 0.73614 1.33065      0.       0        
GRID         156       0 0.42713 1.31858      0.       0        
CQUAD4         1       1       3       4      83      25                
CQUAD4         2       1       2       3      25      24                
CQUAD4         3       1      22      10      50      12                
CQUAD4         4       1       8      22      12      13                
CQUAD4         5       1       2      24       8      13                
CTRIA3         6       1       1       2      13                        
CQUAD4         7       1      83      93      36      25                
CQUAD4         8       1      25      36      35      24                
CQUAD4         9       1      24      35      18       8                
CQUAD4        10       1       8      18      33      22                
CQUAD4        11       1      22      33      71      10                
CQUAD4        12       2      70      71      33      37                
CQUAD4        13       2      69      70      37      38                
CQUAD4        14       2     106      68      69      38                
CTRIA3        15       2      37      33      18                        
CTRIA3        16       2      35      36      39                        
CQUAD4        17       2      38      37      35      39                
CQUAD4        18       2      27     106      38      39                
CQUAD4        19       2      93      27      39      36                
CTRIA3        20       2      37      18      35                        
CQUAD4        21       1      43      42      45      44                
CQUAD4        22       1      42      41      46      45                
CQUAD4        23       1      10      48      51      50                
CQUAD4        24       1      48      47      52      51                
CQUAD4        25       1      46      41      52      47                
CTRIA3        26       1      41      40      52                        
CQUAD4        27       1     132      44      45      75                
CQUAD4        28       1      75      45      46      74                
CQUAD4        29       1      74      46      47      73                
CQUAD4        30       1      73      47      48      58                
CQUAD4        31       1      58      48      10      71                
CQUAD4        32       2      71      70      76      58                
CQUAD4        33       2      70      69      77      76                
CQUAD4        34       2      68      67      77      69                
CTRIA3        35       2      58      76      73                        
CTRIA3        36       2      75      74      78                        
CQUAD4        37       2      76      77      78      74                
CQUAD4        38       2      67      66      78      77                
CQUAD4        39       2      66     132      75      78                
CTRIA3        40       2      73      76      74                        
CQUAD4        41       1       4      81      84      83                
CQUAD4        42       1      81      80     102      84                
CQUAD4        43       1     138      87      90      89                
CQUAD4        44       1      87     101      91      90                
CQUAD4        45       1     102      80      91     101                
CTRIA3        46       1      80      79      91                        
CQUAD4        47       1      93      83      84      94                
CQUAD4        48       1      94      84     102      95                
CQUAD4        49       1      95     102     101      96                
CQUAD4        50       1      96     101      87      97                
CQUAD4        51       1      97      87     138     137                
CQUAD4        52       2     137     109     115      97                
CQUAD4        53       2     109     108     116     115                
CQUAD4        54       2      68     106     116     108                
CTRIA3        55       2      97     115      96                        
CTRIA3        56       2      94      95     117                        
CQUAD4        57       2     115     116     117      95                
CQUAD4        58       2     106      27     117     116                
CQUAD4        59       2      27      93      94     117                
CTRIA3        60       2      96     115      95                        
CQUAD4        61       1     120      43      44     142                
CQUAD4        62       1     119     120     142     141                
CQUAD4        63       1     126     138      89     129                
CQUAD4        64       1     140     126     129     130                
CQUAD4        65       1     119     141     140     130                
CTRIA3        66       1     118     119     130                        
CQUAD4        67       1      44     132     133     142                
CQUAD4        68       1     142     133     134     141                
CQUAD4        69       1     141     134     135     140                
CQUAD4        70       1     140     135     136     126                
CQUAD4        71       1     126     136     137     138                
CQUAD4        72       2     109     137     136     154                
CQUAD4        73       2     108     109     154     155                
CQUAD4        74       2      67      68     108     155                
CTRIA3        75       2     154     136     135                        
CTRIA3        76       2     134     133     156                        
CQUAD4        77       2     155     154     134     156                
CQUAD4        78       2      66      67     155     156                
CQUAD4        79       2     132      66     156     133                
CTRIA3        80       2     154     135     134                        
ENDDATA