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
$   Date       : Thu Jun 14 18:02:57 2007
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
$ FEMAP Load Set 111 : Inclusion BSOURCE 0.3821 0.3821 0.319 0.319 0.
$ FEMAP Property 1 : Matrix
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Fulleren
PSHELL         2       1      0.       1               1              0.
$ FEMAP Material 1 : mat
MAT1           1                              0.      0.      0.        
GRID           3       0 0.25237  0.6953      0.       0        
GRID           5       00.068905 0.51183      0.       0        
GRID           7       0  0.0631  0.3821      0.       0        
GRID           9       0 0.15653 0.60767      0.       0        
GRID          13       0 0.14239 0.62181      0.       0        
GRID          15       0  0.3821  0.7211      0.       0        
GRID          16       0  0.3821  0.7642      0.       0        
GRID          17       0 0.25473  0.7642      0.       0        
GRID          18       0 0.12737  0.7642      0.       0        
GRID          19       0      0.  0.7642      0.       0        
GRID          20       0      0. 0.63683      0.       0        
GRID          21       0      0. 0.50947      0.       0        
GRID          22       0      0.  0.3821      0.       0        
GRID          25       0  0.3821 0.59477      0.       0        
GRID          27       0 0.26002 0.67682      0.       0        
GRID          29       00.087382 0.50418      0.       0        
GRID          33       0 0.18693 0.48119      0.       0        
GRID          34       0 0.28176 0.48226      0.       0        
GRID          35       0 0.28326 0.57752      0.       0        
GRID          36       0 0.20857 0.55539      0.       0        
GRID          39       0 0.252370.068905      0.       0        
GRID          41       00.068905 0.25237      0.       0        
GRID          44       00.087382 0.26002      0.       0        
GRID          45       0 0.15653 0.15653      0.       0        
GRID          46       0 0.260020.087382      0.       0        
GRID          47       0  0.0431  0.3821      0.       0        
GRID          49       0 0.14239 0.14239      0.       0        
GRID          52       0  0.3821      0.      0.       0        
GRID          53       0 0.25473      0.      0.       0        
GRID          54       0 0.12737      0.      0.       0        
GRID          55       0      0.      0.      0.       0        
GRID          56       0      0. 0.12737      0.       0        
GRID          57       0      0. 0.25473      0.       0        
GRID          60       0  0.3821 0.27577      0.       0        
GRID          67       0 0.16943  0.3821      0.       0        
GRID          68       0 0.27577  0.3821      0.       0        
GRID          69       0 0.18693 0.28301      0.       0        
GRID          70       0 0.28176 0.28194      0.       0        
GRID          71       0 0.28326 0.18668      0.       0        
GRID          72       0 0.20857 0.20881      0.       0        
GRID          75       0 0.51183  0.6953      0.       0        
GRID          76       0 0.62181 0.62181      0.       0        
GRID          78       0  0.7211  0.3821      0.       0        
GRID          79       0  0.7011  0.3821      0.       0        
GRID          80       0 0.67682 0.50418      0.       0        
GRID          82       0 0.50418 0.67682      0.       0        
GRID          84       0  0.6953 0.51183      0.       0        
GRID          89       0 0.50947  0.7642      0.       0        
GRID          90       0 0.63683  0.7642      0.       0        
GRID          91       0  0.7642  0.7642      0.       0        
GRID          92       0  0.7642 0.63683      0.       0        
GRID          93       0  0.7642 0.50947      0.       0        
GRID          94       0  0.7642  0.3821      0.       0        
GRID          96       0  0.3821 0.48843      0.       0        
GRID          98       0  0.3821  0.7011      0.       0        
GRID         100       0 0.60767 0.60767      0.       0        
GRID         104       0 0.48843  0.3821      0.       0        
GRID         105       0 0.57727 0.48119      0.       0        
GRID         106       0 0.48244 0.48226      0.       0        
GRID         107       0 0.48094 0.57752      0.       0        
GRID         108       0 0.55563 0.55539      0.       0        
GRID         110       0  0.3821  0.0431      0.       0        
GRID         113       0  0.6953 0.25237      0.       0        
GRID         118       0 0.504180.087382      0.       0        
GRID         121       0 0.62181 0.14239      0.       0        
GRID         122       0 0.511830.068905      0.       0        
GRID         125       0 0.50947      0.      0.       0        
GRID         126       0 0.63683      0.      0.       0        
GRID         127       0  0.7642      0.      0.       0        
GRID         128       0  0.7642 0.12737      0.       0        
GRID         129       0  0.7642 0.25473      0.       0        
GRID         131       0  0.3821  0.3821      0.       0        
GRID         133       0  0.3821 0.16943      0.       0        
GRID         134       0  0.3821  0.0631      0.       0        
GRID         136       0 0.60767 0.15653      0.       0        
GRID         137       0 0.67682 0.26002      0.       0        
GRID         139       0 0.59477  0.3821      0.       0        
GRID         141       0 0.57727 0.28301      0.       0        
GRID         142       0 0.48244 0.28194      0.       0        
GRID         143       0 0.48094 0.18668      0.       0        
GRID         144       0 0.55563 0.20881      0.       0        
CQUAD4         1       1      98      15       3      27                
CQUAD4         2       1      27       3      13       9                
CQUAD4         3       1       9      13       5      29                
CQUAD4         4       1      29       5      47       7                
CQUAD4         5       1       3      15      16      17                
CQUAD4         6       1      13       3      17      18                
CQUAD4         7       1      13      18      19      20                
CQUAD4         8       1       5      13      20      21                
CQUAD4         9       1      47       5      21      22                
CQUAD4        10       2       9      29      33      36                
CQUAD4        11       2      36      33      34      35                
CQUAD4        12       2      27       9      36      35                
CQUAD4        13       2      25      98      27      35                
CQUAD4        14       2      96      25      35      34                
CQUAD4        15       2      33      29       7      67                
CQUAD4        16       2      34      33      67      68                
CQUAD4        17       2     131      96      34      68                
CQUAD4        18       1     110     134      46      39                
CQUAD4        19       1      39      46      45      49                
CQUAD4        20       1      49      45      44      41                
CQUAD4        21       1      41      44       7      47                
CQUAD4        22       1     110      39      53      52                
CQUAD4        23       1      39      49      54      53                
CQUAD4        24       1      54      49      56      55                
CQUAD4        25       1      49      41      57      56                
CQUAD4        26       1      41      47      22      57                
CQUAD4        27       2      44      45      72      69                
CQUAD4        28       2      69      72      71      70                
CQUAD4        29       2      45      46      71      72                
CQUAD4        30       2     134     133      71      46                
CQUAD4        31       2     133      60      70      71                
CQUAD4        32       2      44      69      67       7                
CQUAD4        33       2      69      70      68      67                
CQUAD4        34       2      60     131      68      70                
CQUAD4        35       1      15      98      82      75                
CQUAD4        36       1      75      82     100      76                
CQUAD4        37       1      76     100      80      84                
CQUAD4        38       1      84      80      79      78                
CQUAD4        39       1      15      75      89      16                
CQUAD4        40       1      75      76      90      89                
CQUAD4        41       1      90      76      92      91                
CQUAD4        42       1      76      84      93      92                
CQUAD4        43       1      84      78      94      93                
CQUAD4        44       2      80     100     108     105                
CQUAD4        45       2     105     108     107     106                
CQUAD4        46       2     100      82     107     108                
CQUAD4        47       2      98      25     107      82                
CQUAD4        48       2      25      96     106     107                
CQUAD4        49       2      80     105     139      79                
CQUAD4        50       2     105     106     104     139                
CQUAD4        51       2      96     131     104     106                
CQUAD4        52       1     134     110     122     118                
CQUAD4        53       1     118     122     121     136                
CQUAD4        54       1     136     121     113     137                
CQUAD4        55       1     137     113      78      79                
CQUAD4        56       1     122     110      52     125                
CQUAD4        57       1     121     122     125     126                
CQUAD4        58       1     121     126     127     128                
CQUAD4        59       1     113     121     128     129                
CQUAD4        60       1      78     113     129      94                
CQUAD4        61       2     136     137     141     144                
CQUAD4        62       2     144     141     142     143                
CQUAD4        63       2     118     136     144     143                
CQUAD4        64       2     133     134     118     143                
CQUAD4        65       2      60     133     143     142                
CQUAD4        66       2     141     137      79     139                
CQUAD4        67       2     142     141     139     104                
CQUAD4        68       2     131      60     142     104                
ENDDATA
