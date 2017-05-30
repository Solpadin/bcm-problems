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
$   Date       : Wed Jun 13 11:12:17 2007
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
$ FEMAP Load Set 111 : Inclusion BSOURCE 1. 0.4 0.8 0.2 0.
$ FEMAP Property 1 : Matrix
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Fulleren
PSHELL         2       1      0.       1               1              0.
$ FEMAP Material 1 : mat
MAT1           1                              0.      0.      0.        
GRID           3       0 0.74352 0.61874      0.       0        
GRID           5       0 0.32852 0.53519      0.       0        
GRID           6       0 0.21062 0.47107      0.       0        
GRID          10       0 0.35279 0.51756      0.       0        
GRID          13       0      1.     0.4      0.       0        
GRID          14       0      1.     0.6      0.       0        
GRID          15       0 0.75279 0.59021      0.       0        
GRID          16       0 0.52977  0.5618      0.       0        
GRID          18       0 0.23915  0.4618      0.       0        
GRID          19       0     0.2     0.4      0.       0        
GRID          21       0     0.6     0.4      0.       0        
GRID          22       0     0.8     0.4      0.       0        
GRID          26       0 0.51214 0.58607      0.       0        
GRID          30       0 0.83333     0.8      0.       0        
GRID          31       0 0.66667     0.8      0.       0        
GRID          32       0     0.5     0.8      0.       0        
GRID          33       0 0.33333     0.8      0.       0        
GRID          34       0 0.16667     0.8      0.       0        
GRID          35       0      0.     0.8      0.       0        
GRID          36       0      0.     0.6      0.       0        
GRID          37       0      0.     0.4      0.       0        
GRID          38       0 0.17831 0.69214      0.       0        
GRID          39       0 0.19043 0.58368      0.       0        
GRID          40       0 0.25899 0.60589      0.       0        
GRID          41       0 0.25663 0.69889      0.       0        
GRID          42       0 0.11481 0.56743      0.       0        
GRID          43       00.097494 0.68614      0.       0        
GRID          45       0      1.    0.63      0.       0        
GRID          46       0 1.25648 0.61874      0.       0        
GRID          47       0 1.48786 0.58607      0.       0        
GRID          49       0 1.78938 0.47107      0.       0        
GRID          55       0 1.24721 0.59021      0.       0        
GRID          59       0 1.47023  0.5618      0.       0        
GRID          60       0 1.64721 0.51756      0.       0        
GRID          61       0 1.76085  0.4618      0.       0        
GRID          63       0     1.6     0.4      0.       0        
GRID          64       0     1.4     0.4      0.       0        
GRID          65       0     1.2     0.4      0.       0        
GRID          68       0 1.67148 0.53519      0.       0        
GRID          72       0      1.     0.8      0.       0        
GRID          73       0 1.16667     0.8      0.       0        
GRID          74       0 1.33333     0.8      0.       0        
GRID          75       0     1.5     0.8      0.       0        
GRID          76       0 1.66667     0.8      0.       0        
GRID          77       0 1.83333     0.8      0.       0        
GRID          78       0      2.     0.8      0.       0        
GRID          79       0      2.     0.6      0.       0        
GRID          81       0 1.82169 0.69214      0.       0        
GRID          82       0 1.80957 0.58368      0.       0        
GRID          83       0 1.74101 0.60589      0.       0        
GRID          84       0 1.74337 0.69889      0.       0        
GRID          85       0 1.88519 0.56743      0.       0        
GRID          86       0 1.90251 0.68614      0.       0        
GRID          87       0      1.     0.2      0.       0        
GRID          88       0      1.    0.17      0.       0        
GRID          89       0 0.74352 0.18126      0.       0        
GRID          90       0 0.51214 0.21393      0.       0        
GRID          92       0 0.21062 0.32893      0.       0        
GRID          95       0 0.23915  0.3382      0.       0        
GRID          96       0 0.35279 0.28244      0.       0        
GRID          98       0 0.75279 0.20979      0.       0        
GRID         102       0 0.52977  0.2382      0.       0        
GRID         106       0     0.4     0.4      0.       0        
GRID         109       0    0.17     0.4      0.       0        
GRID         111       0 0.32852 0.26481      0.       0        
GRID         115       0      1.      0.      0.       0        
GRID         116       0 0.83333      0.      0.       0        
GRID         117       0 0.66667      0.      0.       0        
GRID         118       0     0.5      0.      0.       0        
GRID         119       0 0.33333      0.      0.       0        
GRID         120       0 0.16667      0.      0.       0        
GRID         121       0      0.      0.      0.       0        
GRID         122       0      0.     0.2      0.       0        
GRID         124       0 0.17831 0.10786      0.       0        
GRID         125       0 0.19043 0.21632      0.       0        
GRID         126       0 0.25899 0.19411      0.       0        
GRID         127       0 0.25663 0.10111      0.       0        
GRID         128       0 0.11481 0.23257      0.       0        
GRID         129       00.097494 0.11386      0.       0        
GRID         135       0 1.78938 0.32893      0.       0        
GRID         137       0     1.8     0.4      0.       0        
GRID         138       0 1.76085  0.3382      0.       0        
GRID         139       0 1.64721 0.28244      0.       0        
GRID         140       0 1.47023  0.2382      0.       0        
GRID         141       0 1.24721 0.20979      0.       0        
GRID         152       0    1.83     0.4      0.       0        
GRID         154       0 1.67148 0.26481      0.       0        
GRID         155       0 1.48786 0.21393      0.       0        
GRID         156       0 1.25648 0.18126      0.       0        
GRID         159       0 1.16667      0.      0.       0        
GRID         160       0 1.33333      0.      0.       0        
GRID         161       0     1.5      0.      0.       0        
GRID         162       0 1.66667      0.      0.       0        
GRID         163       0 1.83333      0.      0.       0        
GRID         164       0      2.      0.      0.       0        
GRID         165       0      2.     0.2      0.       0        
GRID         166       0      2.     0.4      0.       0        
GRID         167       0 1.82169 0.10786      0.       0        
GRID         168       0 1.80957 0.21632      0.       0        
GRID         169       0 1.74101 0.19411      0.       0        
GRID         170       0 1.74337 0.10111      0.       0        
GRID         171       0 1.88519 0.23257      0.       0        
GRID         172       0 1.90251 0.11386      0.       0        
CQUAD4         1       1      14      45       3      15                
CQUAD4         2       1      15       3      26      16                
CQUAD4         3       1      16      26       5      10                
CQUAD4         4       1      10       5       6      18                
CQUAD4         5       1      18       6     109      19                
CQUAD4         6       2      10      18      19     106                
CQUAD4         7       2      16      10     106      21                
CQUAD4         8       2      15      16      21      22                
CQUAD4         9       2      13      14      15      22                
CQUAD4        10       1       3      45      72      30                
CTRIA3        11       1       3      30      31                        
CQUAD4        12       1      26       3      31      32                
CQUAD4        13       1       5      26      32      33                
CQUAD4        14       1      33      34      38      41                
CQUAD4        15       1      41      38      39      40                
CQUAD4        16       1       5      33      41      40                
CQUAD4        17       1       6       5      40      39                
CQUAD4        18       1      35      36      42      43                
CQUAD4        19       1      38      34      35      43                
CQUAD4        20       1      39      38      43      42                
CQUAD4        21       1     109       6      39      42                
CQUAD4        22       1     109      42      36      37                
CQUAD4        23       1      45      14      55      46                
CQUAD4        24       1      46      55      59      47                
CQUAD4        25       1      47      59      60      68                
CQUAD4        26       1      68      60      61      49                
CQUAD4        27       1      49      61     137     152                
CQUAD4        28       2      61      60      63     137                
CQUAD4        29       2      60      59      64      63                
CQUAD4        30       2      59      55      65      64                
CQUAD4        31       2      14      13      65      55                
CQUAD4        32       1      45      46      73      72                
CTRIA3        33       1      73      46      74                        
CQUAD4        34       1      46      47      75      74                
CQUAD4        35       1      47      68      76      75                
CQUAD4        36       1      77      76      84      81                
CQUAD4        37       1      81      84      83      82                
CQUAD4        38       1      76      68      83      84                
CQUAD4        39       1      68      49      82      83                
CQUAD4        40       1      79      78      86      85                
CQUAD4        41       1      77      81      86      78                
CQUAD4        42       1      81      82      85      86                
CQUAD4        43       1      49     152      85      82                
CQUAD4        44       1      85     152     166      79                
CQUAD4        45       1      88      87      98      89                
CQUAD4        46       1      89      98     102      90                
CQUAD4        47       1      90     102      96     111                
CQUAD4        48       1     111      96      95      92                
CQUAD4        49       1      92      95      19     109                
CQUAD4        50       2      95      96     106      19                
CQUAD4        51       2      96     102      21     106                
CQUAD4        52       2     102      98      22      21                
CQUAD4        53       2      87      13      22      98                
CQUAD4        54       1      88      89     116     115                
CTRIA3        55       1     116      89     117                        
CQUAD4        56       1      89      90     118     117                
CQUAD4        57       1      90     111     119     118                
CQUAD4        58       1     120     119     127     124                
CQUAD4        59       1     124     127     126     125                
CQUAD4        60       1     119     111     126     127                
CQUAD4        61       1     111      92     125     126                
CQUAD4        62       1     122     121     129     128                
CQUAD4        63       1     120     124     129     121                
CQUAD4        64       1     124     125     128     129                
CQUAD4        65       1      92     109     128     125                
CQUAD4        66       1     128     109      37     122                
CQUAD4        67       1      87      88     156     141                
CQUAD4        68       1     141     156     155     140                
CQUAD4        69       1     140     155     154     139                
CQUAD4        70       1     139     154     135     138                
CQUAD4        71       1     138     135     152     137                
CQUAD4        72       2     139     138     137      63                
CQUAD4        73       2     140     139      63      64                
CQUAD4        74       2     141     140      64      65                
CQUAD4        75       2      13      87     141      65                
CQUAD4        76       1     156      88     115     159                
CTRIA3        77       1     156     159     160                        
CQUAD4        78       1     155     156     160     161                
CQUAD4        79       1     154     155     161     162                
CQUAD4        80       1     162     163     167     170                
CQUAD4        81       1     170     167     168     169                
CQUAD4        82       1     154     162     170     169                
CQUAD4        83       1     135     154     169     168                
CQUAD4        84       1     164     165     171     172                
CQUAD4        85       1     167     163     164     172                
CQUAD4        86       1     168     167     172     171                
CQUAD4        87       1     152     135     168     171                
CQUAD4        88       1     152     171     165     166                
ENDDATA
