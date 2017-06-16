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
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\GeoNTS\Model\model_layer.MOD
$   Date       : Thu Mar 27 19:05:47 2008
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
$ FEMAP Property 1 : Phase 1
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Phase 2
PSHELL         2       1      0.       1               1              0.
$ FEMAP Material 1 : Abstract material
MAT1           1                              0.      0.      0.        
GRID           1       0      0.      0.      0.       0        
GRID           4       0      0.      1.      0.       0        
GRID           5       0    0.05      0.      0.       0        
GRID           8       0    0.05      1.      0.       0        
GRID           9       0     0.1      0.      0.       0        
GRID          10       0     0.1      1.      0.       0        
GRID          11       0    0.15      0.      0.       0        
GRID          12       0     0.2      0.      0.       0        
GRID          14       0    0.15      1.      0.       0        
GRID          16       0     0.2      1.      0.       0        
GRID          17       0    0.25      0.      0.       0        
GRID          18       0     0.3      0.      0.       0        
GRID          19       0     0.3      1.      0.       0        
GRID          20       0    0.25      1.      0.       0        
GRID          23       0    0.35      0.      0.       0        
GRID          25       0     0.4      1.      0.       0        
GRID          26       0    0.35      1.      0.       0        
GRID          27       0     0.4      0.      0.       0        
GRID          29       0    0.45      0.      0.       0        
GRID          30       0     0.5      0.      0.       0        
GRID          31       0     0.5      1.      0.       0        
GRID          32       0    0.45      1.      0.       0        
GRID          35       0    0.55      0.      0.       0        
GRID          38       0    0.55      1.      0.       0        
GRID          39       0     0.6      0.      0.       0        
GRID          40       0     0.6      1.      0.       0        
GRID          41       0    0.65      0.      0.       0        
GRID          42       0     0.7      0.      0.       0        
GRID          44       0    0.65      1.      0.       0        
GRID          46       0     0.7      1.      0.       0        
GRID          47       0    0.75      0.      0.       0        
GRID          49       0     0.8      1.      0.       0        
GRID          50       0    0.75      1.      0.       0        
GRID          51       0     0.8      0.      0.       0        
GRID          53       0    0.85      0.      0.       0        
GRID          55       0     0.9      1.      0.       0        
GRID          56       0    0.85      1.      0.       0        
GRID          57       0     0.9      0.      0.       0        
GRID          59       0    0.95      0.      0.       0        
GRID          62       0    0.95      1.      0.       0        
GRID          63       0      1.      0.      0.       0        
GRID          64       0      1.      1.      0.       0        
GRID          65       0    1.05      0.      0.       0        
GRID          66       0     1.1      0.      0.       0        
GRID          67       0     1.1      1.      0.       0        
GRID          68       0    1.05      1.      0.       0        
GRID          71       0    1.15      0.      0.       0        
GRID          74       0    1.15      1.      0.       0        
GRID          75       0     1.2      0.      0.       0        
GRID          76       0     1.2      1.      0.       0        
GRID          77       0    1.25      0.      0.       0        
GRID          80       0    1.25      1.      0.       0        
GRID          81       0     1.3      0.      0.       0        
GRID          82       0     1.3      1.      0.       0        
GRID          83       0    1.35      0.      0.       0        
GRID          84       0     1.4      0.      0.       0        
GRID          86       0    1.35      1.      0.       0        
GRID          88       0     1.4      1.      0.       0        
GRID          89       0    1.45      0.      0.       0        
GRID          90       0     1.5      0.      0.       0        
GRID          92       0    1.45      1.      0.       0        
GRID          94       0     1.5      1.      0.       0        
GRID          95       0    1.55      0.      0.       0        
GRID          97       0     1.6      1.      0.       0        
GRID          98       0    1.55      1.      0.       0        
GRID          99       0     1.6      0.      0.       0        
GRID         101       0    1.65      0.      0.       0        
GRID         104       0    1.65      1.      0.       0        
GRID         105       0     1.7      0.      0.       0        
GRID         106       0     1.7      1.      0.       0        
GRID         107       0    1.75      0.      0.       0        
GRID         108       0     1.8      0.      0.       0        
GRID         109       0     1.8      1.      0.       0        
GRID         110       0    1.75      1.      0.       0        
GRID         113       0    1.85      0.      0.       0        
GRID         115       0     1.9      1.      0.       0        
GRID         116       0    1.85      1.      0.       0        
GRID         117       0     1.9      0.      0.       0        
GRID         119       0    1.95      0.      0.       0        
GRID         120       0      2.      0.      0.       0        
GRID         121       0      2.      1.      0.       0        
GRID         122       0    1.95      1.      0.       0        
CQUAD4         1       1       1       5       8       4                
CQUAD4         2       2       5       9      10       8                
CQUAD4         3       1       9      11      14      10                
CQUAD4         4       2      11      12      16      14                
CQUAD4         5       1      12      17      20      16                
CQUAD4         6       2      17      18      19      20                
CQUAD4         7       1      18      23      26      19                
CQUAD4         8       2      23      27      25      26                
CQUAD4         9       1      27      29      32      25                
CQUAD4        10       2      29      30      31      32                
CQUAD4        11       1      30      35      38      31                
CQUAD4        12       2      35      39      40      38                
CQUAD4        13       1      39      41      44      40                
CQUAD4        14       2      41      42      46      44                
CQUAD4        15       1      42      47      50      46                
CQUAD4        16       2      47      51      49      50                
CQUAD4        17       1      51      53      56      49                
CQUAD4        18       2      53      57      55      56                
CQUAD4        19       1      57      59      62      55                
CQUAD4        20       2      59      63      64      62                
CQUAD4        21       1      63      65      68      64                
CQUAD4        22       2      65      66      67      68                
CQUAD4        23       1      66      71      74      67                
CQUAD4        24       2      71      75      76      74                
CQUAD4        25       1      75      77      80      76                
CQUAD4        26       2      77      81      82      80                
CQUAD4        27       1      81      83      86      82                
CQUAD4        28       2      83      84      88      86                
CQUAD4        29       1      84      89      92      88                
CQUAD4        30       2      89      90      94      92                
CQUAD4        31       1      90      95      98      94                
CQUAD4        32       2      95      99      97      98                
CQUAD4        33       1      99     101     104      97                
CQUAD4        34       2     101     105     106     104                
CQUAD4        35       1     105     107     110     106                
CQUAD4        36       2     107     108     109     110                
CQUAD4        37       1     108     113     116     109                
CQUAD4        38       2     113     117     115     116                
CQUAD4        39       1     117     119     122     115                
CQUAD4        40       2     119     120     121     122                
ENDDATA
