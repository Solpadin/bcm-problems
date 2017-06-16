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
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\Cylinder\Solid\cylinder.MOD
$   Date       : Tue Dec 12 16:14:44 2006
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
$ FEMAP Property 1 : Cylinder
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Inclusion
PSHELL         2       1      0.       1               1              0.
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           1       0   -0.55      0.     0.0       0        
GRID           2       0-0.52308-0.16996     0.0       0        
GRID           3       0-0.44496-0.32328     0.0       0        
GRID           4       0-0.32328-0.44496     0.0       0        
GRID           5       0-0.16996-0.52308     0.0       0        
GRID           6       0      0.   -0.55     0.0       0        
GRID           7       0 0.16996-0.52308     0.0       0        
GRID           8       0 0.32328-0.44496     0.0       0        
GRID           9       0 0.44496-0.32328     0.0       0        
GRID          10       0 0.52308-0.16996     0.0       0        
GRID          11       0    0.55      0.     0.0       0        
GRID          12       0 0.52308 0.16996     0.0       0        
GRID          13       0 0.44496 0.32328     0.0       0        
GRID          14       0 0.32328 0.44496     0.0       0        
GRID          15       0 0.16996 0.52308     0.0       0        
GRID          16       0      0.    0.55     0.0       0        
GRID          17       0-0.16996 0.52308     0.0       0        
GRID          18       0-0.32328 0.44496     0.0       0        
GRID          19       0-0.44496 0.32328     0.0       0        
GRID          20       0-0.52308 0.16996     0.0       0        
GRID          26       0   -0.25      0.     0.0       0        
GRID          30       0 0.20225 0.14695     0.0       0        
GRID          31       0  0.35358.973E-5     0.0       0        
GRID          32       0-0.381360.027721     0.0       0        
GRID          33       00.091212 -0.3699     0.0       0        
GRID          34       0 0.28803-0.24185     0.0       0        
GRID          35       0-0.32403-0.20164     0.0       0        
GRID          36       0-0.14358-0.35354     0.0       0        
GRID          37       00.093306 0.37611     0.0       0        
GRID          38       0 0.28832 0.24267     0.0       0        
GRID          39       0-0.29729 0.24956     0.0       0        
GRID          40       0-0.12902 0.39678     0.0       0        
GRID          42       0-0.20225-0.14695     0.0       0        
GRID          43       0-7.73E-2-0.23776     0.0       0        
GRID          44       00.077254-0.23776     0.0       0        
GRID          45       0 0.20225-0.14695     0.0       0        
GRID          46       0    0.25      0.     0.0       0        
GRID          48       00.077254 0.23776     0.0       0        
GRID          49       0-7.73E-2 0.23776     0.0       0        
GRID          50       0-0.20225 0.14695     0.0       0        
GRID          51       00.093408      0.     0.0       0        
GRID          52       0-9.37E-2      0.     0.0       0        
GRID          53       0-5.63E-5-0.12824     0.0       0        
GRID          54       0-5.63E-5 0.12824     0.0       0        
CTRIA3         1       1      31      46      45                        
CTRIA3         2       1      31      45      34                        
CTRIA3         3       1      10      11      31                        
CTRIA3         4       1      10      31      34                        
CTRIA3         5       1       9      10      34                        
CTRIA3         6       1      34      45      44                        
CTRIA3         7       1      34      44      33                        
CTRIA3         8       1       8       9      34                        
CTRIA3         9       1       8      34      33                        
CTRIA3        10       1       7       8      33                        
CTRIA3        11       1      33      44      43                        
CTRIA3        12       1      33      43      36                        
CTRIA3        13       1       6       7      33                        
CTRIA3        14       1       6      33      36                        
CTRIA3        15       1       5       6      36                        
CTRIA3        16       1      36      43      42                        
CTRIA3        17       1      36      42      35                        
CTRIA3        18       1       4       5      36                        
CTRIA3        19       1       4      36      35                        
CTRIA3        20       1       3       4      35                        
CTRIA3        21       1      35      42      26                        
CTRIA3        22       1      35      26      32                        
CTRIA3        23       1       2       3      35                        
CTRIA3        24       1       2      35      32                        
CTRIA3        25       1       1       2      32                        
CTRIA3        26       1      31      11      12                        
CTRIA3        27       1      12      13      38                        
CTRIA3        28       1      31      12      38                        
CTRIA3        29       1      30      46      31                        
CTRIA3        30       1      30      31      38                        
CTRIA3        31       1      38      13      14                        
CTRIA3        32       1      14      15      37                        
CTRIA3        33       1      38      14      37                        
CTRIA3        34       1      48      30      38                        
CTRIA3        35       1      48      38      37                        
CTRIA3        36       1      37      15      16                        
CTRIA3        37       1      16      17      40                        
CTRIA3        38       1      37      16      40                        
CTRIA3        39       1      49      48      37                        
CTRIA3        40       1      49      37      40                        
CTRIA3        41       1      18      19      39                        
CTRIA3        42       1      40      17      18                        
CTRIA3        43       1      40      18      39                        
CTRIA3        44       1      49      40      39                        
CTRIA3        45       1      50      49      39                        
CTRIA3        46       1      32      26      50                        
CTRIA3        47       1      32      50      39                        
CTRIA3        48       1      39      19      20                        
CTRIA3        49       1      32      39      20                        
CTRIA3        50       1       1      32      20                        
CTRIA3        51       2      51      52      53                        
CTRIA3        52       2      45      46      51                        
CTRIA3        53       2      45      51      53                        
CTRIA3        54       2      44      45      53                        
CTRIA3        55       2      43      44      53                        
CTRIA3        56       2      42      43      53                        
CTRIA3        57       2      42      53      52                        
CTRIA3        58       2      26      42      52                        
CTRIA3        59       2      51      46      30                        
CTRIA3        60       2      30      48      54                        
CTRIA3        61       2      51      30      54                        
CTRIA3        62       2      52      51      54                        
CTRIA3        63       2      54      48      49                        
CTRIA3        64       2      54      49      50                        
CTRIA3        65       2      52      54      50                        
CTRIA3        66       2      26      52      50                        
ENDDATA
