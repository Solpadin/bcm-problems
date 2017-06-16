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
$   Date       : Thu Oct 19 01:53:41 2006
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
$ FEMAP Load Set 100 : Surface1  MARKER1
$ FEMAP Load Set 200 : Surface2  MARKER2
$ FEMAP Property 1 : Untitled
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Cylinder
PSOLID         2       1       0        
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           2       0-1.84776-0.76537      1.       0        
GRID           4       0-0.76537-1.84776      1.       0        
GRID           5       0      0.     -2.      1.       0        
GRID           7       0 1.41421-1.41421      1.       0        
GRID           8       0 1.84776-0.76537      1.       0        
GRID           9       0      2.      0.      1.       0        
GRID          11       0 1.41421 1.41421      1.       0        
GRID          12       0 0.76537 1.84776      1.       0        
GRID          13       0      0.      2.      1.       0        
GRID          15       0-1.41421 1.41421      1.       0        
GRID          18       0 0.70711-0.70711      1.       0        
GRID          20       0-0.70711-0.70711      1.       0        
GRID          21       0     -1.      0.      1.       0        
GRID          22       0-0.70711 0.70711      1.       0        
GRID          25       0     -2.      0.      1.       0        
GRID          27       0-1.41421-1.41421      1.       0        
GRID          30       0 0.76537-1.84776      1.       0        
GRID          34       0 1.84776 0.76537      1.       0        
GRID          38       0-0.76537 1.84776      1.       0        
GRID          40       0-1.84776 0.76537      1.       0        
GRID          41       0      1.      0.      1.       0        
GRID          43       0      0.     -1.      1.       0        
GRID          47       0      0.      1.      1.       0        
GRID          48       0 0.70711 0.70711      1.       0        
GRID          49       0     -2.      0.     -1.       0        
GRID          50       0-1.84776-0.76537     -1.       0        
GRID          51       0-1.41421-1.41421     -1.       0        
GRID          52       0-0.76537-1.84776     -1.       0        
GRID          53       0      0.     -2.     -1.       0        
GRID          54       0 0.76537-1.84776     -1.       0        
GRID          55       0 1.41421-1.41421     -1.       0        
GRID          56       0 1.84776-0.76537     -1.       0        
GRID          57       0      2.      0.     -1.       0        
GRID          58       0 1.84776 0.76537     -1.       0        
GRID          59       0 1.41421 1.41421     -1.       0        
GRID          60       0 0.76537 1.84776     -1.       0        
GRID          61       0      0.      2.     -1.       0        
GRID          62       0-0.76537 1.84776     -1.       0        
GRID          63       0-1.41421 1.41421     -1.       0        
GRID          64       0-1.84776 0.76537     -1.       0        
GRID          65       0      1.      0.     -1.       0        
GRID          66       0 0.70711-0.70711     -1.       0        
GRID          67       0      0.     -1.     -1.       0        
GRID          68       0-0.70711-0.70711     -1.       0        
GRID          69       0     -1.      0.     -1.       0        
GRID          70       0-0.70711 0.70711     -1.       0        
GRID          71       0      0.      1.     -1.       0        
GRID          72       0 0.70711 0.70711     -1.       0        
CPENTA        25       2      67      53      54      43       5      30        
CPENTA        26       2      66      67      54      18      43      30        
CPENTA        27       2      66      54      55      18      30       7        
CPENTA        28       2      66      55      56      18       7       8        
CPENTA        29       2      65      66      56      41      18       8        
CPENTA        30       2      65      56      57      41       8       9        
CPENTA        31       2      65      57      58      41       9      34        
CPENTA        32       2      72      65      58      48      41      34        
CPENTA        33       2      72      58      59      48      34      11        
CPENTA        34       2      72      59      60      48      11      12        
CPENTA        35       2      71      72      60      47      48      12        
CPENTA        36       2      71      60      61      47      12      13        
CPENTA        37       2      52      53      67       4       5      43        
CPENTA        38       2      52      67      68       4      43      20        
CPENTA        39       2      51      52      68      27       4      20        
CPENTA        40       2      50      51      68       2      27      20        
CPENTA        41       2      50      68      69       2      20      21        
CPENTA        42       2      49      50      69      25       2      21        
CPENTA        43       2      71      61      62      47      13      38        
CPENTA        44       2      70      71      62      22      47      38        
CPENTA        45       2      70      62      63      22      38      15        
CPENTA        46       2      70      63      64      22      15      40        
CPENTA        47       2      69      70      64      21      22      40        
CPENTA        48       2      49      69      64      25      21      40        
ENDDATA
