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
$   Date       : Wed Oct 18 12:32:43 2006
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
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : penta
PSOLID         2       1       0        
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           1       0     -2.      0.      1.       0        
GRID           2       0-1.84776-0.76537      1.       0        
GRID           3       0-1.41421-1.41421      1.       0        
GRID           4       0-0.76537-1.84776      1.       0        
GRID           5       0      0.     -2.      1.       0        
GRID           6       0 0.76537-1.84776      1.       0        
GRID           7       0 1.41421-1.41421      1.       0        
GRID           8       0 1.84776-0.76537      1.       0        
GRID           9       0      2.      0.      1.       0        
GRID          10       0 1.84776 0.76537      1.       0        
GRID          11       0 1.41421 1.41421      1.       0        
GRID          12       0 0.76537 1.84776      1.       0        
GRID          13       0      0.      2.      1.       0        
GRID          14       0-0.76537 1.84776      1.       0        
GRID          15       0-1.41421 1.41421      1.       0        
GRID          16       0-1.84776 0.76537      1.       0        
GRID          17       0      1.      0.      1.       0        
GRID          18       0 0.70711-0.70711      1.       0        
GRID          19       0      0.     -1.      1.       0        
GRID          20       0-0.70711-0.70711      1.       0        
GRID          21       0     -1.      0.      1.       0        
GRID          22       0-0.70711 0.70711      1.       0        
GRID          23       0      0.      1.      1.       0        
GRID          24       0 0.70711 0.70711      1.       0        
GRID          25       0     -2.      0.      1.       0        
GRID          26       0-1.84776-0.76537      1.       0        
GRID          27       0-1.41421-1.41421      1.       0        
GRID          28       0-0.76537-1.84776      1.       0        
GRID          29       0      0.     -2.      1.       0        
GRID          30       0 0.76537-1.84776      1.       0        
GRID          31       0 1.41421-1.41421      1.       0        
GRID          32       0 1.84776-0.76537      1.       0        
GRID          33       0      2.      0.      1.       0        
GRID          34       0 1.84776 0.76537      1.       0        
GRID          35       0 1.41421 1.41421      1.       0        
GRID          36       0 0.76537 1.84776      1.       0        
GRID          37       0      0.      2.      1.       0        
GRID          38       0-0.76537 1.84776      1.       0        
GRID          39       0-1.41421 1.41421      1.       0        
GRID          40       0-1.84776 0.76537      1.       0        
GRID          41       0      1.      0.      1.       0        
GRID          42       0 0.70711-0.70711      1.       0        
GRID          43       0      0.     -1.      1.       0        
GRID          44       0-0.70711-0.70711      1.       0        
GRID          45       0     -1.      0.      1.       0        
GRID          46       0-0.70711 0.70711      1.       0        
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
CPENTA        29       2      65      66      56      41      42      32        
CPENTA        30       2      65      56      57      41      32      33        
CPENTA        31       2      65      57      58      41      33      34        
CPENTA        32       2      72      65      58      48      41      34        
ENDDATA
