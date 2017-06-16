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
  LOAD = 1
BEGIN BULK
$ ***************************************************************************
$   Written by : FEMAP
$   Version    : 7.00
$   Translator : MSC/NASTRAN
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\Ono_box\Solid\ono_box.MOD
$   Date       : Sat Mar 29 09:54:28 2008
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
$ FEMAP Load Set 1 : Srcplane Vn -1.0 0.0
PLOAD4         1      65    111.                              80        
PLOAD4         1      66    111.                              80      52
$ FEMAP Load Set 2 : Endwall Absorb 1.0 0.0
PLOAD4         2      49    222.                              60        
PLOAD4         2      50    222.                               1        
PLOAD4         2      51    222.                              89       1
PLOAD4         2      52    222.                               1        
PLOAD4         2      62    222.                              39      71
PLOAD4         2      61    222.                              39        
PLOAD4         2      62    222.                              39      71
$ FEMAP Property 1 : Ono box
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Inclusion
PSHELL         2       1      0.       1               1              0.
$ FEMAP Property 3 : Untitled
PSOLID         3       1       0        
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           1       0      0.      0.     2.5       0        
GRID           3       0      0.     0.5 1.66667       0        
GRID           4       0      0.     0.5 0.83333       0        
GRID           7       0      0.      0. 0.83333       0        
GRID           8       0      0.      0. 1.66667       0        
GRID          11       0     0.6      0. 0.83333       0        
GRID          18       0      0.     0.5     2.5       0        
GRID          24       0      0.     0.5      0.       0        
GRID          30       0      0.      0.      0.       0        
GRID          31       0    0.35     0.5      0.       0        
GRID          39       0     0.6    0.25     2.5       0        
GRID          49       0    0.43    0.25 1.66667       0        
GRID          50       0    0.43    0.25 0.83333       0        
GRID          52       0    0.35    0.33      0.       0        
GRID          55       0     0.6      0.      0.       0        
GRID          58       0     0.6    0.25 1.66667       0        
GRID          60       0     0.6      0.     2.5       0        
GRID          61       0     0.6      0. 1.66667       0        
GRID          66       0    0.35     0.5 0.83333       0        
GRID          71       0    0.35    0.33     2.5       0        
GRID          73       0    0.35     0.5 1.66667       0        
GRID          77       0    0.35    0.33 0.83333       0        
GRID          78       0    0.35    0.33 1.66667       0        
GRID          80       0     0.6    0.25      0.       0        
GRID          88       0     0.6     0.5     2.5       0        
GRID          89       0    0.35     0.5     2.5       0        
GRID          91       0    0.43    0.25     2.5       0        
GRID          93       0     0.6     0.5      0.       0        
GRID         102       0    0.43    0.25      0.       0        
GRID         108       0     0.6     0.5 1.66667       0        
GRID         112       0     0.6    0.25 0.83333       0        
GRID         116       0     0.6     0.5 0.83333       0        
CPENTA        49       3      60      91      39      61      49      58        
CPENTA        50       3       1      91      60       8      49      61        
CHEXA         51       3      89      71       1      18      73      78+EL   1F
+EL   1F       8       3                                                        
CPENTA        52       3       1      71      91       8      78      49        
CPENTA        53       3      61      49      58      11      50     112        
CPENTA        54       3       8      49      61       7      50      11        
CHEXA         55       3      73      78       8       3      66      77+EL   1J
+EL   1J       7       4                                                        
CPENTA        56       3       8      78      49       7      77      50        
CPENTA        57       3      11      50     112      55     102      80        
CPENTA        58       3       7      50      11      30     102      55        
CHEXA         59       3      66      77       7       4      31      52+EL   1N
+EL   1N      30      24                                                        
CPENTA        60       3       7      77      50      30      52     102        
CPENTA        61       3      39      89      88      58      73     108        
CHEXA         62       3      39      91      71      89      58      49+EL   1Q
+EL   1Q      78      73                                                        
CPENTA        63       3      58      73     108     112      66     116        
CHEXA         64       3      58      49      78      73     112      50+EL   1S
+EL   1S      77      66                                                        
CPENTA        65       3     112      66     116      80      31      93        
CHEXA         66       3     112      50      77      66      80     102+EL   1U
+EL   1U      52      31                                                        
ENDDATA
