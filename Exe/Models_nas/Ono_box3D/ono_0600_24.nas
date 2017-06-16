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
$   Date       : Sat Mar 29 09:55:37 2008
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
PLOAD4         1      83    111.                             112        
PLOAD4         1      84    111.                             112      61
$ FEMAP Load Set 2 : Endwall Absorb 1.0 0.0
PLOAD4         2      61    222.                              71        
PLOAD4         2      62    222.                               1        
PLOAD4         2      63    222.                              76       1
PLOAD4         2      64    222.                               1        
PLOAD4         2      78    222.                              25      27
PLOAD4         2      77    222.                              25        
PLOAD4         2      78    222.                              25      27
$ FEMAP Property 1 : Ono box
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Inclusion
PSHELL         2       1      0.       1               1              0.
$ FEMAP Property 3 : Untitled
PSOLID         3       1       0        
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           1       0      0.      0.     2.5       0        
GRID           4       0      0.     0.5    1.25       0        
GRID           5       0      0.     0.5   0.625       0        
GRID           6       0      0.     0.5      0.       0        
GRID          11       0      0.      0.      0.       0        
GRID          14       0     0.6      0.    1.25       0        
GRID          15       0     0.6      0.   1.875       0        
GRID          18       0      0.      0.   1.875       0        
GRID          19       0      0.      0.    1.25       0        
GRID          20       0      0.      0.   0.625       0        
GRID          22       0      0.     0.5     2.5       0        
GRID          25       0     0.6    0.25     2.5       0        
GRID          27       0    0.35    0.33     2.5       0        
GRID          33       0     0.6      0.      0.       0        
GRID          36       0    0.35     0.5   0.625       0        
GRID          41       0    0.35    0.33   1.875       0        
GRID          46       0     0.6    0.25   1.875       0        
GRID          48       0     0.6    0.25   0.625       0        
GRID          51       0    0.43    0.25   0.625       0        
GRID          52       0    0.43    0.25    1.25       0        
GRID          57       0    0.43    0.25   1.875       0        
GRID          61       0    0.35    0.33      0.       0        
GRID          63       0    0.35    0.33    1.25       0        
GRID          71       0     0.6      0.     2.5       0        
GRID          74       0     0.6      0.   0.625       0        
GRID          76       0    0.35     0.5     2.5       0        
GRID          77       0    0.35     0.5   1.875       0        
GRID          78       0    0.35     0.5    1.25       0        
GRID          84       0      0.     0.5   1.875       0        
GRID         109       0    0.43    0.25     2.5       0        
GRID         110       0    0.35     0.5      0.       0        
GRID         112       0     0.6    0.25      0.       0        
GRID         113       0    0.43    0.25      0.       0        
GRID         119       0    0.35    0.33   0.625       0        
GRID         126       0     0.6     0.5      0.       0        
GRID         128       0     0.6     0.5    1.25       0        
GRID         133       0     0.6    0.25    1.25       0        
GRID         136       0     0.6     0.5     2.5       0        
GRID         137       0     0.6     0.5   1.875       0        
GRID         139       0     0.6     0.5   0.625       0        
CPENTA        61       3      71     109      25      15      57      46        
CPENTA        62       3       1     109      71      18      57      15        
CHEXA         63       3      76      27       1      22      77      41+EL   1R
+EL   1R      18      84                                                        
CPENTA        64       3       1      27     109      18      41      57        
CPENTA        65       3      15      57      46      14      52     133        
CPENTA        66       3      18      57      15      19      52      14        
CHEXA         67       3      77      41      18      84      78      63+EL   1V
+EL   1V      19       4                                                        
CPENTA        68       3      18      41      57      19      63      52        
CPENTA        69       3      14      52     133      74      51      48        
CPENTA        70       3      19      52      14      20      51      74        
CHEXA         71       3      78      63      19       4      36     119+EL   1Z
+EL   1Z      20       5                                                        
CPENTA        72       3      19      63      52      20     119      51        
CPENTA        73       3      74      51      48      33     113     112        
CPENTA        74       3      20      51      74      11     113      33        
CHEXA         75       3      36     119      20       5     110      61+EL   23
+EL   23      11       6                                                        
CPENTA        76       3      20     119      51      11      61     113        
CPENTA        77       3      25      76     136      46      77     137        
CHEXA         78       3      25     109      27      76      46      57+EL   26
+EL   26      41      77                                                        
CPENTA        79       3      46      77     137     133      78     128        
CHEXA         80       3      46      57      41      77     133      52+EL   28
+EL   28      63      78                                                        
CPENTA        81       3     133      78     128      48      36     139        
CHEXA         82       3     133      52      63      78      48      51+EL   2A
+EL   2A     119      36                                                        
CPENTA        83       3      48      36     139     112     110     126        
CHEXA         84       3      48      51     119      36     112     113+EL   2C
+EL   2C      61     110                                                        
ENDDATA
