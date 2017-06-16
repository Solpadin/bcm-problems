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
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\Ono_box\Solid\ono_box0.MOD
$   Date       : Sat Mar 29 11:21:23 2008
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
PLOAD4         1      60    111.                              64      66
$ FEMAP Load Set 2 : Endwall Absorb 1.0 0.0
PLOAD4         2      57    222.                               4      88
PLOAD4         2      62    222.                              88      74
PLOAD4         2      65    222.                              88     119
PLOAD4         2      66    222.                              88     125
$ FEMAP Property 1 : ONO BOX
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : INCLUSION
PSHELL         2       1      0.       1               1              0.
$ FEMAP Property 3 : Untitled
PSOLID         3       1       0        
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           4       0      0.      0.     2.5       0        
GRID           8       0      0.      0.      0.       0        
GRID          11       0      0.    0.25 1.66667       0        
GRID          12       0      0.    0.25 0.83333       0        
GRID          15       0      0.      0. 0.83333       0        
GRID          23       0      0.      0. 1.66667       0        
GRID          28       0    0.35      0.      0.       0        
GRID          44       0    0.35     0.5 0.83333       0        
GRID          45       0    0.35     0.5      0.       0        
GRID          48       0    0.35    0.25 1.66667       0        
GRID          64       0    0.35    0.25      0.       0        
GRID          65       0     0.6    0.25      0.       0        
GRID          66       0     0.6     0.5      0.       0        
GRID          67       0     0.6     0.5 0.83333       0        
GRID          71       0     0.6    0.25 1.66667       0        
GRID          73       0    0.35     0.5     2.5       0        
GRID          74       0     0.6     0.5     2.5       0        
GRID          75       0     0.6     0.5 1.66667       0        
GRID          81       0     0.6    0.25     2.5       0        
GRID          83       0     0.6    0.25 0.83333       0        
GRID          88       0    0.35    0.25     2.5       0        
GRID          95       0     0.6      0. 1.66667       0        
GRID          98       0    0.35      0. 0.83333       0        
GRID         103       0    0.35    0.25 0.83333       0        
GRID         106       0    0.35      0. 1.66667       0        
GRID         110       0     0.6      0. 0.83333       0        
GRID         113       0     0.6      0.      0.       0        
GRID         118       0    0.35      0.     2.5       0        
GRID         119       0     0.6      0.     2.5       0        
GRID         122       0      0.     0.5      0.       0        
GRID         125       0      0.     0.5     2.5       0        
GRID         136       0    0.35     0.5 1.66667       0        
GRID         138       0      0.    0.25     2.5       0        
GRID         153       0      0.    0.25      0.       0        
GRID         158       0      0.     0.5 1.66667       0        
GRID         159       0      0.     0.5 0.83333       0        
CHEXA         57       3       4      23      11     138     118     106+EL   1L
+EL   1L      48      88                                                        
CHEXA         58       3      23      15      12      11     106      98+EL   1M
+EL   1M     103      48                                                        
CHEXA         59       3      15       8     153      12      98      28+EL   1N
+EL   1N      64     103                                                        
CHEXA         60       3      64     103      83      65      45      44+EL   1O
+EL   1O      67      66                                                        
CHEXA         61       3     103      48      71      83      44     136+EL   1P
+EL   1P      75      67                                                        
CHEXA         62       3      48      88      81      71     136      73+EL   1Q
+EL   1Q      74      75                                                        
CHEXA         63       3      28      64     103      98     113      65+EL   1R
+EL   1R      83     110                                                        
CHEXA         64       3      98     103      48     106     110      83+EL   1S
+EL   1S      71      95                                                        
CHEXA         65       3     106      48      88     118      95      71+EL   1T
+EL   1T      81     119                                                        
CHEXA         66       3      88      48      11     138      73     136+EL   1U
+EL   1U     158     125                                                        
CHEXA         67       3      48     103      12      11     136      44+EL   1V
+EL   1V     159     158                                                        
CHEXA         68       3     103      64     153      12      44      45+EL   1W
+EL   1W     122     159                                                        
ENDDATA
