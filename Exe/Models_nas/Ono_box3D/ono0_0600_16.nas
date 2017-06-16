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
$   Date       : Sat Mar 29 11:25:24 2008
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
PLOAD4         1      77    111.                             139      78
$ FEMAP Load Set 2 : Endwall Absorb 1.0 0.0
PLOAD4         2      73    222.                              25     122
PLOAD4         2      80    222.                             122      82
PLOAD4         2      84    222.                             122     143
PLOAD4         2      85    222.                             122     188
$ FEMAP Property 1 : ONO BOX
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : INCLUSION
PSHELL         2       1      0.       1               1              0.
$ FEMAP Property 3 : Untitled
PSOLID         3       1       0        
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID          16       0      0.      0.   0.625       0        
GRID          18       0      0.      0.   1.875       0        
GRID          19       0      0.      0.      0.       0        
GRID          21       0    0.35      0.   0.625       0        
GRID          23       0    0.35      0.   1.875       0        
GRID          25       0      0.      0.     2.5       0        
GRID          27       0      0.      0.    1.25       0        
GRID          39       0      0.    0.25     2.5       0        
GRID          43       0    0.35    0.25   0.625       0        
GRID          48       0      0.    0.25   1.875       0        
GRID          53       0    0.35     0.5   0.625       0        
GRID          58       0    0.35    0.25   1.875       0        
GRID          62       0     0.6    0.25    1.25       0        
GRID          64       0     0.6    0.25     2.5       0        
GRID          78       0     0.6     0.5      0.       0        
GRID          80       0     0.6     0.5    1.25       0        
GRID          81       0     0.6     0.5   1.875       0        
GRID          82       0     0.6     0.5     2.5       0        
GRID          84       0     0.6    0.25   1.875       0        
GRID          86       0     0.6    0.25   0.625       0        
GRID          91       0     0.6     0.5   0.625       0        
GRID         108       0     0.6    0.25      0.       0        
GRID         116       0     0.6      0.   0.625       0        
GRID         117       0    0.35      0.      0.       0        
GRID         121       0    0.35      0.     2.5       0        
GRID         122       0    0.35    0.25     2.5       0        
GRID         129       0    0.35      0.    1.25       0        
GRID         134       0     0.6      0.    1.25       0        
GRID         135       0     0.6      0.   1.875       0        
GRID         137       0     0.6      0.      0.       0        
GRID         139       0    0.35    0.25      0.       0        
GRID         143       0     0.6      0.     2.5       0        
GRID         145       0    0.35     0.5      0.       0        
GRID         146       0      0.     0.5      0.       0        
GRID         148       0      0.     0.5    1.25       0        
GRID         151       0    0.35     0.5     2.5       0        
GRID         152       0    0.35     0.5   1.875       0        
GRID         153       0    0.35     0.5    1.25       0        
GRID         158       0    0.35    0.25    1.25       0        
GRID         169       0      0.    0.25   0.625       0        
GRID         183       0      0.    0.25      0.       0        
GRID         185       0      0.    0.25    1.25       0        
GRID         188       0      0.     0.5     2.5       0        
GRID         189       0      0.     0.5   1.875       0        
GRID         191       0      0.     0.5   0.625       0        
CHEXA         73       3      25      18      48      39     121      23+EL   21
+EL   21      58     122                                                        
CHEXA         74       3      18      27     185      48      23     129+EL   22
+EL   22     158      58                                                        
CHEXA         75       3      27      16     169     185     129      21+EL   23
+EL   23      43     158                                                        
CHEXA         76       3      16      19     183     169      21     117+EL   24
+EL   24     139      43                                                        
CHEXA         77       3     139      43      86     108     145      53+EL   25
+EL   25      91      78                                                        
CHEXA         78       3      43     158      62      86      53     153+EL   26
+EL   26      80      91                                                        
CHEXA         79       3     158      58      84      62     153     152+EL   27
+EL   27      81      80                                                        
CHEXA         80       3      58     122      64      84     152     151+EL   28
+EL   28      82      81                                                        
CHEXA         81       3     117     139      43      21     137     108+EL   29
+EL   29      86     116                                                        
CHEXA         82       3      21      43     158     129     116      86+EL   2A
+EL   2A      62     134                                                        
CHEXA         83       3     129     158      58      23     134      62+EL   2B
+EL   2B      84     135                                                        
CHEXA         84       3      23      58     122     121     135      84+EL   2C
+EL   2C      64     143                                                        
CHEXA         85       3     122      58      48      39     151     152+EL   2D
+EL   2D     189     188                                                        
CHEXA         86       3      58     158     185      48     152     153+EL   2E
+EL   2E     148     189                                                        
CHEXA         87       3     158      43     169     185     153      53+EL   2F
+EL   2F     191     148                                                        
CHEXA         88       3      43     139     183     169      53     145+EL   2G
+EL   2G     146     191                                                        
ENDDATA
