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
$   Date       : Sat Mar 29 09:56:56 2008
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
PLOAD4         1     135    111.                             154        
PLOAD4         1     136    111.                             154     146
$ FEMAP Load Set 2 : Endwall Absorb 1.0 0.0
PLOAD4         2     102    222.                              27     161
PLOAD4         2     103    222.                               1     173
PLOAD4         2     104    222.                              64       1
PLOAD4         2     124    222.                             161      48
PLOAD4         2     123    222.                             161        
PLOAD4         2     124    222.                             161      48
$ FEMAP Property 1 : Ono box
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Inclusion
PSHELL         2       1      0.       1               1              0.
$ FEMAP Property 3 : Untitled
PSOLID         3       1       0        
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           1       0      0.      0.     2.5       0        
GRID           6       0      0.     0.5 1.07143       0        
GRID          12       0      0.      0. 0.71429       0        
GRID          13       0      0.      0. 1.07143       0        
GRID          15       0      0.      0. 1.78571       0        
GRID          24       0     0.6      0. 1.78571       0        
GRID          27       0     0.3      0.     2.5       0        
GRID          29       0      0.      0. 2.14286       0        
GRID          31       0      0.      0. 1.42857       0        
GRID          34       0      0.      0. 0.35714       0        
GRID          35       0     0.3      0. 0.35714       0        
GRID          36       0     0.3      0. 0.71429       0        
GRID          37       0     0.3      0. 1.07143       0        
GRID          38       0     0.3      0. 1.42857       0        
GRID          39       0     0.3      0. 1.78571       0        
GRID          40       0     0.3      0. 2.14286       0        
GRID          42       0      0.     0.5     2.5       0        
GRID          45       0     0.6      0.     2.5       0        
GRID          48       0    0.35    0.33     2.5       0        
GRID          55       0     0.3      0.      0.       0        
GRID          56       0      0.      0.      0.       0        
GRID          57       0    0.35     0.5      0.       0        
GRID          59       0    0.35     0.5 0.71429       0        
GRID          61       0    0.35     0.5 1.42857       0        
GRID          64       0    0.35     0.5     2.5       0        
GRID          87       0    0.43    0.25 2.14286       0        
GRID          93       0    0.43    0.25 1.42857       0        
GRID          94       0    0.43    0.25 1.07143       0        
GRID          96       0    0.43    0.25 0.35714       0        
GRID         100       0    0.35    0.33 0.71429       0        
GRID         101       0    0.35    0.33 1.07143       0        
GRID         105       0     0.6      0.      0.       0        
GRID         108       0     0.6    0.25 0.71429       0        
GRID         109       0     0.6    0.25 1.07143       0        
GRID         115       0     0.6      0. 2.14286       0        
GRID         117       0     0.6      0. 1.42857       0        
GRID         118       0     0.6      0. 1.07143       0        
GRID         119       0     0.6      0. 0.71429       0        
GRID         120       0     0.6      0. 0.35714       0        
GRID         126       0    0.35     0.5 1.07143       0        
GRID         128       0    0.35     0.5 0.35714       0        
GRID         130       0      0.     0.5      0.       0        
GRID         131       0      0.     0.5 0.35714       0        
GRID         132       0      0.     0.5 0.71429       0        
GRID         134       0      0.     0.5 1.42857       0        
GRID         135       0      0.     0.5 1.78571       0        
GRID         136       0      0.     0.5 2.14286       0        
GRID         139       0    0.35     0.5 2.14286       0        
GRID         140       0    0.35     0.5 1.78571       0        
GRID         146       0    0.35    0.33      0.       0        
GRID         147       0    0.35    0.33 0.35714       0        
GRID         150       0    0.35    0.33 1.42857       0        
GRID         154       0     0.6    0.25      0.       0        
GRID         155       0     0.6    0.25 0.35714       0        
GRID         159       0     0.6    0.25 1.78571       0        
GRID         161       0     0.6    0.25     2.5       0        
GRID         164       0    0.43    0.25 1.78571       0        
GRID         167       0    0.43    0.25 0.71429       0        
GRID         170       0     0.6     0.5     2.5       0        
GRID         173       0    0.43    0.25     2.5       0        
GRID         177       0    0.43    0.25      0.       0        
GRID         181       0    0.35    0.33 2.14286       0        
GRID         182       0    0.35    0.33 1.78571       0        
GRID         200       0     0.6     0.5 1.42857       0        
GRID         205       0     0.6    0.25 2.14286       0        
GRID         207       0     0.6    0.25 1.42857       0        
GRID         213       0     0.6     0.5 2.14286       0        
GRID         214       0     0.6     0.5 1.78571       0        
GRID         216       0     0.6     0.5 1.07143       0        
GRID         217       0     0.6     0.5 0.71429       0        
GRID         218       0     0.6     0.5 0.35714       0        
GRID         219       0     0.6     0.5      0.       0        
CHEXA        102       3      27     173     161      45      40      87+EL   2U
+EL   2U     205     115                                                        
CHEXA        103       3       1      48     173      27      29     181+EL   2V
+EL   2V      87      40                                                        
CHEXA        104       3      64      48       1      42     139     181+EL   2W
+EL   2W      29     136                                                        
CHEXA        105       3      40      87     205     115      39     164+EL   2X
+EL   2X     159      24                                                        
CHEXA        106       3      29     181      87      40      15     182+EL   2Y
+EL   2Y     164      39                                                        
CHEXA        107       3     139     181      29     136     140     182+EL   2Z
+EL   2Z      15     135                                                        
CHEXA        108       3      39     164     159      24      38      93+EL   30
+EL   30     207     117                                                        
CHEXA        109       3      15     182     164      39      31     150+EL   31
+EL   31      93      38                                                        
CHEXA        110       3     140     182      15     135      61     150+EL   32
+EL   32      31     134                                                        
CHEXA        111       3      38      93     207     117      37      94+EL   33
+EL   33     109     118                                                        
CHEXA        112       3      31     150      93      38      13     101+EL   34
+EL   34      94      37                                                        
CHEXA        113       3      61     150      31     134     126     101+EL   35
+EL   35      13       6                                                        
CHEXA        114       3      37      94     109     118      36     167+EL   36
+EL   36     108     119                                                        
CHEXA        115       3      13     101      94      37      12     100+EL   37
+EL   37     167      36                                                        
CHEXA        116       3     126     101      13       6      59     100+EL   38
+EL   38      12     132                                                        
CHEXA        117       3      36     167     108     119      35      96+EL   39
+EL   39     155     120                                                        
CHEXA        118       3      12     100     167      36      34     147+EL   3A
+EL   3A      96      35                                                        
CHEXA        119       3      59     100      12     132     128     147+EL   3B
+EL   3B      34     131                                                        
CHEXA        120       3      35      96     155     120      55     177+EL   3C
+EL   3C     154     105                                                        
CHEXA        121       3      34     147      96      35      56     146+EL   3D
+EL   3D     177      55                                                        
CHEXA        122       3     128     147      34     131      57     146+EL   3E
+EL   3E      56     130                                                        
CPENTA       123       3     161      64     170     205     139     213        
CHEXA        124       3     161     173      48      64     205      87+EL   3G
+EL   3G     181     139                                                        
CPENTA       125       3     205     139     213     159     140     214        
CHEXA        126       3     205      87     181     139     159     164+EL   3I
+EL   3I     182     140                                                        
CPENTA       127       3     159     140     214     207      61     200        
CHEXA        128       3     159     164     182     140     207      93+EL   3K
+EL   3K     150      61                                                        
CPENTA       129       3     207      61     200     109     126     216        
CHEXA        130       3     207      93     150      61     109      94+EL   3M
+EL   3M     101     126                                                        
CPENTA       131       3     109     126     216     108      59     217        
CHEXA        132       3     109      94     101     126     108     167+EL   3O
+EL   3O     100      59                                                        
CPENTA       133       3     108      59     217     155     128     218        
CHEXA        134       3     108     167     100      59     155      96+EL   3Q
+EL   3Q     147     128                                                        
CPENTA       135       3     155     128     218     154      57     219        
CHEXA        136       3     155      96     147     128     154     177+EL   3S
+EL   3S     146      57                                                        
ENDDATA
