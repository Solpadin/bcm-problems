ID FEMAP,FEMAP
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
$   From Model : 
$   Date       : Tue Feb 05 21:51:59 2008
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
PSHELL         1       1      0.       1      1.       1 0.83333      0.
$ FEMAP Property 2 : Untitled
PSHELL         2       1      0.       1      1.       1 0.83333      0.
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           4       0   0.725    0.25      0.       0        
GRID           6       0     1.2      0.      0.       0        
GRID           7       0    1.08      0.      0.       0        
GRID           8       0    0.96      0.      0.       0        
GRID           9       0    0.84      0.      0.       0        
GRID          10       0    0.72      0.      0.       0        
GRID          12       0     0.6   0.125      0.       0        
GRID          13       0     0.6    0.25      0.       0        
GRID          14       0    0.85    0.25      0.       0        
GRID          15       0    0.85   0.375      0.       0        
GRID          20       0     1.2   0.375      0.       0        
GRID          21       0     1.2    0.25      0.       0        
GRID          22       0     1.2   0.125      0.       0        
GRID          23       0 0.96594 0.37562      0.       0        
GRID          24       0  0.9649 0.25072      0.       0        
GRID          25       0 0.96224 0.12557      0.       0        
GRID          26       0 0.84369 0.12494      0.       0        
GRID          27       0 0.72213 0.12495      0.       0        
GRID          28       0 1.08266 0.37562      0.       0        
GRID          29       0 1.08199 0.25088      0.       0        
GRID          30       0 1.08099 0.12546      0.       0        
GRID          34       0   0.725    0.75      0.       0        
GRID          36       0     1.2      1.      0.       0        
GRID          37       0    1.08      1.      0.       0        
GRID          38       0    0.96      1.      0.       0        
GRID          39       0    0.84      1.      0.       0        
GRID          40       0    0.72      1.      0.       0        
GRID          43       0     0.6    0.75      0.       0        
GRID          44       0    0.85    0.75      0.       0        
GRID          45       0    0.85   0.625      0.       0        
GRID          46       0    0.85     0.5      0.       0        
GRID          47       0 0.96667     0.5      0.       0        
GRID          48       0 1.08333     0.5      0.       0        
GRID          49       0     1.2     0.5      0.       0        
GRID          50       0     1.2   0.625      0.       0        
GRID          51       0     1.2    0.75      0.       0        
GRID          52       0     1.2   0.875      0.       0        
GRID          53       0 0.96594 0.62438      0.       0        
GRID          54       0  0.9649 0.74928      0.       0        
GRID          55       0 0.96224 0.87443      0.       0        
GRID          56       0 0.84369 0.87506      0.       0        
GRID          57       0 0.72213 0.87505      0.       0        
GRID          58       0 1.08266 0.62438      0.       0        
GRID          59       0 1.08199 0.74912      0.       0        
GRID          60       0 1.08099 0.87454      0.       0        
GRID          64       0   0.475    0.75      0.       0        
GRID          66       0      0.      1.      0.       0        
GRID          67       0    0.12      1.      0.       0        
GRID          68       0    0.24      1.      0.       0        
GRID          69       0    0.36      1.      0.       0        
GRID          70       0    0.48      1.      0.       0        
GRID          71       0     0.6      1.      0.       0        
GRID          72       0     0.6   0.875      0.       0        
GRID          74       0    0.35    0.75      0.       0        
GRID          75       0    0.35   0.625      0.       0        
GRID          76       0    0.35     0.5      0.       0        
GRID          77       0 0.23333     0.5      0.       0        
GRID          78       0 0.11667     0.5      0.       0        
GRID          80       0      0.   0.625      0.       0        
GRID          81       0      0.    0.75      0.       0        
GRID          82       0      0.   0.875      0.       0        
GRID          83       0 0.23406 0.62438      0.       0        
GRID          84       0  0.2351 0.74928      0.       0        
GRID          85       0 0.23776 0.87443      0.       0        
GRID          86       0 0.35631 0.87506      0.       0        
GRID          87       0 0.47787 0.87505      0.       0        
GRID          88       0 0.11734 0.62438      0.       0        
GRID          89       0 0.11801 0.74912      0.       0        
GRID          90       0 0.11901 0.87454      0.       0        
GRID         124       0   0.475    0.25      0.       0        
GRID         126       0      0.      0.      0.       0        
GRID         127       0    0.12      0.      0.       0        
GRID         128       0    0.24      0.      0.       0        
GRID         129       0    0.36      0.      0.       0        
GRID         130       0    0.48      0.      0.       0        
GRID         131       0     0.6      0.      0.       0        
GRID         135       0    0.35    0.25      0.       0        
GRID         136       0    0.35   0.375      0.       0        
GRID         140       0      0.     0.5      0.       0        
GRID         141       0      0.   0.375      0.       0        
GRID         142       0      0.    0.25      0.       0        
GRID         143       0      0.   0.125      0.       0        
GRID         144       0 0.23406 0.37562      0.       0        
GRID         145       0  0.2351 0.25072      0.       0        
GRID         146       0 0.23776 0.12557      0.       0        
GRID         147       0 0.35631 0.12494      0.       0        
GRID         148       0 0.47787 0.12495      0.       0        
GRID         149       0 0.11734 0.37562      0.       0        
GRID         150       0 0.11801 0.25088      0.       0        
GRID         151       0 0.11901 0.12546      0.       0        
CQUAD4         5       1      13      12      27       4                
CQUAD4         6       1     131      10      27      12                
CQUAD4         7       1       4      27      26      14                
CQUAD4         8       1      10       9      26      27                
CQUAD4         9       1      46      15      23      47                
CQUAD4        10       1      15      14      24      23                
CQUAD4        11       1      14      26      25      24                
CQUAD4        12       1       9       8      25      26                
CQUAD4        13       1      47      23      28      48                
CQUAD4        14       1      23      24      29      28                
CQUAD4        15       1      24      25      30      29                
CQUAD4        16       1       8       7      30      25                
CQUAD4        17       1      48      28      20      49                
CQUAD4        18       1      28      29      21      20                
CQUAD4        19       1      29      30      22      21                
CQUAD4        20       1       7       6      22      30                
CQUAD4        25       1      72      43      34      57                
CQUAD4        26       1      40      71      72      57                
CQUAD4        27       1      57      34      44      56                
CQUAD4        28       1      39      40      57      56                
CQUAD4        29       1      45      46      47      53                
CQUAD4        30       1      44      45      53      54                
CQUAD4        31       1      56      44      54      55                
CQUAD4        32       1      38      39      56      55                
CQUAD4        33       1      53      47      48      58                
CQUAD4        34       1      54      53      58      59                
CQUAD4        35       1      55      54      59      60                
CQUAD4        36       1      37      38      55      60                
CQUAD4        37       1      58      48      49      50                
CQUAD4        38       1      59      58      50      51                
CQUAD4        39       1      60      59      51      52                
CQUAD4        40       1      36      37      60      52                
CQUAD4        45       1      43      72      87      64                
CQUAD4        46       1      71      70      87      72                
CQUAD4        47       1      64      87      86      74                
CQUAD4        48       1      70      69      86      87                
CQUAD4        49       1      76      75      83      77                
CQUAD4        50       1      75      74      84      83                
CQUAD4        55       1      12      13     124     148                
CQUAD4        56       1     130     131      12     148                
CQUAD4        57       1     148     124     135     147                
CQUAD4        58       1     129     130     148     147                
CQUAD4        59       1     136      76      77     144                
CQUAD4        60       1     135     136     144     145                
CQUAD4        61       1     147     135     145     146                
CQUAD4        62       1     128     129     147     146                
CQUAD4        63       1     144      77      78     149                
CQUAD4        64       1     145     144     149     150                
CQUAD4        65       1     146     145     150     151                
CQUAD4        66       1     127     128     146     151                
CQUAD4        67       1     149      78     140     141                
CQUAD4        68       1     150     149     141     142                
CQUAD4        69       1     151     150     142     143                
CQUAD4        70       1     126     127     151     143                
CQUAD4        71       1      74      86      85      84                
CQUAD4        72       1      69      68      85      86                
CQUAD4        73       1      77      83      88      78                
CQUAD4        74       1      83      84      89      88                
CQUAD4        75       1      84      85      90      89                
CQUAD4        76       1      68      67      90      85                
CQUAD4        77       1      78      88      80     140                
CQUAD4        78       1      88      89      81      80                
CQUAD4        79       1      89      90      82      81                
CQUAD4        80       1      67      66      82      90                
ENDDATA
