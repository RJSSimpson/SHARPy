G95 module created on Tue Apr 19 20:21:29 2011 from ../../src/lib_fem.f90
If you edit this, you'll get what you deserve.
module-version 8
(() () () () () () () () () () () () () () () () () () () () ())

()

()

()

()

(2 'fem_1d_gauss_val' 'lib_fem' 1 ((PROCEDURE UNKNOWN MODULE-PROC DECL
NONE NONE SUBROUTINE) (PROCEDURE 0) 0 0 (3 NONE 4 NONE 5 NONE) () () ''
() ())
6 'fem_1d_shapefun' 'lib_fem' 1 ((PROCEDURE UNKNOWN MODULE-PROC DECL
NONE NONE SUBROUTINE) (PROCEDURE 0) 0 0 (7 NONE 8 NONE 9 NONE 10 NONE) ()
() '' () ())
11 'fem_2d_gauss_num' 'lib_fem' 1 ((PROCEDURE UNKNOWN MODULE-PROC DECL
NONE NONE SUBROUTINE) (PROCEDURE 0) 0 0 (12 NONE 13 NONE 14 NONE) () () ''
() ())
15 'fem_2d_gauss_val' 'lib_fem' 1 ((PROCEDURE UNKNOWN MODULE-PROC DECL
NONE NONE SUBROUTINE) (PROCEDURE 0) 0 0 (16 NONE 17 NONE 18 NONE 19 NONE
20 NONE) () () '' () ())
21 'fem_2d_shapefun' 'lib_fem' 1 ((PROCEDURE UNKNOWN MODULE-PROC DECL
NONE NONE SUBROUTINE) (PROCEDURE 0) 0 0 (22 NONE 23 NONE 24 NONE 25 NONE
26 NONE 27 NONE) () () '' () ())
28 'fem_bandmat' 'lib_fem' 1 ((PROCEDURE UNKNOWN MODULE-PROC DECL NONE
NONE DIMENSION FUNCTION) (REAL 8) 0 0 (29 NONE 30 NONE 31 NONE) (2
EXPLICIT (CONSTANT (INTEGER 8) 0 '1') (FUNCTION (INTEGER 8) 0 0 (('' 4 (
VARIABLE (REAL 8) 2 30 ((ARRAY (FULL 0))))) ('dim' 0 (CONSTANT (INTEGER
8) 0 '1')) ('' 0 ())) '_g95_size_8' 0 'size') (CONSTANT (INTEGER 8) 0 '1')
(FUNCTION (INTEGER 8) 0 0 (('' 4 (VARIABLE (REAL 8) 2 31 ((ARRAY (FULL 0)))))
('dim' 0 (CONSTANT (INTEGER 8) 0 '2')) ('' 0 ())) '_g95_size_8' 0 'size'))
() '' () ())
32 'fem_bandvec' 'lib_fem' 1 ((PROCEDURE UNKNOWN MODULE-PROC DECL NONE
NONE DIMENSION FUNCTION INVOKED) (REAL 8) 0 0 (33 NONE 34 NONE 35 NONE)
(1 EXPLICIT (CONSTANT (INTEGER 8) 0 '1') (FUNCTION (INTEGER 8) 0 0 ((''
4 (VARIABLE (REAL 8) 1 35 ((ARRAY (FULL 0))))) ('' 0 ()) ('' 0 ()))
'_g95_size_8' 0 'size')) () '' () ())
36 'fem_choleski_dcmp' 'lib_fem' 1 ((PROCEDURE UNKNOWN MODULE-PROC DECL
NONE NONE SUBROUTINE) (PROCEDURE 0) 0 0 (37 NONE 38 NONE) () () '' () ())
39 'fem_choleski_solv' 'lib_fem' 1 ((PROCEDURE UNKNOWN MODULE-PROC DECL
NONE NONE SUBROUTINE) (PROCEDURE 0) 0 0 (40 NONE 41 NONE) () () '' () ())
42 'fem_elmdof' 'lib_fem' 1 ((PROCEDURE UNKNOWN MODULE-PROC DECL NONE
NONE SUBROUTINE) (PROCEDURE 0) 0 0 (43 NONE 44 NONE 45 NONE) () () '' ()
())
46 'fem_glob2loc_extract' 'lib_fem' 1 ((PROCEDURE UNKNOWN MODULE-PROC
DECL NONE NONE SUBROUTINE) (PROCEDURE 0) 0 0 (47 NONE 48 NONE 49 NONE 50
NONE) () () '' () ())
51 'fem_glob2loc_map' 'lib_fem' 1 ((PROCEDURE UNKNOWN MODULE-PROC DECL
NONE NONE SUBROUTINE) (PROCEDURE 0) 0 0 (52 NONE 53 NONE 54 NONE) () () ''
() ())
55 'fem_m2v' 'lib_fem' 1 ((PROCEDURE UNKNOWN MODULE-PROC DECL NONE NONE
DIMENSION FUNCTION) (REAL 8) 0 0 (56 NONE 57 NONE 58 NONE) (1 EXPLICIT (
CONSTANT (INTEGER 8) 0 '1') (VARIABLE (INTEGER 8) 0 57 ())) () '' () ())
59 'fem_n2d' 'lib_fem' 1 ((PROCEDURE UNKNOWN MODULE-PROC DECL NONE NONE
DIMENSION FUNCTION) (INTEGER 8) 0 0 (60 NONE 61 NONE 62 NONE) (1
EXPLICIT (CONSTANT (INTEGER 8) 0 '1') (OP (INTEGER 8) 0 TIMES (VARIABLE
(INTEGER 8) 0 60 ()) (VARIABLE (INTEGER 8) 0 61 ()))) () '' () ())
63 'fem_v2m' 'lib_fem' 1 ((PROCEDURE UNKNOWN MODULE-PROC DECL NONE NONE
DIMENSION FUNCTION) (REAL 8) 0 0 (64 NONE 65 NONE 66 NONE) (2 EXPLICIT (
CONSTANT (INTEGER 8) 0 '1') (VARIABLE (INTEGER 8) 0 65 ()) (CONSTANT (
INTEGER 8) 0 '1') (VARIABLE (INTEGER 8) 0 66 ())) () '' () ())
67 'lib_fem' 'lib_fem' 1 ((MODULE UNKNOWN UNKNOWN UNKNOWN NONE NONE) (
UNKNOWN) 0 0 () () () '' () ())
66 'ncol' '' 68 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (INTEGER
8) 0 0 () () () '' () ())
65 'nrow' '' 68 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (INTEGER
8) 0 0 () () () '' () ())
64 'vector' '' 68 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION
DUMMY) (REAL 8) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
62 'nodelist' '' 69 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION
DUMMY) (INTEGER 8) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
61 'numdofnode' '' 69 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (
INTEGER 8) 0 0 () () () '' () ())
60 'numnodes' '' 69 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (
INTEGER 8) 0 0 () () () '' () ())
58 'filter' '' 70 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION
OPTIONAL DUMMY) (INTEGER 8) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
57 'n' '' 70 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (INTEGER 8)
0 0 () () () '' () ())
56 'matrix' '' 70 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION
DUMMY) (REAL 8) 0 0 () (2 ASSUMED_SHAPE () () () ()) () '' () ())
54 'g' '' 71 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY) (
INTEGER 8) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
53 'local' '' 71 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DIMENSION
DUMMY) (REAL 8) 0 0 () (2 ASSUMED_SHAPE () () () ()) () '' () ())
52 'global' '' 71 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION
DUMMY) (REAL 8) 0 0 () (2 ASSUMED_SHAPE () () () ()) () '' () ())
50 'numnodeelem' '' 72 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DUMMY) (
INTEGER 8) 0 0 () () () '' () ())
49 'locvector' '' 72 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DIMENSION
DUMMY) (REAL 8) 0 0 () (2 ASSUMED_SHAPE () () () ()) () '' () ())
48 'globvector' '' 72 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION
DUMMY) (REAL 8) 0 0 () (2 ASSUMED_SHAPE () () () ()) () '' () ())
47 'elemnodes' '' 72 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION
DUMMY) (INTEGER 8) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
45 'g' '' 73 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY) (
INTEGER 8) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
44 'nodelist' '' 73 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION
DUMMY) (INTEGER 8) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
43 'numdofnode' '' 73 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (
INTEGER 8) 0 0 () () () '' () ())
41 'b' '' 74 ((VARIABLE INOUT UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY)
(REAL 8) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
40 'a' '' 74 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY) (
REAL 8) 0 0 () (2 ASSUMED_SHAPE () () () ()) () '' () ())
38 'error' '' 75 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DUMMY) (
INTEGER 8) 0 0 () () () '' () ())
37 'a' '' 75 ((VARIABLE INOUT UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY)
(REAL 8) 0 0 () (2 ASSUMED_SHAPE () () () ()) () '' () ())
35 'b' '' 76 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY) (
REAL 8) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
34 'a' '' 76 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY) (
REAL 8) 0 0 () (2 ASSUMED_SHAPE () () () ()) () '' () ())
33 'symmetry' '' 76 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (
CHARACTER 1 ((CONSTANT (INTEGER 8) 0 '1'))) 0 0 () () () '' () ())
31 'b' '' 77 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY) (
REAL 8) 0 0 () (2 ASSUMED_SHAPE () () () ()) () '' () ())
30 'a' '' 77 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY) (
REAL 8) 0 0 () (2 ASSUMED_SHAPE () () () ()) () '' () ())
29 'symmetry' '' 77 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (
CHARACTER 1 ((CONSTANT (INTEGER 8) 0 '1'))) 0 0 () () () '' () ())
27 'der' '' 78 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY)
(REAL 8) 0 0 () (2 ASSUMED_SHAPE () () () ()) () '' () ())
26 'fun' '' 78 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY)
(REAL 8) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
25 'z2' '' 78 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (REAL 8) 0
0 () () () '' () ())
24 'z1' '' 78 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (REAL 8) 0
0 () () () '' () ())
23 'node' '' 78 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY)
(INTEGER 8) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
22 'flag' '' 78 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (INTEGER
8) 0 0 () () () '' () ())
20 'ierr' '' 79 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DUMMY) (
INTEGER 8) 0 0 () () () '' () ())
19 'weight' '' 79 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DIMENSION
DUMMY) (REAL 8) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
18 'coords' '' 79 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DIMENSION
DUMMY) (REAL 8) 0 0 () (2 ASSUMED_SHAPE () () () ()) () '' () ())
17 'numgauss' '' 79 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (
INTEGER 8) 0 0 () () () '' () ())
16 'flag' '' 79 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (INTEGER
8) 0 0 () () () '' () ())
14 'numgauss' '' 80 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DUMMY) (
INTEGER 8) 0 0 () () () '' () ())
13 'numnode_elem' '' 80 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (
INTEGER 8) 0 0 () () () '' () ())
12 'flag' '' 80 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (INTEGER
8) 0 0 () () () '' () ())
10 'der' '' 81 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY)
(REAL 8) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
9 'fun' '' 81 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY)
(REAL 8) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
8 'z' '' 81 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (REAL 8) 0 0
() () () '' () ())
7 'numnodeelem' '' 81 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (
INTEGER 8) 0 0 () () () '' () ())
5 'weight' '' 82 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DIMENSION
DUMMY) (REAL 8) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
4 'coords' '' 82 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DIMENSION
DUMMY) (REAL 8) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
3 'numgauss' '' 82 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (
INTEGER 8) 0 0 () () () '' () ())
)

('fem_1d_gauss_val' 0 2 'fem_1d_shapefun' 0 6 'fem_2d_gauss_num' 0 11
'fem_2d_gauss_val' 0 15 'fem_2d_shapefun' 0 21 'fem_bandmat' 0 28
'fem_bandvec' 0 32 'fem_choleski_dcmp' 0 36 'fem_choleski_solv' 0 39
'fem_elmdof' 0 42 'fem_glob2loc_extract' 0 46 'fem_glob2loc_map' 0 51
'fem_m2v' 0 55 'fem_n2d' 0 59 'fem_v2m' 0 63 'lib_fem' 0 67)
