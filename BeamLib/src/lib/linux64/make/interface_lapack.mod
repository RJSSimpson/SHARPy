G95 module created on Tue Apr 19 20:21:29 2011 from ../../src/interface_lapack_v3.f90
If you edit this, you'll get what you deserve.
module-version 8
(() () () () () () () () () () () () () () () () () () () () ())

()

()

()

()

(2 'interface_lapack' 'interface_lapack' 1 ((MODULE UNKNOWN UNKNOWN
UNKNOWN NONE NONE) (UNKNOWN) 0 0 () () () '' () ())
3 'lapack_band_inv' 'interface_lapack' 1 ((PROCEDURE UNKNOWN MODULE-PROC
DECL NONE NONE SUBROUTINE) (PROCEDURE 0) 0 0 (4 NONE 5 NONE 6 NONE 7
NONE 8 NONE 9 NONE) () () '' () ())
10 'lapack_band_luback' 'interface_lapack' 1 ((PROCEDURE UNKNOWN
MODULE-PROC DECL NONE NONE SUBROUTINE INVOKED) (PROCEDURE 0) 0 0 (11
NONE 12 NONE 13 NONE 14 NONE 15 NONE 16 NONE 17 NONE) () () '' () ())
18 'lapack_band_lufact' 'interface_lapack' 1 ((PROCEDURE UNKNOWN
MODULE-PROC DECL NONE NONE SUBROUTINE INVOKED) (PROCEDURE 0) 0 0 (19
NONE 20 NONE 21 NONE 22 NONE 23 NONE 24 NONE) () () '' () ())
25 'lapack_inv' 'interface_lapack' 1 ((PROCEDURE UNKNOWN MODULE-PROC
DECL NONE NONE SUBROUTINE) (PROCEDURE 0) 0 0 (26 NONE 27 NONE 28 NONE 29
NONE) () () '' () ())
30 'lapack_luback' 'interface_lapack' 1 ((PROCEDURE UNKNOWN MODULE-PROC
DECL NONE NONE SUBROUTINE) (PROCEDURE 0) 0 0 (31 NONE 32 NONE 33 NONE 34
NONE 35 NONE) () () '' () ())
36 'lapack_lufact' 'interface_lapack' 1 ((PROCEDURE UNKNOWN MODULE-PROC
DECL NONE NONE SUBROUTINE) (PROCEDURE 0) 0 0 (37 NONE 38 NONE 39 NONE 40
NONE) () () '' () ())
41 'lapack_nonsymeigv' 'interface_lapack' 1 ((PROCEDURE UNKNOWN
MODULE-PROC DECL NONE NONE SUBROUTINE) (PROCEDURE 0) 0 0 (42 NONE 43
NONE 44 NONE 45 NONE 46 NONE 47 NONE 48 NONE 49 NONE 50 NONE) () () '' ()
())
51 'lapack_sparse' 'interface_lapack' 1 ((PROCEDURE UNKNOWN MODULE-PROC
DECL NONE NONE SUBROUTINE) (PROCEDURE 0) 0 0 (52 NONE 53 NONE 54 NONE 55
NONE) () () '' () ())
56 'lapack_sparse_inv' 'interface_lapack' 1 ((PROCEDURE UNKNOWN
MODULE-PROC DECL NONE NONE SUBROUTINE) (PROCEDURE 0) 0 0 (57 NONE 58
NONE 59 NONE 60 NONE 61 NONE) () () '' () ())
61 'ainv' '' 62 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY)
(DERIVED 63) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
63 'sparse' 'lib_sparse' 62 ((DERIVED UNKNOWN UNKNOWN UNKNOWN NONE NONE)
(UNKNOWN) 0 0 () () () '' ((64 'i' (INTEGER 8) () () 0 0 0 ()) (65 'j' (
INTEGER 8) () () 0 0 0 ()) (66 'a' (REAL 8) () () 0 0 0 ())) PUBLIC ())
60 'dimainv' '' 62 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DUMMY) (
INTEGER 8) 0 0 () () () '' () ())
59 'a' '' 62 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY) (
DERIVED 63) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
58 'dima' '' 62 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (INTEGER
8) 0 0 () () () '' () ())
57 'n' '' 62 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (INTEGER 8)
0 0 () () () '' () ())
55 'x' '' 67 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY) (
REAL 8) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
54 'b' '' 67 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY) (
REAL 8) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
53 'sprmat' '' 67 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION
DUMMY) (DERIVED 68) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
68 'sparse' 'lib_sparse' 67 ((DERIVED UNKNOWN UNKNOWN UNKNOWN NONE NONE)
(UNKNOWN) 0 0 () () () '' ((69 'i' (INTEGER 8) () () 0 0 0 ()) (70 'j' (
INTEGER 8) () () 0 0 0 ()) (71 'a' (REAL 8) () () 0 0 0 ())) PUBLIC ())
52 'dimarray' '' 67 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (
INTEGER 8) 0 0 () () () '' () ())
50 'info' '' 72 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DUMMY) (
INTEGER 8) 0 0 () () () '' () ())
49 'vectors' '' 72 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DIMENSION
DUMMY) (COMPLEX 8) 0 0 () (2 ASSUMED_SHAPE () () () ()) () '' () ())
48 'lambda' '' 72 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DIMENSION
DUMMY) (COMPLEX 8) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
47 'numlambda' '' 72 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (
INTEGER 8) 0 0 () () () '' () ())
46 'b' '' 72 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY) (
DERIVED 73) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
73 'sparse' 'lib_sparse' 72 ((DERIVED UNKNOWN UNKNOWN UNKNOWN NONE NONE)
(UNKNOWN) 0 0 () () () '' ((74 'i' (INTEGER 8) () () 0 0 0 ()) (75 'j' (
INTEGER 8) () () 0 0 0 ()) (76 'a' (REAL 8) () () 0 0 0 ())) PUBLIC ())
45 'lenb' '' 72 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (INTEGER
8) 0 0 () () () '' () ())
44 'a' '' 72 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY) (
DERIVED 73) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
43 'lena' '' 72 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (INTEGER
8) 0 0 () () () '' () ())
42 'n' '' 72 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (INTEGER 8)
0 0 () () () '' () ())
40 'info' '' 77 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DUMMY) (
INTEGER 8) 0 0 () () () '' () ())
39 'lupivot' '' 77 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DIMENSION
DUMMY) (INTEGER 8) 0 0 () (1 EXPLICIT (CONSTANT (INTEGER 8) 0 '1') (
VARIABLE (INTEGER 8) 0 37 ())) () '' () ())
38 'a' '' 77 ((VARIABLE INOUT UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY)
(REAL 8) 0 0 () (2 EXPLICIT (CONSTANT (INTEGER 8) 0 '1') (VARIABLE (
INTEGER 8) 0 37 ()) (CONSTANT (INTEGER 8) 0 '1') (VARIABLE (INTEGER 8) 0
37 ())) () '' () ())
37 'n' '' 77 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (INTEGER 8)
0 0 () () () '' () ())
35 'info' '' 78 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DUMMY) (
INTEGER 8) 0 0 () () () '' () ())
34 'b' '' 78 ((VARIABLE INOUT UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY)
(REAL 8) 0 0 () (1 EXPLICIT (CONSTANT (INTEGER 8) 0 '1') (VARIABLE (
INTEGER 8) 0 31 ())) () '' () ())
33 'lupivot' '' 78 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION
DUMMY) (INTEGER 8) 0 0 () (1 EXPLICIT (CONSTANT (INTEGER 8) 0 '1') (
VARIABLE (INTEGER 8) 0 31 ())) () '' () ())
32 'a' '' 78 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY) (
REAL 8) 0 0 () (2 EXPLICIT (CONSTANT (INTEGER 8) 0 '1') (VARIABLE (
INTEGER 8) 0 31 ()) (CONSTANT (INTEGER 8) 0 '1') (VARIABLE (INTEGER 8) 0
31 ())) () '' () ())
31 'n' '' 78 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (INTEGER 8)
0 0 () () () '' () ())
29 'info' '' 79 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DUMMY) (
INTEGER 8) 0 0 () () () '' () ())
28 'inva' '' 79 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY)
(REAL 8) 0 0 () (2 ASSUMED_SHAPE () () () ()) () '' () ())
27 'a' '' 79 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY) (
REAL 8) 0 0 () (2 ASSUMED_SHAPE () () () ()) () '' () ())
26 'n' '' 79 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (INTEGER 8)
0 0 () () () '' () ())
24 'info' '' 80 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DUMMY) (
INTEGER 8) 0 0 () () () '' () ())
23 'lupivot' '' 80 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DIMENSION
DUMMY) (INTEGER 8) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
22 'a' '' 80 ((VARIABLE INOUT UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY)
(REAL 8) 0 0 () (2 ASSUMED_SHAPE () () () ()) () '' () ())
21 'ku' '' 80 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (INTEGER 8)
0 0 () () () '' () ())
20 'kl' '' 80 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (INTEGER 8)
0 0 () () () '' () ())
19 'n' '' 80 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (INTEGER 8)
0 0 () () () '' () ())
17 'info' '' 81 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DUMMY) (
INTEGER 8) 0 0 () () () '' () ())
16 'b' '' 81 ((VARIABLE INOUT UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY)
(REAL 8) 0 0 () (1 EXPLICIT (CONSTANT (INTEGER 8) 0 '1') (VARIABLE (
INTEGER 8) 0 11 ())) () '' () ())
15 'lupivot' '' 81 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION
DUMMY) (INTEGER 8) 0 0 () (1 EXPLICIT (CONSTANT (INTEGER 8) 0 '1') (
VARIABLE (INTEGER 8) 0 11 ())) () '' () ())
14 'a' '' 81 ((VARIABLE INOUT UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY)
(REAL 8) 0 0 () (2 ASSUMED_SHAPE () () () ()) () '' () ())
13 'ku' '' 81 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (INTEGER 8)
0 0 () () () '' () ())
12 'kl' '' 81 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (INTEGER 8)
0 0 () () () '' () ())
11 'n' '' 81 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (INTEGER 8)
0 0 () () () '' () ())
9 'info' '' 82 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DUMMY) (INTEGER
8) 0 0 () () () '' () ())
8 'inva' '' 82 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY)
(REAL 8) 0 0 () (2 ASSUMED_SHAPE () () () ()) () '' () ())
7 'a' '' 82 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY) (
REAL 8) 0 0 () (2 ASSUMED_SHAPE () () () ()) () '' () ())
6 'ku' '' 82 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (INTEGER 8)
0 0 () () () '' () ())
5 'kl' '' 82 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (INTEGER 8)
0 0 () () () '' () ())
4 'n' '' 82 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (INTEGER 8) 0
0 () () () '' () ())
)

('interface_lapack' 0 2 'lapack_band_inv' 0 3 'lapack_band_luback' 0 10
'lapack_band_lufact' 0 18 'lapack_inv' 0 25 'lapack_luback' 0 30
'lapack_lufact' 0 36 'lapack_nonsymeigv' 0 41 'lapack_sparse' 0 51
'lapack_sparse_inv' 0 56)
