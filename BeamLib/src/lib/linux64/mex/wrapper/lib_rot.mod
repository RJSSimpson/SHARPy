G95 module created on Wed Aug 31 18:57:50 2011 from ../../src/lib_rot.f90
If you edit this, you'll get what you deserve.
module-version 8
(() () () () () () () () () () () () () () () () () () () () ())

()

()

()

()

(2 'lib_rot' 'lib_rot' 1 ((MODULE UNKNOWN UNKNOWN UNKNOWN NONE NONE) (
UNKNOWN) 0 0 () () () '' () ())
3 'reshape' '(intrinsic)' 1 ((PROCEDURE UNKNOWN INTRINSIC UNKNOWN NONE
NONE FUNCTION) (UNKNOWN) 0 0 () () () '' () ())
4 'rot_beta2mat' 'lib_rot' 1 ((PROCEDURE UNKNOWN MODULE-PROC DECL NONE
NONE SUBROUTINE) (PROCEDURE 0) 0 0 (5 NONE 6 NONE) () () '' () ())
7 'rot_cross' 'lib_rot' 1 ((PROCEDURE UNKNOWN MODULE-PROC DECL NONE NONE
DIMENSION FUNCTION) (REAL 8) 0 0 (8 NONE 9 NONE) (1 EXPLICIT (CONSTANT (
INTEGER 8) 0 '1') (CONSTANT (INTEGER 8) 0 '3')) () '' () ())
10 'rot_epsilon' 'lib_rot' 1 ((PARAMETER UNKNOWN UNKNOWN UNKNOWN NONE
NONE) (REAL 8) 0 0 () (CONSTANT (REAL 8) 0 '4611686018427388' 0 1013) ()
() '' () ())
11 'rot_epsilonmin' 'lib_rot' 1 ((PARAMETER UNKNOWN UNKNOWN UNKNOWN NONE
NONE) (REAL 8) 0 0 () (CONSTANT (REAL 8) 0 '4835703278458517' 0 993) ()
() '' () ())
12 'rot_outprod' 'lib_rot' 1 ((PROCEDURE UNKNOWN MODULE-PROC DECL NONE
NONE DIMENSION FUNCTION) (REAL 8) 0 0 (13 NONE 14 NONE) (2 EXPLICIT (
CONSTANT (INTEGER 8) 0 '1') (CONSTANT (INTEGER 8) 0 '3') (CONSTANT (
INTEGER 8) 0 '1') (CONSTANT (INTEGER 8) 0 '3')) () '' () ())
15 'rot_phi2mat' 'lib_rot' 1 ((PROCEDURE UNKNOWN MODULE-PROC DECL NONE
NONE DIMENSION FUNCTION) (REAL 8) 0 0 (16 NONE) (2 EXPLICIT (CONSTANT (
INTEGER 8) 0 '1') (CONSTANT (INTEGER 8) 0 '3') (CONSTANT (INTEGER 8) 0 '1')
(CONSTANT (INTEGER 8) 0 '3')) () '' () ())
17 'rot_points2mat' 'lib_rot' 1 ((PROCEDURE UNKNOWN MODULE-PROC DECL
NONE NONE SUBROUTINE) (PROCEDURE 0) 0 0 (18 NONE 19 NONE 20 NONE) () () ''
() ())
21 'rot_skew' 'lib_rot' 1 ((PROCEDURE UNKNOWN MODULE-PROC DECL NONE NONE
DIMENSION FUNCTION INVOKED) (REAL 8) 0 0 (22 NONE) (2 EXPLICIT (
CONSTANT (INTEGER 8) 0 '1') (CONSTANT (INTEGER 8) 0 '3') (CONSTANT (
INTEGER 8) 0 '1') (CONSTANT (INTEGER 8) 0 '3')) () '' () ())
23 'rot_vect' 'lib_rot' 1 ((PROCEDURE UNKNOWN MODULE-PROC DECL NONE NONE
DIMENSION FUNCTION) (REAL 8) 0 0 (24 NONE) (1 EXPLICIT (CONSTANT (
INTEGER 8) 0 '1') (CONSTANT (INTEGER 8) 0 '3')) () '' () ())
24 'skewmatrix' '' 25 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION
DUMMY) (REAL 8) 0 0 () (2 EXPLICIT (CONSTANT (INTEGER 8) 0 '1') (
CONSTANT (INTEGER 8) 0 '3') (CONSTANT (INTEGER 8) 0 '1') (CONSTANT (
INTEGER 8) 0 '3')) () '' () ())
22 'vector' '' 26 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION
DUMMY) (REAL 8) 0 0 () (1 EXPLICIT (CONSTANT (INTEGER 8) 0 '1') (
CONSTANT (INTEGER 8) 0 '3')) () '' () ())
20 'rotmatrix' '' 27 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DIMENSION
DUMMY) (REAL 8) 0 0 () (2 EXPLICIT (CONSTANT (INTEGER 8) 0 '1') (
CONSTANT (INTEGER 8) 0 '3') (CONSTANT (INTEGER 8) 0 '1') (CONSTANT (
INTEGER 8) 0 '3')) () '' () ())
19 'point2' '' 27 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION
DUMMY) (REAL 8) 0 0 () (1 EXPLICIT (CONSTANT (INTEGER 8) 0 '1') (
CONSTANT (INTEGER 8) 0 '3')) () '' () ())
18 'point1' '' 27 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION
DUMMY) (REAL 8) 0 0 () (1 EXPLICIT (CONSTANT (INTEGER 8) 0 '1') (
CONSTANT (INTEGER 8) 0 '3')) () '' () ())
16 'phi' '' 28 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY)
(REAL 8) 0 0 () (1 EXPLICIT (CONSTANT (INTEGER 8) 0 '1') (CONSTANT (
INTEGER 8) 0 '3')) () '' () ())
14 'b' '' 29 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY) (
REAL 8) 0 0 () (1 EXPLICIT (CONSTANT (INTEGER 8) 0 '1') (CONSTANT (
INTEGER 8) 0 '3')) () '' () ())
13 'a' '' 29 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY) (
REAL 8) 0 0 () (1 EXPLICIT (CONSTANT (INTEGER 8) 0 '1') (CONSTANT (
INTEGER 8) 0 '3')) () '' () ())
9 'vec2' '' 30 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY)
(REAL 8) 0 0 () (1 EXPLICIT (CONSTANT (INTEGER 8) 0 '1') (CONSTANT (
INTEGER 8) 0 '3')) () '' () ())
8 'vec1' '' 30 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY)
(REAL 8) 0 0 () (1 EXPLICIT (CONSTANT (INTEGER 8) 0 '1') (CONSTANT (
INTEGER 8) 0 '3')) () '' () ())
6 'rotmatrix' '' 31 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DIMENSION
DUMMY) (REAL 8) 0 0 () (2 ASSUMED_SHAPE () () () ()) () '' () ())
5 'beta' '' 31 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY)
(REAL 8) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
)

('lib_rot' 0 2 'reshape' 0 3 'rot_beta2mat' 0 4 'rot_cross' 0 7
'rot_epsilon' 0 10 'rot_epsilonmin' 0 11 'rot_outprod' 0 12 'rot_phi2mat'
0 15 'rot_points2mat' 0 17 'rot_skew' 0 21 'rot_vect' 0 23)
