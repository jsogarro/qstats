/ test_linalg.q — Unit tests for linear algebra utilities
/ Compare against numpy reference values

-1 "\n--- Linear Algebra Unit Tests ---";

/ Load reference data
ref:.j.k raze read0 `:tests/reference/linalg.json;

/ Test matrix 4x4
-1 "\nTesting 4x4 matrix operations:";
mat4x4_data:ref`matrix_4x4;
mat4x4:mat4x4_data`A;
diag_expected:mat4x4_data`diag;
trace_expected:mat4x4_data`trace;
det_expected:mat4x4_data`det;
cond_expected:mat4x4_data`cond;

/ Test diag extraction
diag_actual:.la.diag mat4x4;
.tst.assert_approx["diag(A)";diag_actual;diag_expected;1e-10];

/ Test trace
trace_actual:.la.trace mat4x4;
.tst.assert_approx["trace(A)";trace_actual;trace_expected;1e-10];

/ Test determinant
det_actual:.la.det mat4x4;
.tst.assert_approx["det(A)";det_actual;det_expected;1e-6];

/ Test condition number
cond_actual:.la.cond mat4x4;
.tst.assert_approx["cond(A)";cond_actual;cond_expected;1e-6];

/ Test solve 3x3
-1 "\nTesting solve 3x3:";
solve_data:ref`solve_3x3;
solve_a:solve_data`A;
solve_b:solve_data`b;
solve_x_expected:solve_data`x;
solve_x_actual:.la.solve[solve_a;solve_b];
.tst.assert_approx["solve(A,b)";solve_x_actual;solve_x_expected;1e-10];

/ Test crossprod
-1 "\nTesting crossprod:";
test_x:(2 3f;4 5f;6 7f);
crossprod_actual:.la.crossprod test_x;
crossprod_expected:(flip test_x) mmu test_x;
.tst.assert_equal["crossprod(X)";crossprod_actual;crossprod_expected];

/ Test eye
-1 "\nTesting eye:";
eye3:.la.eye 3;
eye3_expected:(1 0 0f;0 1 0f;0 0 1f);
.tst.assert_equal["eye(3)";eye3;eye3_expected];

/ Test is_symmetric
-1 "\nTesting is_symmetric:";
sym_mat:(1 2 3f;2 4 5f;3 5 6f);
nonsym_mat:(1 2 3f;4 5 6f;7 8 9f);
.tst.assert["is_symmetric(symmetric)";.la.is_symmetric[sym_mat;1e-10]];
.tst.assert["is_symmetric(non-symmetric)";not .la.is_symmetric[nonsym_mat;1e-10]];

/ Test outer product
-1 "\nTesting outer:";
outer_x:1 2 3f;
outer_y:4 5f;
outer_actual:.la.outer[outer_x;outer_y];
outer_expected:(4 5f;8 10f;12 15f);
.tst.assert_equal["outer([1,2,3],[4,5])";outer_actual;outer_expected];

/ Test diag_matrix
-1 "\nTesting diag_matrix:";
diag_vec:1 2 3f;
diag_mat_actual:.la.diag_matrix diag_vec;
diag_mat_expected:(1 0 0f;0 2 0f;0 0 3f);
.tst.assert_equal["diag_matrix([1,2,3])";diag_mat_actual;diag_mat_expected];
