/ linalg.q — Linear Algebra Utilities for qstats
/ qstats: Statistical Testing & Diagnostics for Q/kdb+
/ License: MIT
/ Namespace: .la

/ @desc Extract diagonal of a matrix
/ @param mat:float[][] — input matrix
/ @return float[] — diagonal elements
.la.diag:{[mat]
  n:count mat;
  result:n#0f;
  i:0;
  while[i<n;
    result[i]:mat[i;i];
    i+:1;
  ];
  result
 };

/ @desc Create diagonal matrix from vector
/ @param vec:float[] — diagonal elements
/ @return float[][] — diagonal matrix
.la.diag_matrix:{[vec]
  n:count vec;
  mat:n#enlist n#0f;
  i:0;
  while[i<n;
    mat[i;i]:vec i;
    i+:1;
  ];
  mat
 };

/ @desc Trace of a matrix (sum of diagonal)
/ @param mat:float[][] — input matrix
/ @return float — trace
.la.trace:{[mat]
  sum .la.diag mat
 };

/ @desc Determinant of a matrix
/ @param mat:float[][] — input matrix
/ @return float — determinant
.la.det:{[mat]
  / Use prd over diagonal of upper triangular matrix from QR or LU decomposition
  / For simplicity, use inv and recover det from product formula
  / det(A) = 1/det(inv(A))
  / But this is numerically unstable. Better approach: use eigenvalues
  / For small matrices, use cofactor expansion
  n:count mat;
  /
  / Base cases
  if[n=1;:first first mat];
  if[n=2;:(mat[0;0]*mat[1;1]) - mat[0;1]*mat[1;0]];
  /
  / For larger matrices, use LU-like approach via inv
  / Since q doesn't have built-in LU, use formula: det(A) * det(inv(A)) = 1
  / So det(A) = 1 / det(inv(A))
  / But this is circular. Instead, use eigenvalue product or numerical method.
  /
  / Simple approach: use inv to get inverse, then compute det via trace of log
  / Actually, we can use the identity: det(A) = prd(eigenvalues)
  / But q doesn't have eig. Let's use a numerical trick:
  / det(A) can be approximated by prd diag of upper triangular from QR
  /
  / For now, use recursive cofactor expansion (slow but correct)
  d:0f;
  j:0;
  while[j<n;
    / Cofactor expansion along first row
    / det(A) = sum_j (-1)^j * A[0,j] * det(minor(A,0,j))
    idx:(til n) where not (til n)=j;
    minor:(1_mat)[;idx];
    cofactor:mat[0;j] * .la.det minor;
    d+:$[0=j mod 2;cofactor;neg cofactor];
    j+:1;
  ];
  d
 };

/ @desc Solve linear system Ax = b
/ @param aa:float[][] — coefficient matrix A
/ @param bb:float[] — right-hand side vector b
/ @return float[] — solution vector x
.la.solve:{[aa;bb]
  / Use built-in inv: x = inv(A) mmu b
  (inv aa) mmu bb
 };

/ @desc Power iteration: largest eigenvalue (in absolute value) of A.
/ For symmetric A (in particular A^T A), returns the dominant eigenvalue.
.la.power_iter:{[A;maxiter;tol]
  n:count A;
  v:`float$1+til n;
  v%:sqrt sum v*v;
  lam:0f;
  i:0;
  while[i<maxiter;
    vn:A mmu v;
    lam_new:sqrt sum vn*vn;
    vn%:lam_new;
    if[tol>abs lam_new-lam; i:maxiter];
    lam:lam_new;
    v:vn;
    i+:1];
  lam
 };

/ @desc 2-norm condition number: sigma_max / sigma_min.
/ Computed via power iteration on M = A^T A (whose eigenvalues are the squared
/ singular values of A): sigma_max = sqrt(lambda_max(M)),
/ sigma_min = sqrt(1 / lambda_max(M^-1)).
/ The earlier impl used Frobenius cond (||A||_F * ||A^-1||_F), which is a
/ different metric. 2-norm cond is what numpy/scipy/MATLAB return when called
/ without a norm argument and is the standard convention.
/ @param mat:float[][] — input matrix (square, invertible)
/ @return float — 2-norm condition number
.la.cond:{[mat]
  M:"f"$(flip mat) mmu mat;
  lam_max:.la.power_iter[M;500;1e-14];
  Minv:"f"$inv M;
  lam_max_inv:.la.power_iter[Minv;500;1e-14];
  sqrt lam_max*lam_max_inv
 };

/ @desc Cross-product X'X (efficiently)
/ @param xx:float[][] — input matrix X
/ @return float[][] — X'X
.la.crossprod:{[xx]
  (flip xx) mmu xx
 };

/ @desc Outer product of two vectors
/ @param xx:float[] — first vector
/ @param yy:float[] — second vector
/ @return float[][] — outer product matrix x*y'
.la.outer:{[xx;yy]
  xx*\:yy
 };

/ @desc Identity matrix
/ @param nn:long — matrix dimension
/ @return float[][] — n×n identity matrix
.la.eye:{[nn]
  .la.diag_matrix nn#1f
 };

/ @desc Check if matrix is symmetric
/ @param mat:float[][] — input matrix
/ @param tol:float — tolerance for symmetry check
/ @return boolean — 1b if symmetric, 0b otherwise
.la.is_symmetric:{[mat;tol]
  diff:mat - flip mat;
  all raze (abs diff)<=tol
 };

/ @desc Eigenvalues of a symmetric matrix via Jacobi iteration
/ @param A:float[][] — symmetric matrix (n x n)
/ @return float[] — eigenvalues (unsorted)
/ Uses the classical Jacobi algorithm: iteratively zero out off-diagonal elements
/ via Givens rotations until convergence. Suitable for small-to-medium matrices.
.la.eigen_jacobi:{[A]
  nn:count A;
  / Initialize V = I (eigenvectors, not needed for DW but Jacobi computes both)
  V:(nn;nn)#(nn*nn)?0f; ii:0; while[ii<nn; V[ii;ii]:1f; ii+:1];
  / Jacobi iteration: zero out off-diagonal elements
  / Stopping criterion: max off-diagonal < tol
  tol_conv:1e-12; maxiter:100; iter_ct:0;
  while[iter_ct<maxiter;
    / Find largest off-diagonal element
    max_off:0f; pp:0; qq:1;
    ii:0; while[ii<(nn)-1;
      jj:(ii)+1; while[jj<nn;
        if[max_off<abs A[ii;jj]; max_off:abs A[ii;jj]; pp:ii; qq:jj];
        jj+:1];
      ii+:1];
    / If max off-diagonal is small, converged
    if[max_off<tol_conv; iter_ct:maxiter];
    / Compute rotation angle theta
    if[max_off>=tol_conv;
      theta_ang:$[1e-12>abs (A[pp;pp])-A[qq;qq];
        .special.PI%4f;  / 45 degrees if diagonal elements equal
        0.5*atan (2f*A[pp;qq])%(A[pp;pp])-A[qq;qq]];
      cc:cos theta_ang; sn:sin theta_ang;
      / Apply Givens rotation to A (A' = G^T A G)
      / Update A[p,p], A[q,q], A[p,q], and all A[i,p], A[i,q]
      app:A[pp;pp]; aqq:A[qq;qq]; apq:A[pp;qq];
      A[pp;pp]:((cc*cc*app)+(sn*sn*aqq))-(2f*cc*sn*apq);
      A[qq;qq]:((sn*sn*app)+(cc*cc*aqq))+(2f*cc*sn*apq);
      A[pp;qq]:0f; A[qq;pp]:0f;
      ii:0; while[ii<nn;
        if[(ii<>pp) and ii<>qq;
          aip:A[ii;pp]; aiq:A[ii;qq];
          A[ii;pp]:(cc*aip)-(sn*aiq); A[pp;ii]:A[ii;pp];
          A[ii;qq]:(sn*aip)+(cc*aiq); A[qq;ii]:A[ii;qq]];
        ii+:1];
      iter_ct+:1]];
  / Extract eigenvalues from diagonal of A
  {[A;x] A[x;x]}[A;] each til nn
 };
