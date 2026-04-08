/ test_special.q — Tie-out tests for special functions
/ Compare against scipy reference values

-1 "\n--- Special Functions Tie-Out Tests ---";

/ Load reference data
ref_json:read0 `:tests/reference/special.json;
ref:.j.k raze ref_json;

/ Test lgamma
-1 "\nTesting .special.lgamma:";
lgamma_data:ref`lgamma;
lgamma_x:lgamma_data[`inputs][`x];
lgamma_expected:lgamma_data[`outputs][`y];
lgamma_actual:.special.lgamma lgamma_x;
{.tst.assert_approx["lgamma(x=",string[x],")";y;z;1e-10]}'[lgamma_x;lgamma_actual;lgamma_expected];

/ Test betainc
-1 "\nTesting .special.betainc:";
betainc_cases:ref`betainc;
{[case]
  aa:case`a;
  bb:case`b;
  x_vals:case`x;
  y_expected:case`y;
  {[aa;bb;xv;ye]
    act:.special.betainc[xv;aa;bb];
    nm:"betainc(x=",string[xv],", a=",string[aa],", b=",string[bb],")";
    .tst.assert_approx[nm;act;ye;1e-8];
  }[aa;bb;;] '[x_vals;y_expected];
 } each betainc_cases;

/ Test gammainc
-1 "\nTesting .special.gammainc:";
gammainc_cases:ref`gammainc;
{[case]
  aa:case`a;
  x_vals:case`x;
  y_expected:case`y;
  {[aa;xv;ye]
    act:.special.gammainc[aa;xv];
    nm:"gammainc(a=",string[aa],", x=",string[xv],")";
    .tst.assert_approx[nm;act;ye;1e-8];
  }[aa;;] '[x_vals;y_expected];
 } each gammainc_cases;
