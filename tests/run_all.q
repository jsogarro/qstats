/ run_all.q — Master test runner for qstats
/ Usage: q tests/run_all.q
/ Exit code: 0 = all pass, 1 = any fail

/ Load library
\l src/load.q

/ Test framework
.tst.pass:0;
.tst.fail:0;
.tst.failures:();

.tst.assert_approx:{[nm;actual;expected;tol]
  / Guard against null/NaN: q treats 0n as smaller than any float, so without
  / this check `tol>abs(actual-expected)` is vacuously true and NaN passes.
  has_null:any null actual;
  ok:(not has_null) and all (tol>abs actual-expected) | (actual=expected);
  $[ok;
    [.tst.pass+:1; -1 "  PASS: ",nm];
    [.tst.fail+:1; .tst.failures,:enlist nm; -1 "  FAIL: ",nm," actual=",(string actual)," expected=",string expected]
  ];
 };

.tst.assert_equal:{[nm;actual;expected]
  ok:actual~expected;
  $[ok;
    [.tst.pass+:1; -1 "  PASS: ",nm];
    [.tst.fail+:1; .tst.failures,:enlist nm; -1 "  FAIL: ",nm," actual=",(-3!actual)," expected=",-3!expected]
  ];
 };

.tst.assert:{[nm;cond]
  $[cond;
    [.tst.pass+:1; -1 "  PASS: ",nm];
    [.tst.fail+:1; .tst.failures,:enlist nm; -1 "  FAIL: ",nm]
  ];
 };

/ Load and run tie-out tests (reference JSON must exist)
-1 "\n=== qstats Test Suite ===\n";

/ Check if reference files exist
ref_exists:0<count key `:tests/reference/special.json;
if[not ref_exists; -1 "WARNING: Reference files not found. Run: cd tests/reference && python gen_all.py"; -1 "Skipping tie-out tests.\n"];

/ Load tie-out tests (system commands must be at top level, not in if)
\l tests/tieout/test_special.q
\l tests/tieout/test_distributions.q
\l tests/unit/test_linalg.q
\l tests/unit/test_validation.q
\l tests/unit/test_random_variates.q

/ Summary
-1 "\n============================================================";
-1 "  RESULTS: ",string[.tst.pass]," passed, ",string[.tst.fail]," failed";
if[0<.tst.fail; -1 "  FAILURES: ",", " sv .tst.failures];
-1 "============================================================";

exit $[0=.tst.fail;0;1];
