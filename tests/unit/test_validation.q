/ test_validation.q — Input validation tests
/ Verify that invalid parameters raise 'invalid_arg signals

-1 "\n--- Input Validation Tests ---";

/ Helper to test that a function throws an error with "invalid_arg" in the message
.tst.assert_throws:{[nm;f]
  result:@[{x[]; 0b}; f; {1b}];
  $[result;
    [.tst.pass+:1; -1 "  PASS: ",nm];
    [.tst.fail+:1; .tst.failures,:enlist nm; -1 "  FAIL: ",nm," (expected signal)"]
  ];
 };

-1 "\nTesting normal distribution validation:";
.tst.assert_throws["dnorm with sigma=0"; {.dist.dnorm[0f;0f;0f]}];
.tst.assert_throws["pnorm with sigma=-1"; {.dist.pnorm[0f;0f;-1f]}];
.tst.assert_throws["qnorm with p=0"; {.dist.qnorm[0f;0f;1f]}];
.tst.assert_throws["qnorm with p=1"; {.dist.qnorm[1f;0f;1f]}];
.tst.assert_throws["qnorm with p=-0.1"; {.dist.qnorm[-0.1;0f;1f]}];
.tst.assert_throws["qnorm with p=1.5"; {.dist.qnorm[1.5;0f;1f]}];
.tst.assert_throws["rnorm with n=-1"; {.dist.rnorm[-1;0f;1f]}];
.tst.assert_throws["rnorm with n=0.5 (float)"; {.dist.rnorm[0.5;0f;1f]}];

-1 "\nTesting chi-squared distribution validation:";
.tst.assert_throws["dchisq with df=0"; {.dist.dchisq[1f;0f]}];
.tst.assert_throws["pchisq with df=-1"; {.dist.pchisq[1f;-1f]}];
.tst.assert_throws["qchisq with p=0"; {.dist.qchisq[0f;1f]}];
.tst.assert_throws["qchisq with df=0"; {.dist.qchisq[0.5;0f]}];
.tst.assert_throws["rchisq with n=-5"; {.dist.rchisq[-5;1f]}];

-1 "\nTesting t-distribution validation:";
.tst.assert_throws["dt with df=0"; {.dist.dt[0f;0f]}];
.tst.assert_throws["pt with df=-1"; {.dist.pt[0f;-1f]}];
.tst.assert_throws["qt with p=0"; {.dist.qt[0f;1f]}];
.tst.assert_throws["qt with df=0"; {.dist.qt[0.5;0f]}];
.tst.assert_throws["rt with n=-1"; {.dist.rt[-1;1f]}];

-1 "\nTesting F-distribution validation:";
.tst.assert_throws["df with df1=0"; {.dist.df[1f;0f;1f]}];
.tst.assert_throws["pf with df2=0"; {.dist.pf[1f;1f;0f]}];
.tst.assert_throws["qf with p=0"; {.dist.qf[0f;1f;1f]}];
.tst.assert_throws["qf with df1=0"; {.dist.qf[0.5;0f;1f]}];
.tst.assert_throws["rf with n=-1"; {.dist.rf[-1;1f;1f]}];

-1 "\nTesting uniform distribution validation:";
.tst.assert_throws["dunif with b=a"; {.dist.dunif[0.5;1f;1f]}];
.tst.assert_throws["punif with b<a"; {.dist.punif[0.5;2f;1f]}];
.tst.assert_throws["qunif with p=-0.1"; {.dist.qunif[-0.1;0f;1f]}];
.tst.assert_throws["qunif with p=1.5"; {.dist.qunif[1.5;0f;1f]}];
.tst.assert_throws["runif with n=-1"; {.dist.runif[-1;0f;1f]}];
