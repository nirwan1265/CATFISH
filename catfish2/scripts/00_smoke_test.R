#!/usr/bin/env Rscript
## =============================================================================
## 00_smoke_test.R -- checks the package-backed component layer and times it,
## before any long run is launched. Writes nothing.
##
## NOTE on TFisher: it now contributes Wo (package column tfisher_p_min), a
## STATISTIC, not a p-value. Smaller = more extreme, but it is not uniform
## under the null, so do not read it as a p-value.
## =============================================================================
suppressPackageStartupMessages(library(CATFISH))
source("00_lib.R")
set.seed(1)

G <- 3000L; BLOCK <- 50L; RHO <- 0.3; K <- 150L; POOL <- 600L; M <- 20L
uni   <- make_universe(G, BLOCK, RHO)
paths <- make_pathways(K, M, POOL, G)
pchr  <- paths_to_chr(paths, G)

cat("[1] one null draw through the CATFISH component wrappers\n")
z  <- draw_universe_z(uni)
lp <- stats::pnorm(z, lower.tail = FALSE, log.p = TRUE)
t0 <- Sys.time()
mm <- components_from_draw(lp, z, pchr, G)
el <- as.numeric(diffti      draws 60/200
      draws 80/200
      draws 100/200
      draws 120/200
      draws 140/200
      draws 160/200
      draws 180/200
      draws 200/200
    evaluation draws
      draws 20/200
      draws 40/200
      draws 60/200
      draws 80/200
      draws 100/200
      draws 120/200
      draws 140/200
      draws 160/200
      draws 180/200
      draws 200/200
    drop=<none>    marginal mvn=0.0498 | FWER wy=0.075 bonf=0.050 | thr 4.20e-04 vs 3.33e-04
    drop=ACAT      marginal mvn=0.0499 | FWER wy=0.070 bonf=0.065 | thr 3.50e-04 vs 3.33e-04
    drop=Fisher    marginal mvn=0.0483 | FWER wy=0.045 bonf=0.035 | thr 3.81e-04 vs 3.33e-04
    drop=Stouffer  marginal mvn=0.0486 | FWER wy=0.045 bonf=0.040 | thr 3.44e-04 vs 3.33e-04
    drop=minP      marginal mvn=0.0502 | FWER wy=0.070 bonf=0.065 | thr 3.49e-04 vs 3.33e-04
    drop=TFisher   marginal mvn=0.0498 | FWER wy=0.070 bonf=0.060 | thr 4.10e-04 vs 3.33e-04

  LD=independent m=50 : 12.5 pathways per gene
    calibration draws
      draws 20/200
      draws 40/200
      draws 60/200
      draws 80/200
      draws 100/200
      draws 120/200
      draws 140/200
      draws 160/200
      draws 180/200
      draws 200/200
    evaluation draws
      draws 20/200
      draws 40/200
      draws 60/200
      draws 80/200
      draws 100/200
      draws 120/200
      draws 140/200
      draws 160/200
      draws 180/200
      draws 200/200
    drop=<none>    marginal mvn=0.0443 | FWER wy=0.080 bonf=0.045 | thr 7.57e-04 vs 3.33e-04
    drop=ACAT      marginal mvn=0.0437 | FWER wy=0.080 bonf=0.055 | thr 6.08e-04 vs 3.33e-04
    drop=Fisher    marginal mvn=0.0431 | FWER wy=0.070 bonf=0.045 | thr 8.54e-04 vs 3.33e-04
    drop=Stouffer  marginal mvn=0.0422 | FWER wy=0.090 bonf=0.030 | thr 8.31e-04 vs 3.33e-04
    drop=minP      marginal mvn=0.0439 | FWER wy=0.080 bonf=0.055 | thr 6.07e-04 vs 3.33e-04
    drop=TFisher   marginal mvn=0.0442 | FWER wy=0.085 bonf=0.030 | thr 7.46e-04 vs 3.33e-04

  LD=moderate m=20 : 5.0 pathways per gene
    calibration draws
      draws 20/200
      draws 40/200
      draws 60/200
      draws 80/200
      draws 100/200
      draws 120/200
      draws 140/200
      draws 160/200
      draws 180/200
      draws 200/200
    evaluation draws
      draws 20/200
      draws 40/200
      draws 60/200
      draws 80/200
      draws 100/200
      draws 120/200
      draws 140/200
      draws 160/200
      draws 180/200
      draws 200/200
    drop=<none>    marginal mvn=0.0521 | FWER wy=0.050 bonf=0.115 | thr 6.03e-05 vs 3.33e-04
    drop=ACAT      marginal mvn=0.0521 | FWER wy=0.050 bonf=0.130 | thr 4.83e-05 vs 3.33e-04
    drop=Fisher    marginal mvn=0.0526 | FWER wy=0.045 bonf=0.100 | thr 6.21e-05 vs 3.33e-04
    drop=Stouffer  marginal mvn=0.0520 | FWER wy=0.070 bonf=0.100 | thr 1.30e-04 vs 3.33e-04
    drop=minP      marginal mvn=0.0524 | FWER wy=0.050 bonf=0.130 | thr 4.83e-05 vs 3.33e-04
    drop=TFisher   marginal mvn=0.0519 | FWER wy=0.045 bonf=0.110 | thr 5.18e-05 vs 3.33e-04

  LD=moderate m=50 : 12.5 pathways per gene
    calibration draws
      draws 20/200
      draws 40/200
      draws 60/200
      draws 80/200
      draws 100/200
      draws 120/200
      draws 140/200
      draws 160/200
      draws 180/200
      draws 200/200
    evaluation draws
      draws 20/200
      draws 40/200
      draws 60/200
      draws 80/200
      draws 100/200
      draws 120/200
      draws 140/200
      draws 160/200
      draws 180/200
      draws 200/200
    drop=<none>    marginal mvn=0.0466 | FWER wy=0.070 bonf=0.155 | thr 1.16e-05 vs 3.33e-04
    drop=ACAT      marginal mvn=0.0470 | FWER wy=0.070 bonf=0.165 | thr 9.28e-06 vs 3.33e-04
    drop=Fisher    marginal mvn=0.0458 | FWER wy=0.070 bonf=0.155 | thr 9.63e-06 vs 3.33e-04
    drop=Stouffer  marginal mvn=0.0461 | FWER wy=0.055 bonf=0.125 | thr 2.01e-05 vs 3.33e-04
    drop=minP      marginal mvn=0.0468 | FWER wy=0.070 bonf=0.165 | thr 9.28e-06 vs 3.33e-04
    drop=TFisher   marginal mvn=0.0455 | FWER wy=0.065 bonf=0.155 | thr 9.28e-06 vs 3.33e-04

  LD=strong m=20 : 5.0 pathways per gene
    calibration draws
      draws 20/200
      draws 40/200
      draws 60/200
      draws 80/200
      draws 100/200
      draws 120/200
      draws 140/200
      draws 160/200
      draws 180/200
      draws 200/200
    evaluation draws
      draws 20/200
      draws 40/200
      draws 60/200
      draws 80/200
      draws 100/200
      draws 120/200
      draws 140/200
      draws 160/200
      draws 180/200
      draws 200/200
    drop=<none>    marginal mvn=0.0555 | FWER wy=0.040 bonf=0.180 | thr 6.98e-06 vs 3.33e-04
    drop=ACAT      marginal mvn=0.0555 | FWER wy=0.040 bonf=0.190 | thr 5.59e-06 vs 3.33e-04
    drop=Fisher    marginal mvn=0.0564 | FWER wy=0.045 bonf=0.180 | thr 1.10e-05 vs 3.33e-04
    drop=Stouffer  marginal mvn=0.0575 | FWER wy=0.040 bonf=0.150 | thr 1.08e-05 vs 3.33e-04
    drop=minP      marginal mvn=0.0552 | FWER wy=0.040 bonf=0.190 | thr 5.59e-06 vs 3.33e-04
    drop=TFisher   marginal mvn=0.0538 | FWER wy=0.045 bonf=0.175 | thr 7.21e-06 vs 3.33e-04

  LD=strong m=50 : 12.5 pathways per gene
    calibration draws
      draws 20/200
      draws 40/200
      draws 60/200
      draws 80/200
      draws 100/200
      draws 120/200
      draws 140/200
      draws 160/200
      draws 180/200
      draws 200/200
    evaluation draws
      draws 20/200
      draws 40/200
      draws 60/200
      draws 80/200
      draws 100/200
      draws 120/200
      draws 140/200
      draws 160/200
      draws 180/200
      draws 200/200
    drop=<none>    marginal mvn=0.0156 | FWER wy=0.010 bonf=0.220 | thr 3.13e-08 vs 3.33e-04
    drop=ACAT      marginal mvn=0.0156 | FWER wy=0.010 bonf=0.235 | thr 2.51e-08 vs 3.33e-04
    drop=Fisher    marginal mvn=0.0146 | FWER wy=0.010 bonf=0.205 | thr 2.75e-08 vs 3.33e-04
    drop=Stouffer  marginal mvn=0.0175 | FWER wy=0.010 bonf=0.175 | thr 5.92e-08 vs 3.33e-04
    drop=minP      marginal mvn=0.0156 | FWER wy=0.010 bonf=0.235 | thr 2.51e-08 vs 3.33e-04
    drop=TFisher   marginal mvn=0.0149 | FWER wy=0.005 bonf=0.210 | thr 5.45e-08 vs 3.33e-04

wrote results/03_adaptive_loo.csv
          ld  m  dropped     t1_mvn fwer_wy fwer_bonferroni
 independent 20   <none> 0.04976667   0.075           0.050
 independent 20     ACAT 0.04993333   0.070           0.065
 independent 20   Fisher 0.04833333   0.045           0.035
 independent 20 Stouffer 0.04856667   0.045           0.040
 independent 20     minP 0.05023333   0.070           0.065
 independent 20  TFisher 0.04983333   0.070           0.060
 independent 50   <none> 0.04426667   0.080           0.045
 independent 50     ACAT 0.04370000   0.080           0.055
 independent 50   Fisher 0.04310000   0.070           0.045
 independent 50 Stouffer 0.04216667   0.090           0.030
 independent 50     minP 0.04390000   0.080           0.055
 independent 50  TFisher 0.04420000   0.085           0.030
    moderate 20   <none> 0.05213333   0.050           0.115
    moderate 20     ACAT 0.05213333   0.050           0.130
    moderate 20   Fisher 0.05256667   0.045           0.100
    moderate 20 Stouffer 0.05203333   0.070           0.100
    moderate 20     minP 0.05236667   0.050           0.130
    moderate 20  TFisher 0.05193333   0.045           0.110
    moderate 50   <none> 0.04660000   0.070           0.155
    moderate 50     ACAT 0.04696667   0.070           0.165
    moderate 50   Fisher 0.04583333   0.070           0.155
    moderate 50 Stouffer 0.04606667   0.055           0.125
    moderate 50     minP 0.04680000   0.070           0.165
    moderate 50  TFisher 0.04546667   0.065           0.155
      strong 20   <none> 0.05553333   0.040           0.180
      strong 20     ACAT 0.05550000   0.040           0.190
      strong 20   Fisher 0.05636667   0.045           0.180
      strong 20 Stouffer 0.05746667   0.040           0.150
      strong 20     minP 0.05516667   0.040           0.190
      strong 20  TFisher 0.05376667   0.045           0.175
      strong 50   <none> 0.01560000   0.010           0.220
      strong 50     ACAT 0.01563333   0.010           0.235
      strong 50   Fisher 0.01456667   0.010           0.205
      strong 50 Stouffer 0.01753333   0.010           0.175
      strong 50     minP 0.01563333   0.010           0.235
      strong 50  TFisher 0.01490000   0.005           0.210
[1m[7m%[27m[1m[0m                                                                                                ]7;file://Nirwans-MBP.lan/Users/nirwantandukar/Documents/Github/MAGCAT/catfish2/scripts[0m[27m[24m[J(base) nirwantandukar@Nirwans-MBP scripts % [K[?2004h[7mRscript 00_smoke_test.R[27m[23D[27mR[27ms[27mc[27mr[27mi[27mp[27mt[27m [27m0[27m0[27m_[27ms[27mm[27mo[27mk[27me[27m_[27mt[27me[27ms[27mt[27m.[27mR[?2004l
[1] one null draw through the CATFISH component wrappers
Error: unexpected symbol in "el <- as.numeric(diffti      draws"
Execution halted
[1m[7m%[27m[1m[0m                                                                                                ]7;file://Nirwans-MBP.lan/Users/nirwantandukar/Documents/Github/MAGCAT/catfish2/scripts[0m[27m[24m[J(base) nirwantandukar@Nirwans-MBP scripts % [K[?2004h[7mRscript 00_smoke_test.R[27m[23D[27mR[27ms[27mc[27mr[27mi[27mp[27mt[27m [27m0[27m0[27m_[27ms[27mm[27mo[27mk[27me[27m_[27mt[27me[27ms[27mt[27m.[27mR[?2004l
[1] one null draw through the CATFISH component wrappers
Error: unexpected symbol in "el <- as.numeric(diffti      draws"
Execution halted
[1m[7m%[27m[1m[0m                                                                                                ]7;file://Nirwans-MBP.lan/Users/nirwantandukar/Documents/Github/MAGCAT/catfish2/scripts[0m[27m[24m[J(base) nirwantandukar@Nirwans-MBP scripts % [K[?2004h[7mRscript smoke_check.R[27m[21D[27mR[27ms[27mc[27mr[27mi[27mp[27mt[27m [27ms[27mm[27mo[27mk[27me[27m_[27mc[27mh[27me[27mc[27mk[27m.[27mR[?2004l
[1] one null draw through the CATFISH component wrappers
    150 pathways x 6 stats in 2.22 s
          ACAT  Fisher Stouffer    minP TFisher Omnibus
PW0001 0.64677 0.84368  0.75415 0.77187 0.78656 0.77358
PW0002 0.38049 0.44552  0.32452 0.45141 0.46788 0.41092
PW0003 0.58951 0.05669  0.11460 0.37416 0.11068 0.13472
    non-finite cells: 0

[2] TFisher direction and saturation
    null      min gene p=7.12e-02  TFisher Wo=7.866e-01  Omnibus=7.736e-01
    moderate  min gene p=5.30e-03  TFisher Wo=2.890e-02  Omnibus=3.079e-02
    strong    min gene p=1.39e-08  TFisher Wo=1.000e-50  Omnibus=5.000e-50
    Wo must FALL as signal rises, and must NOT sit at 1e-50

[3] how much of a null bank saturates at the floor
    ACAT      at floor=0.0000   min=   -6.30
    Fisher    at floor=0.0000   min=  -10.78
    Stouffer  at floor=0.0000   min=  -10.95
    minP      at floor=0.0000   min=   -6.29
    TFisher   at floor=0.0000   min=  -11.18
    Omnibus   at floor=0.0000   min=  -10.05

[4] timing
    1.432 s per draw (K=150, m=20)
    B=500    -> 11.9 min per setting
    B=2000   -> 47.7 min per setting
    B=15000  -> 358.0 min per setting
[1m[7m%[27m[1m[0m                                                                                                ]7;file://Nirwans-MBP.lan/Users/nirwantandukar/Documents/Github/MAGCAT/catfish2/scripts[0m[27m[24m[J(base) nirwantandukar@Nirwans-MBP scripts % [K[?2004hppwd[?2004l
/Users/nirwantandukar/Documents/Github/MAGCAT/catfish2/scripts
[1m[7m%[27m[1m[0m                                                                                                ]7;file://Nirwans-MBP.lan/Users/nirwantandukar/Documents/Github/MAGCAT/catfish2/scripts[0m[27m[24m[J(base) nirwantandukar@Nirwans-MBP scripts % [K[?2004h[7mRscript 09_component_calibration.R[27m[34D[27mR[27ms[27mc[27mr[27mi[27mp[27mt[27m [27m0[27m9[27m_[27mc[27mo[27mm[27mp[27mo[27mn[27me[27mn[27mt[27m_[27mc[27ma[27ml[27mi[27mb[27mr[27ma[27mt[27mi[27mo[27mn[27m.[27mR[?2004l
component calibration: K=150 pathways, B_cal=500 B_eval=500

  LD=independent m=20
    calibration draws
      draws 50/500
      draws 100/500
      draws 150/500
      draws 200/500
      draws 250/500
      draws 300/500
      draws 350/500
      draws 400/500
      draws 450/500
      draws 500/500
    evaluation draws
      draws 50/500
      draws 100/500
      draws 150/500
      draws 200/500
      draws 250/500
      draws 300/500
      draws 350/500
      draws 400/500
      draws 450/500
      draws 500/500
    ACAT      analytic=0.0479 (lambda 0.99)  calibrated=0.0461  wy=0.026
    Fisher    analytic=0.0482 (lambda 0.97)  calibrated=0.0468  wy=0.048
    Stouffer  analytic=0.0483 (lambda 0.99)  calibrated=0.0477  wy=0.062
    minP      analytic=0.0477 (lambda 0.99)  calibrated=0.0472  wy=0.026
    TFisher   analytic=0.0714 (lambda 1.19)  calibrated=0.0465  wy=0.030
    Omnibus   analytic=0.0586 (lambda 0.89)  calibrated=0.0467  wy=0.040

  LD=independent m=50
    calibration draws
      draws 50/500
      draws 100/500
      draws 150/500
      draws 200/500
      draws 250/500
      draws 300/500
      draws 350/500
      draws 400/500
      draws 450/500
      draws 500/500
    evaluation draws
      draws 50/500
      draws 100/500
      draws 150/500
      draws 200/500
      draws 250/500
      draws 300/500
      draws 350/500
      draws 400/500
      draws 450/500
      draws 500/500
    ACAT      analytic=0.0497 (lambda 1.00)  calibrated=0.0502  wy=0.048
    Fisher    analytic=0.0505 (lambda 1.05)  calibrated=0.0497  wy=0.024
    Stouffer  analytic=0.0516 (lambda 1.04)  calibrated=0.0537  wy=0.026
    minP      analytic=0.0506 (lambda 1.02)  calibrated=0.0523  wy=0.054
    TFisher   analytic=0.0842 (lambda 1.57)  calibrated=0.0466  wy=0.058
    Omnibus   analytic=0.0621 (lambda 0.98)  calibrated=0.0472  wy=0.044

  LD=moderate m=20
    calibration draws
      draws 50/500
      draws 100/500
      draws 150/500
      draws 200/500
      draws 250/500
      draws 300/500
      draws 350/500
      draws 400/500
      draws 450/500
      draws 500/500
    evaluation draws
      draws 50/500
      draws 100/500
      draws 150/500
      draws 200/500
      draws 250/500
      draws 300/500
      draws 350/500
      draws 400/500
      draws 450/500
      draws 500/500
    ACAT      analytic=0.0514 (lambda 0.94)  calibrated=0.0502  wy=0.046
    Fisher    analytic=0.0720 (lambda 0.85)  calibrated=0.0438  wy=0.056
    Stouffer  analytic=0.0786 (lambda 0.91)  calibrated=0.0415  wy=0.046
    minP      analytic=0.0477 (lambda 0.87)  calibrated=0.0508  wy=0.046
    TFisher   analytic=0.0772 (lambda 1.05)  calibrated=0.0488  wy=0.044
    Omnibus   analytic=0.0771 (lambda 0.78)  calibrated=0.0452  wy=0.044

  LD=moderate m=50
    calibration draws
      draws 50/500
      draws 100/500
      draws 150/500
      draws 200/500
