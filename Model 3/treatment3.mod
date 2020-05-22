param num_beams; #number of avaiable beamlets
param num_rows >= 1;
param num_cols >= 1;

set ROWS:= 1..num_rows;
set COLS:= 1..num_cols;
set BEAM:= 1..num_beams;

param beam_value {BEAM, ROWS, COLS} >= 0;

param tumor_value {ROWS, COLS} >= 0;

param critical_value {ROWS, COLS} >= 0;

param critical_max;
param tumor_min;

param a {j in ROWS, k in COLS} = if critical_value[j,k] > 0 then 1 else 0;
param e {j in ROWS, k in COLS} = if tumor_value[j,k] > 0 then 1 else 0;

var X {i in BEAM} >= 0;

# slack variables
var P {j in ROWS, k in COLS} >= 0;
var Q {j in ROWS, k in COLS} >= 0;

# minimize total slacks
minimize total_slack: sum {j in ROWS, k in COLS} ( (P[j,k] + Q[j,k]) + sum {i in BEAM} sum {m in -1 .. 1} sum {n in -1 .. 1} (1 - a[j,k]) * (a[min(max(j+m,1), num_rows),min(max(k+n,1),num_cols)]) * X[i] * beam_value[i,j,k]);

minimize treatment_dosage: (sum{i in BEAM} sum {j in ROWS, k in COLS} a[j,k] * X[i] * beam_value[i,j,k]);

subject to tumor_limit {j in ROWS, k in COLS} : sum {i in BEAM} X[i] * beam_value[i,j,k] * e[j,k] >= (tumor_min - Q[j,k]) * e[j,k] ;

subject to critical_limit {j in ROWS, k in COLS} : sum {i in BEAM} a[j,k] * X[i] * beam_value[i,j,k] <= critical_max + P[j,k];





