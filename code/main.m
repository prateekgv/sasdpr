% Run six demos. Uncomment the run command to exec the demos

% SASD example (simulated data)
fprintf('Running EX1: SASD on simulated dataset\n')
run ./ex1/ex1

% IIR filters as zero-phase matrices
fprintf('Running EX2: Zero-phase filters demo\n')
run ./ex2/demo_filt_matrices

% SASD example (real ECG data)
fprintf('Running EX3: SASD on real ECG dataset\n')
run ./ex3/ex3

% SASDPR example (simulated data)
fprintf('Running EX4: SASDPR on simulated dataset\n')
run ./ex4/ex4

% SASDPR example (sleep spindle)
fprintf('Running EX5: SASDPR for sleep spindle detection\n')
run ./ex5/ex5

% SAPR example (K-complex)
fprintf('Running EX6: SAPR for K-complex detection\n')
run ./ex6/ex6