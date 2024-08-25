%% Define Constant

gamma = 0.5
delta = 0.0001
D = 10

point_set = [-46.68, 16.01;
             -39.23, -15.55;
             84.08, -4.41;
             -77.49, -3.84;
             -121.06, -14.72];

%% Make input
[h, r, epsion] = make_input(point_set, 1);


%% Run Algorithm
iterative_localization(h, r, epsion, gamma, D, delta);