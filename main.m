%% Define points
point_set = [-46.68, 16.01;
             -39.23, -15.55;
             84.08, -4.41;
             -77.49, -3.84;
             -121.06, -14.72];

%% Make input
[h, r, epsion] = make_input(point_set, 1);


n = height(h) / 2;
m = height(r);

%% Define Optimziation Constant

gamma = 0.5
delta = 0.0001
D = 10

Q1 = (5).^2 * eye(2*n)
Q2 = (0.001).^2 * eye(m)


%% Run Algorithm
iterative_localization(h, r, epsion, gamma, D, delta, Q1, Q2);