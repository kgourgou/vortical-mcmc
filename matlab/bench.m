% Monte - Carlo with vortices

% proposal matrix
Q = [[0,2/3.0,1/3.0];
     [1/3.0,0,2/3.0];
     [2/3.0,1/3.0,0]];

% invariant distribution of Q
q = [1/3.0,1/3.0,1/3.0];


% Target invariant distribution 

p(1) = 0.15;
p(2) = 0.33;
p(3) = 1-p(1)-p(2);

N = 10000; % samples to compute

x0 = 1; % Initial point with which to start

% Easy pick of the constant (see Bierken's paper)
a = min(p/q);

x = vortMCMC(Q,q,p,x0,a,N);


% Frequencies 

fprintf('Frequency of picking state 1 : %f\n',mean(x==1));
fprintf('Frequency of picking state 2 : %f\n',mean(x==2));