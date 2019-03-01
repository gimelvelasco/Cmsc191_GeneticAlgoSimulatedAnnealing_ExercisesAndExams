%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VELASCO, Gimel David F.
%2012-58922
%Cmsc 191
%Genetic Algorithm and Simulated Annealing
%Long Exam 1
%Part 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%Initialization Section%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%GA
popsize = 1000;                         %Population Size
pc = 0.9;                               %crossover rate
pm = 0.01;                               %mutation rate
maxgens = 1000;
tol = 0.000025;
%%SA
initial_T = 1000;                       %initial temperature of the system
cooling_stop = 0.000001;      %cooling stops when it reaches this temperature
cooling_ratio = 0.8;                    %sets the cooling ratio to 0.8  i.e. 0.7 < 0.8 < 0.9
num_neigh = 100;                 %number of neighbors to consider
ca = 0.99;                              %cooling constant a in Geometric cooling
cb = 0.99;                              %cooling constant b in Another cooling
ci = 3;                                 %cooling implementation:  | 1 - Geometric Cooling | 2 - Another Cooling w/ small b | ow - Exponential Cooling |
%%HSA
popsize2 = 1000;                         %Population Size
pc2 = 0.9;                               %crossover rate
pm2 = 0.01;                               %mutation rate
initial_T2 = 1000;                       %initial temperature of the system
cooling_stop2 = 0.01;      %cooling stops when it reaches this temperature
cooling_ratio2 = 0.9;                    %sets the cooling ratio to 0.8  i.e. 0.7 < 0.8 < 0.9
num_neigh2 = popsize*pm;                 %number of neighbors to consider
ca2 = 0.8;                              %cooling constant a in Geometric cooling
cb2 = 0.99;                              %cooling constant b in Another cooling
ci2 = 3;                                 %cooling implementation:  | 1 - Geometric Cooling | 2 - Another Cooling w/ small b | ow - Exponential Cooling |
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_func = 1;                          % | '0' - 21-tuple test function | 'ow' - 3-tuple quadratic equation |
tpl = 3;        %number of coefficients
ulb = 100;        %upper and lower bound
maxruns = 5;                            %number of runs
%%
y_sim = zeros(1);
x = [-1 -2 0 3 4 5];
y = [2 3 3 18 27 38];
xmin = -1;
xmax = 5;
%%
a = zeros(tpl);
e = zeros(1);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for runs=1:maxruns
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%CALL ALGORITHM HERE%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %only run one at a time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[a(runs,:), e(runs)] = runSA(x,y,xmax,xmin,cooling_ratio,tpl,test_func,initial_T,cooling_stop,ulb,num_neigh,runs,ci,ca,cb);
%[a(runs,:), e(runs)] = runHSA(x,y,xmax,xmin,ulb,tpl,popsize2,pc2,cooling_ratio2,num_neigh2,initial_T2,cooling_stop2,test_func,runs,ci2,ca2,cb2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SOLUTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
fprintf('====================LONG EXAM 1 PART 2=============================\n');
fprintf('Coefficients from %d runs:\n',maxruns)
disp(a)

a_aveg = zeros(tpl);
e_aveg = 0;
for i=1:tpl
    for runs=1:maxruns
        a_aveg(1,i) = a_aveg(1,i) + a(runs,i);
        e_aveg = e_aveg + e(runs);
    end
end
for i=1:tpl
    a_aveg(1,i) = a_aveg(1,i)/maxruns;
end
e_aveg = e_aveg/maxruns;

fprintf('Averaged Resulting Coefficients of the function:\n');
disp(a_aveg(1,:))
fprintf('Averaged Energy Value: \t\t %.16f\n',e_aveg);