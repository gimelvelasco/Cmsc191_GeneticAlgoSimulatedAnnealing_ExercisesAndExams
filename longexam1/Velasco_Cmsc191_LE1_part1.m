function Velasco_Cmsc191_LE1_part1
%%
%VELASCO, Gimel David F.
%2012-58922
%Cmsc 191
%Genetic Algorithm and Simulated Annealing
%Long Exam 1
%Part 1
%%
%Test Functions to Solve:
%name:          G I M  E L   V  E L  A S  C O
%equivalent: => 7 9 13 5 12  22 5 12 1 19 3 15
%distinct:   => 1 3 5 7 9 12 13 15 19 22
%mod20:      => 1 2 3 5 7 9 12 13 15 19
%test functions to solve:   => 1 2 3 5 7 9 12 13 15 19
%%
gimelvelasco = [1 2 3 5 7 9 12 13 15 19];
restrictions = [5.12 5.12 65.536 5.12 600 32.768 12.5 100 3 5.12];
dimensions = [2 2 2 2 2 2 2 2 2 2];
%%
popsize_tests = [100 1000];
pc_tests = [1.0 0.9 0.8];
pm_tests = [0.05 0.01];
cooling_tests = [0.7 0.9];
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%Initialization Section%%%%%%%%%%%%%%%%%%%%%%%%%%%
initial_T = 1000;                       %initial temperature of the system
cooling_stop = 0.0000000000000001;      %cooling stops when it reaches this temperature
maxruns = 1;                            %number of runs
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
fprintf('==================HYBRID SIMULATED ANNEALING=======================\n');
for num=1:1   %1 to 10 for problems 1 to 10 (all have 2 dimensions)
    test_func = gimelvelasco(num);  %sets the number of w/c test function to be solved
    if test_func == 1
    fprintf('====================DE JONGS FUNCTION==============================\n');
    elseif test_func == 2
    fprintf('==============AXIS PARALLEL HYPER-ELLIPSOID FUNCTION===============\n');
    elseif test_func == 3
    fprintf('===============ROTATED HYPER-ELLIPSOID FUNCTION====================\n');
    elseif test_func == 5
    fprintf('====================RASTRIGINS FUNCTION============================\n');
    elseif test_func == 7
    fprintf('====================GRIEWANGKS FUNCTION============================\n');
    elseif test_func == 9
    fprintf('=====================ACKLEYS FUNCTION==============================\n');
    elseif test_func == 12
    fprintf('=====================BRANINS FUNCTION==============================\n');
    elseif test_func == 13
    fprintf('======================EASOMS FUNCTION==============================\n');
    elseif test_func == 15
    fprintf('=================SIX-HUMP CAMEL BACK FUNCTION======================\n');
    else
    fprintf('====================SHUBERTS FUNCTION==============================\n');
    end
for popsize_index=2:2
    popsize = popsize_tests(popsize_index);                         %Population Size
for pc_index=2:2
    pc = pc_tests(pc_index);                               %crossover rate
for pm_index=1:1
    pm = pm_tests(pm_index);                               %mutation rate
for cooling_index=1:1
    cooling_ratio = cooling_tests(cooling_index);                    %sets the cooling ratio to 0.8  i.e. 0.7 < 0.8 < 0.9
    %fprintf('Population Size: %d\tCrossover Rate: %.2f\tMutation Rate: %.2f\tCooling Ratio: %.1f\n',popsize,pc,pm,cooling_ratio);
    ulb = restrictions(num);        %upper and lower bound
    tpl = dimensions(num);        %dimensions
    num_neigh = popsize*pm;                 %number of neighbors to consider
    energy_array = zeros(1);
    runtime_array = zeros(1);
for runs=1:maxruns
tic;

%%%%%%%%%%%%%%%%%%%%%%%%%HYBRID SIMULATED ANNEALING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%I
cooling_sched = zeros(1);               %pre-allocation for speed
cooling_sched(1) = initial_T;                 %initializes the cooling schedule T0
%II
Chromosome = zeros(popsize,tpl);
for i=1:popsize
    Chromosome(i,:) = 2*ulb*(rand(1,tpl)-0.5);  %initializing first generation
end
%%
sched = 1;                              %index / iteration
while cooling_sched(sched) > cooling_stop     %iteration will stop if the cooling temperature reached less than 0.00000001
    T = cooling_sched(sched);               %sets the value of the temperature T
    %III.a. Do N/2 times
    for j=1:popsize/2
        %III.a.i. Select two parents at random
        red = abs(floor(popsize*rand(1))) + 1;  %two random chromosomes to compete
        blu = abs(floor(popsize*rand(1))) + 1;  %red and blu hold the index of the parents
        %III.a.ii. Generate two offsprings
        %%Recombination Operator (CROSSOVER)
        pc_trial = rand(1);
        if pc_trial < pc     %if trial made it in the crossover rate
            cp = floor(abs((tpl-1)*rand(1)))+1;   %random crossover point
            Child_Chromosome(1,:) = CROSSOVER(Chromosome(red,:),Chromosome(blu,:),cp,tpl);     %crossover red and blu
            Child_Chromosome(2,:) = CROSSOVER(Chromosome(blu,:),Chromosome(red,:),cp,tpl);     %they will have two children
            %%Neighborhood Operator (MUTATION)
            for k=1:2
                x_sol = Child_Chromosome(k,:);               
                for i=1:num_neigh
                    adrs = abs(floor(popsize*rand(1))) + 1; %gets a random address of a neighbor within the population
                    x_tmp = Chromosome(adrs,:);    %selects a random neighbor for comparison. with a decreasing amount of randomness
                    if OBJFUNC(x_tmp,tpl,test_func) < OBJFUNC(x_sol,tpl,test_func)  %if the neighbor is better, change the solution
                        x_sol = x_tmp;
                    elseif OBJFUNC(x_tmp,tpl,test_func) > OBJFUNC(x_sol,tpl,test_func)  %if not, change the solution if it is lucky
                        delta = OBJFUNC(x_tmp,tpl,test_func) - OBJFUNC(x_sol,tpl,test_func);
                        p = P(delta,T);
                        q = rand(1);
                        if q <= p
                            x_sol = x_tmp; 
                        end
                    end
                end
                Child_Chromosome(k,:) = x_sol;           %will overwrite the child based on the neighborhood operator
            end
            %%III.a.iii. Boltzman Trials
            ARpossibility = rand(1);                    % <0.5 - Single Acceptance/Rejection | >=0.5 - Double Acceptance/Rejection
            if ARpossibility < 0.5 %%Case 1: Double Acceptance/Rejection
                E1 = OBJFUNC(Chromosome(red,:),tpl,test_func) + OBJFUNC(Chromosome(blu,:),tpl,test_func);
                E2 = OBJFUNC(Child_Chromosome(1,:),tpl,test_func) + OBJFUNC(Child_Chromosome(2,:),tpl,test_func);
                bp = BOLTZMAN(E1,E2,T);
                bp_trial = rand(1);
                if bp_trial >= bp
                    %%III.a.iv. Overwrite Parents with the Trial Winner
                    Chromosome(red,:) = Child_Chromosome(1,:);
                    Chromosome(blu,:) = Child_Chromosome(2,:);
                end
            else %%Case 2: Single Acceptance/Rejection
                E1 = OBJFUNC(Chromosome(red,:),tpl,test_func);
                E2 = OBJFUNC(Child_Chromosome(1,:),tpl,test_func);
                bp = BOLTZMAN(E1,E2,T);
                bp_trial = rand(1);
                if bp_trial >= bp
                    %%III.a.iv. Overwrite Parents with the Trial Winner
                    Chromosome(red,:) = Child_Chromosome(1,:);  %offsprings wins the trial
                end

                E1 = OBJFUNC(Chromosome(red,:),tpl,test_func);
                E2 = OBJFUNC(Child_Chromosome(2,:),tpl,test_func);
                bp = BOLTZMAN(E1,E2,T);
                bp_trial = rand(1);
                if bp_trial >= bp
                    %%III.a.iv. Overwrite Parents with the Trial Winner
                    Chromosome(red,:) = Child_Chromosome(2,:);  %offsprings wins the trial
                end

                E1 = OBJFUNC(Chromosome(blu,:),tpl,test_func);
                E2 = OBJFUNC(Child_Chromosome(1,:),tpl,test_func);
                bp = BOLTZMAN(E1,E2,T);
                bp_trial = rand(1);
                if bp_trial >= bp
                    %%III.a.iv. Overwrite Parents with the Trial Winner
                    Chromosome(blu,:) = Child_Chromosome(1,:);  %offsprings wins the trial
                end

                E1 = OBJFUNC(Chromosome(blu,:),tpl,test_func);
                E2 = OBJFUNC(Child_Chromosome(2,:),tpl,test_func);
                bp = BOLTZMAN(E1,E2,T);
                bp_trial = rand(1);
                if bp_trial >= bp
                    %%III.a.iv. Overwrite Parents with the Trial Winner
                    Chromosome(blu,:) = Child_Chromosome(2,:);  %offsprings wins the trial
                end
            end
            %%
        else                    %if the whole trial did not make it inside the crossover rate, it will have a tournament
            if OBJFUNC(Chromosome(red,:),tpl,test_func) > OBJFUNC(Chromosome(blu,:),tpl,test_func)    %competition
                Chromosome(red,:) = Chromosome(blu,:);      %Blue Wins the tournament and overwrites Red
            else
                Chromosome(blu,:) = Chromosome(red,:);      %Red Wins the tournament and overwrites Blue
            end
        end
    end
    %III.b. Periodically Lower T
    cooling_sched(sched+1) = T*(cooling_ratio)^sched;
    sched = sched+1;
end
%%
%Post-Evaluation
F_obj = zeros(1);
for i=1:popsize
   F_obj(i) = OBJFUNC(Chromosome(i,:),tpl,test_func);
end
%%
fittest = F_obj(1);
for i=1:popsize
    %fi = 1;
    if fittest > F_obj(i)
       fittest = F_obj(i);
       %fi = i;
    end
end
%%
%SOLUTION
%uncomment this section to see results
%%
%fprintf('===============================RUN %d===============================\n',runs);
%%
%fprintf('Minimum Energy:\t\t\t\t\t\t\t\t\t%.16f\nTotal Runtime:\t\t\t\t\t\t\t\t\t%f seconds\nFinal Temperature:\t\t\t\t\t\t\t\t%.16f\n',fittest,toc,cooling_sched(sched));
%fprintf('Global Minimum is at:\t\t\t\t\t\t');
%disp(Chromosome(fi,:))
%fprintf('===================================================================\n');
%%
energy_array(runs) = fittest;
runtime_array(runs) = toc;
%%
end
energy_aveg = 0;
runtime_aveg = 0;
for a=1:maxruns
    energy_aveg = energy_aveg + energy_array(a);
    runtime_aveg = runtime_aveg + runtime_array(a);
end
energy_aveg = energy_aveg/maxruns;
runtime_aveg = runtime_aveg/maxruns;
fprintf('Population Size: %d\tCrossover Rate: %.2f\tMutation Rate: %.2f\tCooling Ratio: %.1f\n',popsize,pc,pm,cooling_ratio);
fprintf('Energy Average: %.16f\t\tRuntime Average: %.16f\n\n',energy_aveg,runtime_aveg);
end
end
end
end
fprintf('\n\n\n');
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F_obj_ret = OBJFUNC(x,tpl,test_func)   %Objective: The Fitness Function must be very close if not equal to 0
    if test_func == 1       %DeJong
        F_obj_ret = 0;
        for i=1:tpl
           F_obj_ret = F_obj_ret + x(i)^2;
        end
        F_obj_ret = abs(F_obj_ret);
    elseif test_func == 2   %Axis Parallel Hyper-ellipsoid
        F_obj_ret = 0;
        for i=1:tpl
           F_obj_ret = F_obj_ret + i*x(i)^2;
        end
        F_obj_ret = abs(F_obj_ret);
    elseif test_func == 3   %Rotated Hyper-ellipsoid
        F_obj_ret = 0;
        for i=1:tpl
            for j=1:i
                F_obj_ret = F_obj_ret + x(j)^2;
            end
        end
        F_obj_ret = abs(F_obj_ret);
    elseif test_func == 5   %Rastrigin
        F_obj_ret = 0;
        for i=1:tpl
           F_obj_ret = F_obj_ret + x(i)^2 - 10*cos(2*pi*x(i));
        end
        F_obj_ret = F_obj_ret + 10*tpl;
        F_obj_ret = abs(F_obj_ret);
    elseif test_func == 7   %Griewangk
        F_obj_reta = 0;
        F_obj_retb = 0;
        for i=1:tpl
           F_obj_reta = F_obj_reta + x(i)^2;
           F_obj_retb = F_obj_retb*(cos(x(i)/sqrt(i)));
        end
        F_obj_reta = (1/4000)*F_obj_reta;
        F_obj_retb = F_obj_retb + 1;
        F_obj_ret = F_obj_reta - F_obj_retb;
        F_obj_ret = abs(F_obj_ret);
    elseif test_func == 9   %Ackley
        a = 20;
        b = 0.2;
        c = 2*pi;
        F_obj_ret1 = 0;
        F_obj_ret2 = 0;
        for i=1:tpl
            F_obj_ret1 = F_obj_ret1 + x(i)^2;
            F_obj_ret2 = F_obj_ret2 + (cos(c*x(i)));
        end
        F_obj_ret1 = (-a)*exp((-b)*sqrt((1/tpl)*F_obj_ret1));
        F_obj_ret2 = exp((1/tpl)*F_obj_ret2);
        F_obj_ret1 = F_obj_ret1 - F_obj_ret2;
        F_obj_ret2 = a + exp(1);
        F_obj_ret = F_obj_ret1 + F_obj_ret2;
        F_obj_ret = abs(F_obj_ret);
    elseif test_func == 12  %Branin
        %Fitness = abs(f(x1,x2) - 0.397887)
        a = 1;
        b = 5.1/(4*(pi^2));
        c = 5/pi;
        d = 6;
        e = 10;
        f = 1/(8*pi);
        F_obj_ret = a*(x(2) - b*x(1)^2 + c*x(1) - d)^2 + e*(1-f)*cos(x(1)) + e;
        F_obj_ret = F_obj_ret - 0.397887;
        F_obj_ret = abs(F_obj_ret);
    elseif test_func == 13  %Easom
        %Fitness = abs(f(x1,x2) + 1)
        expr1 = -cos(x(1))*(cos(x(2)));
        expra = x(1) - pi;
        expra = expra^2;
        expra = -expra;
        exprb = x(2) - pi;
        exprb = exprb^2;
        exprb = -exprb;
        expr2 = exp(expra + exprb);
        F_obj_ret = (expr1*expr2);
        F_obj_ret = F_obj_ret + 1;
        F_obj_ret = abs(F_obj_ret);
    elseif test_func == 15  %Six-hump Camel back
        %Fitness = abs(f(x1,x2) + 1.0316)
        F_obj_ret = (4 - 2.1*x(1)^2 + (x(1)^4)/3)*x(1)^2 + x(1)*x(2) + (-4 + 4*x(2)^2)*x(2)^2;
        F_obj_ret = F_obj_ret + 1.0316;
        F_obj_ret = abs(F_obj_ret);
    else    %Shubert
        F_obj_ret1 = 0;
        F_obj_ret2 = 0;
        for i=1:5
            F_obj_ret1 = F_obj_ret1 + i*cos((i+1)*x(1) + 1);
            F_obj_ret2 = F_obj_ret2 + i*cos((i+1)*x(2) + 1);
        end
        F_obj_ret1 = -F_obj_ret1;
        F_obj_ret = F_obj_ret1*F_obj_ret2;
        F_obj_ret = F_obj_ret + 186.7309;
        %F_obj_ret = abs(F_obj_ret);
    end
end

function Chromosome_ret = CROSSOVER(PChromosome1,PChromosome2,cp_rec,tpl)
    Chromosome_ret = PChromosome1;
    for i=cp_rec+1:tpl    %crossovers after the crossover point
        Chromosome_ret(i) = PChromosome2(i);
    end
end

function p_ret = P(delta,T)
    p_ret = exp(-(delta)/T);
end

function bp_ret = BOLTZMAN(Ei,Ej,T)
    bp_ret = 1/(1 + exp((Ei-Ej)/T));
end