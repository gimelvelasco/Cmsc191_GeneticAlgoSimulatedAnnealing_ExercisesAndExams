function ga_sa_10testfunctions_ex6
%%
%VELASCO, Gimel David F.
%2012-58922
%Cmsc 191
%Genetic Algorithm and Simulated Annealing
%Exercise 6
%%
%Test Functions to Solve:
%name:          G I M  E L   V  E L  A S  C O
%equivalent: => 7 9 13 5 12  22 5 12 1 19 3 15
%distinct:   => 1 3 5 7 9 12 13 15 19 22
%mod20:      => 1 2 3 5 7 9 12 13 15 19
%test functions to solve:   => 1 2 3 5 7 9 12 13 15 19
%%
gimelvelasco = [1 1 2 2 3 3 5 5 7 7 9 9 12 13 15 19];
restrictions = [5.12 5.12 5.12 5.12 65.536 65.536 5.12 5.12 600 600 32.768 32.768 12.5 100 3 5.12];
dimensions = [5 10 5 10 5 10 5 10 5 10 5 10 2 2 2 2];
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for num=1:1    %1 to 16 for problems 1 to 10 (contains problems with 5, 10 and 2 dimensions)
for run=1:1
tic;
                       %%%%INITIALIZATION SECTION%%%
test_func = gimelvelasco(num);  %sets the number of w/c test function to be solved
ulb = restrictions(num);        %upper and lower bound
tpl = dimensions(num);        %dimensions
%GENETIC ALGORITHM INITIALIZATION
%STEP 1: Initialization
popsize = tpl*1000;   %Population Size
maxgens = popsize*100;   %Maximum Generations
pc = 0.25;      %crossover rate
pm = 0.1;       %mutation rate
tol = 0.0000025;  %tolerance
Chromosome = zeros(popsize,tpl);
for i=1:popsize
    Chromosome(i,:) = 2*ulb*(rand(1,tpl)-0.5);  %initializing first generation
end
%%
%SIMULATED ANNEALING INITIALIZATION
x_sol = 2*ulb*(rand(1,tpl)-0.5);               %initial guess of the solution x
cooling_ratio = 0.7;                    %sets the cooling ratio to 0.8  i.e. 0.7 < 0.8 < 0.9
num_neigh = 10000;                      %initializes the size of the random neighbors
cooling_sched = zeros(1);               %pre-allocation for speed
cooling_sched(1) = 100;                 %initializes the cooling schedule T0
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%GENETIC ALGORITHM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
for generation=1:maxgens   %LOOPS STEPS 2 to 5
%STEP 2: Selection
%Selection by Tournament
for i=1:popsize/4   %Tournament process
    red = abs(floor(popsize*rand(1))) + 1;  %two random chromosomes to compete
    blu = abs(floor(popsize*rand(1))) + 1;
    if OBJFUNC(Chromosome(red,:),tpl,test_func) > OBJFUNC(Chromosome(blu,:),tpl,test_func)    %competition
        Chromosome(red,:) = Chromosome(blu,:);      %Blue Wins the tournament and overwrites Red
    else
        Chromosome(blu,:) = Chromosome(red,:);      %Red Wins the tournament and overwrites Blue
    end
end
R = rand(1,popsize);
k = zeros(1);
PChromosome = zeros(tpl);
cp = zeros(1);
ctr = 0;    %holds number of parents
for i=1:popsize
    if R(i) < pc
        %select parent
        ctr = ctr + 1;
        k(ctr) = i; %will save the positions of the parent chromosomes
        PChromosome(ctr,:) = Chromosome(i,:);
    end
end
if ctr == 0 %if no parents were selected for the next generation
    continue;
end
%%
%STEP 3: Cross-Over
for i=1:ctr
    cp(i) = floor(abs((tpl-1)*rand(1)))+1;   %crossover points
end
for i=1:ctr-1
   Chromosome(k(i),:) = CROSSOVER(PChromosome(i,:),PChromosome(i+1,:),cp(i),tpl);   %crossover ci and ci+1
end
Chromosome(k(ctr),:) = CROSSOVER(PChromosome(ctr,:),PChromosome(1,:),cp(ctr),tpl);    %crossover ck and c1
%%
%STEP 4: Mutation
%Per Chromosome mutation
mu = round(pm*popsize); %#ofchromosomestomutate = mutationrate*populationsize
for i=1:mu
    cngn = abs(floor(popsize*rand(1))) + 1; %random popsize number
    q = OBJFUNC(Chromosome(cngn,:),tpl,test_func);
    if q < 1
        Chromosome(cngn,:) = 2*ulb*q*(rand(1,tpl)-0.5);   %mutation
    else
        Chromosome(cngn,:) = 2*ulb*(rand(1,tpl)-0.5);
    end
end
%%
%STEP 5: Post-Evaluation
F_obj = zeros(1);
for i=1:popsize
   F_obj(i) = OBJFUNC(Chromosome(i,:),tpl,test_func);
end
%%
fittest = F_obj(1);
for i=1:popsize
    fi = 1;
    if fittest > F_obj(i)
       fittest = F_obj(i);
       fi = i;
    end
end
%fprintf('Fittest: %.16f\t\tRuntime: %.2f seconds\n',fittest,toc);
%disp(Chromosome(fi,:))
if fittest < tol %&& abs(DECODE(Chromosome(fi,:),tpl) - round(DECODE(Chromosome(fi,:),tpl))) < 0.00005
    break;
end
%Step 6: Repeat Generation Iteration
end
%%
%STEP 7: Solution (Best Chromosomes)
if test_func == 1
fprintf('====================DE JONGS FUNCTION=============================\n');
elseif test_func == 2
fprintf('==============AXIS PARALLEL HYPER-ELLIPSOID FUNCTION==============\n');
elseif test_func == 3
fprintf('===============ROTATED HYPER-ELLIPSOID FUNCTION====================\n');
elseif test_func == 5
fprintf('====================RASTRIGINS FUNCTION===========================\n');
elseif test_func == 7
fprintf('====================GRIEWANGKS FUNCTION===========================\n');
elseif test_func == 9
fprintf('=====================ACKLEYS FUNCTION=============================\n');
elseif test_func == 12
fprintf('=====================BRANINS FUNCTION=============================\n');
elseif test_func == 13
fprintf('======================EASOMS FUNCTION=============================\n');
elseif test_func == 15
fprintf('=================SIX-HUMP CAMEL BACK FUNCTION=====================\n');
else
fprintf('====================SHUBERTS FUNCTION=============================\n');
end
fprintf('========================RUN %d====================================\n',run);
fprintf('========================%d DIMENSIONS=============================\n',tpl);
fprintf('====================GENETIC ALGORITHM=============================\n');
fprintf('===================Generations: %d\tPopulation: %d===========\n',maxgens,popsize);
%fprintf('Final Generation\n');
%Chromosome = SORT(Chromosome,popsize,tpl,test_func);
%disp(Chromosome)
fprintf('Generation Number: %d\nPopulation Size: %d\nCrossover Rate: %.2f\nMutation Rate: %.2f\n',maxgens,popsize,pc,pm);
%fprintf('\nThe Fittest Chromosome is\n');
%disp(Chromosome(fi,:)')
%fprintf('Using Tournament Selection\n');
fprintf('Fitness Function Value of %.16f\nTotal Runtime of %f seconds\n',fittest,toc);
fprintf('The Root for Test Function %d is:\n',test_func);
disp(Chromosome(fi,:))
fprintf('==================================================================\n');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%SIMULATED ANNEALING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
sched = 1;                                  %index
while cooling_sched(sched) > 0.000000000000001     %iteration will stop if the cooling temperature reached less than 0.00000001
    T = cooling_sched(sched);               %sets the value of the temperature T
    for j=1:num_neigh
        r  = (cooling_ratio)^sched;             %is used so that the randomness of selecting a neighbor becomes narrower
        x_tmp = 2*ulb*r*(rand(1,tpl)-0.5);    %selects a random neighbor for comparison. with a decreasing amount of randomness
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
    cooling_sched(sched+1) = T*(cooling_ratio)^sched;
    sched = sched+1;
end

%SOLUTION
fprintf('==================SIMULATED ANNEALING=============================\n');
fprintf('With the Objective Function Value of %.16f\nTotal Runtime of %f seconds\nAnd Final Cooling Temperature of %.16f\n',OBJFUNC(x_sol,tpl,test_func),toc,cooling_sched(sched));
fprintf('The Root for Test Function %d is\n',test_func);
disp(x_sol)
fprintf('==================================================================\n');
end
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
        F_obj_ret = abs(F_obj_ret - 0.397887);
    elseif test_func == 13  %Easom
        %Fitness = abs(f(x1,x2) + 1)
        F_obj_ret = (-cos(x(1)))*(cos(x(2)))*exp(-(x(1) - pi)^2 - (x(2) - pi)^2);
        F_obj_ret = abs(F_obj_ret + 1);
    elseif test_func == 15  %Six-hump Camel back
        %Fitness = abs(f(x1,x2) + 1.0316)
        F_obj_ret = (4 - 2.1*x(1)^2 + (x(1)^4)/3)*x(1)^2 + x(1)*x(2) + (-4 + 4*x(2)^2)*x(2)^2;
        F_obj_ret = abs(F_obj_ret + 1.0316);
    else    %Shubert
        F_obj_ret1 = 0;
        F_obj_ret2 = 0;
        for i=1:5
            F_obj_ret1 = F_obj_ret1 + i*cos((i+1)*x(1) + 1);
            F_obj_ret2 = F_obj_ret2 + i*cos((i+1)*x(2) + 1);
        end
        F_obj_ret1 = -F_obj_ret1;
        F_obj_ret = F_obj_ret1*F_obj_ret2;
        F_obj_ret = abs(F_obj_ret + 186.7309);
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