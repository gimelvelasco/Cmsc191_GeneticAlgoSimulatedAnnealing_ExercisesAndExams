function geneticalgo_ex4_dejong
%%
%VELASCO, Gimel David F.
%2012-58922
%Cmsc 191
%Genetic Algorithm
%Exercise 4
%%
%De Jong's Function is:
%           f(x) = x(1)^2 + x(2)^2 + ... + x(n)^2;   for n = 1,2,...
%where a,b is in [-5.12,5.12]
%Objective: The Fitness Function must be very close if not equal to 0
%%
clc;    %for ease of documentation
tic;
%STEP 1: Initialization
maxgens = 1000000;   %Maximum Generations
popsize = 100;   %Population Size
pc = 0.25;      %crossover rate
pm = 0.1;       %mutation rate
tpl = 5;        %tpl = n | i.e. n-Tuple De Jong's Function
slctn_proc = 4; %1 - Roulette Wheel | 2 - Ranking | 3 - Elitist | 4 - Tournament
Chromosome = zeros(popsize,tpl);
%%
pnum = zeros(1);
for i=1:popsize
    Chromosome(i,:) = 10.24*(rand(1,tpl)-0.5);  %initializing first generation
    pnum(i) = i;
end
%%
for generation=1:maxgens    %LOOPS STEPS 2 to 5
%%
%STEP 2: Evaluation
F_obj = zeros(1);
for i=1:popsize
   F_obj(i) = OBJFUNC_DEJONG(Chromosome(i,:),tpl);
end
%%
%STEP 3: Selection
if slctn_proc == 1  %Selection by Roulette Wheel
    Total = 0;
    Fitness = zeros(1);
    P = zeros(1);
    C = zeros(1);
    R = zeros(1);
    NewChromosome = zeros(tpl);
    for i=1:popsize
       Fitness(i) = 1/(1+F_obj(i));
       Total = Total + Fitness(i);
    end
    for i=1:popsize
        P(i) = Fitness(i)/Total;
    end
    ctr = 0;
    for i=1:popsize
        ctr = ctr + P(i);
        C(i) = ctr;
    end
    for i=1:popsize
        R(i) = rand(1);
    end
    for i=1:popsize
        NewChromosome(i,:) = zeros(1,tpl);
    for j=1:popsize-1
        if R(i) > C(j) && R(i) <= C(j+1)
            NewChromosome(i,:) = Chromosome(j,:);
        end
    end
    if NewChromosome(i,:) == zeros(1,tpl)
        NewChromosome(i,:) = Chromosome(1,:);
    end
    end
    for i=1:popsize
        Chromosome(i,:) = NewChromosome(i,:);
    end
elseif slctn_proc == 2  %Selection by Ranking
    Chromosome = SORT(Chromosome,popsize,tpl);
    Total = 0;
    Fitness = zeros(1);
    P = zeros(1);
    C = zeros(1);
    R = zeros(1);
    NewChromosome = zeros(tpl);
    for i=1:popsize
       Fitness(i) = 1/(1+F_obj(i));
       Total = Total + Fitness(i);
    end
    for i=1:popsize
        P(i) = Fitness(i)/Total;
    end
    ctr = 0;
    for i=1:popsize
        ctr = ctr + P(i);
        C(i) = ctr;
    end
    for i=1:popsize
        R(i) = rand(1);
    end
    for i=1:popsize
        NewChromosome(i,:) = zeros(1,tpl);
    for j=1:popsize-1
        if R(i) > C(j) && R(i) <= C(j+1)
            NewChromosome(i,:) = Chromosome(j,:);
        end
    end
    if NewChromosome(i,:) == zeros(1,tpl)
        NewChromosome(i,:) = Chromosome(1,:);
    end
    end
    for i=1:popsize
        Chromosome(i,:) = NewChromosome(i,:);
    end
elseif slctn_proc == 3  %Selection by Elitist Strategy
    Chromosome = SORT(Chromosome,popsize,tpl);
    Total = 0;
    Fitness = zeros(1);
    P = zeros(1);
    C = zeros(1);
    R = zeros(1);
    for i=1:popsize
       Fitness(i) = 1/(1+F_obj(i));
       Total = Total + Fitness(i);
    end
    for i=1:popsize
        P(i) = Fitness(i)/Total;
    end
    ctr = 0;
    for i=1:popsize
        ctr = ctr + P(i);
        C(i) = ctr;
    end
    for i=1:popsize
        R(i) = rand(1);
    end
    for i=1:popsize/2
       Chromosome(i + popsize/2,:) = Chromosome(i,:);
    end
elseif slctn_proc == 4  %Selection by Tournament\
    Total = 0;
    Fitness = zeros(1);
    P = zeros(1);
    C = zeros(1);
    R = zeros(1);
    for i=1:popsize
       Fitness(i) = 1/(1+F_obj(i));
       Total = Total + Fitness(i);
    end
    for i=1:popsize
        P(i) = Fitness(i)/Total;
    end
    ctr = 0;
    for i=1:popsize
        ctr = ctr + P(i);
        C(i) = ctr;
    end
    for i=1:popsize
        R(i) = rand(1);
    end
    for i=1:popsize/2   %Tournament process
        red = abs(floor(popsize*rand(1))) + 1;  %two random chromosomes to compete
        blu = abs(floor(popsize*rand(1))) + 1;
        if OBJFUNC_DEJONG(Chromosome(red,:),tpl) > OBJFUNC_DEJONG(Chromosome(blu,:),tpl)    %competition
            Chromosome(red,:) = Chromosome(blu,:);      %Blue Wins the tournament and overwrites Red
        else
            Chromosome(blu,:) = Chromosome(red,:);      %Red Wins the tournament and overwrites Blue
        end
    end
else
    fprintf('No Selection Process Chosen.\n');
    return;
end
%%
%STEP 4: Cross-Over
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
if ctr == 0
    continue;
end
for i=1:ctr
    cp(i) = floor(abs(tpl*rand(1)))+1;   %crossover points
end
for i=1:ctr-1
   Chromosome(k(i),:) = CROSSOVER(PChromosome(i,:),PChromosome(i+1,:),cp(i),tpl);   %crossover ci and ci+1
end
Chromosome(k(ctr),:) = CROSSOVER(PChromosome(ctr,:),PChromosome(1,:),cp(ctr),tpl);    %crossover ck and c1
%%

%STEP 5: Mutation
mu = round(pm*popsize); %#ofmutations = mutationrate*tuples*populationsize
for i=1:mu
    cngn = abs(floor(popsize*rand(1))) + 1; %random popsize number
    Chromosome(cngn,:) = 10.24*(rand(1,tpl)-0.5);   %mutation
end

totalgen = tpl*popsize;
mu = round(pm*totalgen); %#ofmutations = mutationrate*tuples*populationsize
for i=1:mu
    cngn = round(totalgen*rand(1)); %random totalgen number
    cn = floor(cngn/tpl) + 1;          %chromosome number
    gn = mod(cngn,tpl) + 1;            %generation number
    rmn = 10.24*(rand(1)-0.5);     %random mutation number
    Chromosome(cn,gn) = rmn;     %actual mutation event
end
%%
%STEP 6: Solution (Best Chromosomes)
%%
%Final Evaluation
for i=1:popsize
   F_obj(i) = OBJFUNC_DEJONG(Chromosome(i,:),tpl);
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
%Step 7: Repeat Generation Iteration
fprintf('Fittest: %f\t\tRuntime: %f seconds\n',fittest,toc);
if fittest < 0.00005
    break;
end
%pause
end
%%
fprintf('===================Generations: %d\tPopulation: %d===================\n',maxgens,popsize);
fprintf('=======================%d-Tuple De Jongs Function======================\n',tpl);
fprintf('Final Generation\n');
Chromosome = SORT(Chromosome,popsize,tpl);
disp(Chromosome)
fprintf('*Note: Dump Chromosome stays in the final row.\n\n');
fprintf('Generation Number: %d\nPopulation Size: %d\nCrossover Rate: %.2f\nMutation Rate: %.2f\n',maxgens,popsize,pc,pm);
fprintf('\nThe Fittest Chromosome is\n');
disp(Chromosome(fi,:))
if slctn_proc == 1
    fprintf('Using Roulette Wheel Selection\n');
elseif slctn_proc == 2
    fprintf('Using Ranking Selection\n');
elseif slctn_proc == 3
    fprintf('Using Elitist Selection\n');
elseif slctn_proc == 4
    fprintf('Using Tournament Selection\n');
end
fprintf('With the Fitness Value of %f\nand Total Runtime of %f seconds\n',fittest,toc);
fprintf('=====================================================================\n');
%%
end

function F_obj_ret = OBJFUNC_DEJONG(x,tpl)
    F_obj_ret = 0;
    for i=1:tpl
       F_obj_ret = F_obj_ret + x(i)^2;
    end
    %F_obj_ret = abs(F_obj_ret);
    %Objective: The Fitness Function must be very close if not equal to 0
end

function Chromosome_ret = CROSSOVER(PChromosome1,PChromosome2,cp_rec,tpl)
    Chromosome_ret = PChromosome1;
    for i=cp_rec+1:tpl    %crossovers after the crossover point
        Chromosome_ret(i) = PChromosome2(i);
    end
end

function Chromosome_ret = SORT(Chromosome,popsize_rec,tpl)
    %bubble sort
    for i=1:popsize_rec-1
       for j=1:popsize_rec-1
          if OBJFUNC_DEJONG(Chromosome(j,:),tpl) > OBJFUNC_DEJONG(Chromosome(j+1,:),tpl)
              %swap
              temp = Chromosome(j,:);
              Chromosome(j,:) = Chromosome(j+1,:);
              Chromosome(j+1,:) = temp;
          end
       end
    end
    Chromosome_ret = Chromosome;
end