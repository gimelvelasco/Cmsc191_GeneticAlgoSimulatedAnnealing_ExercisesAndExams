function geneticalgo_ex4
%%
%VELASCO, Gimel David F.
%2012-58922
%Cmsc 191
%Genetic Algorithm
%Exercise 4
%%
%Find the roots with 4-digit accuracy of the following functions:
%1.
%           f(x) = x - cos(x)           on [-5,5]
%2.
%           f(x) = exp(-x)*(x-2)        on [-5,5]
%3.
%           f(x) = x^2 - x - 12         on [-5,5]
%%
%The Fitness function of each is
%
%               abs(f(x))
%
%Objective: The Fitness Function must be very close if not equal to 0 (Minimization)
%%
%Each Chromosome is represented by 5 values of which their total is the value of the root
%Per allele is a real number on [-1,1]
%%
clc;    %for ease of documentation
tic;
%%%%INITIALIZATION%%%
%STEP 1: Initialization
maxgens = 100000;   %Maximum Generations
popsize = 100;   %Population Size
pc = 0.25;      %crossover rate
pm = 0.1;       %mutation rate
ulb = 5;        %upper and lower bound
tpl = 5;        %default is 5. 5 alleles per chromosome
tol = 0.00005;  %tolerance
prob_num = 3;   %1 - Problem #1 | 2 - Problem #2 | o.w. - Problem #3 
Chromosome = zeros(popsize,tpl);
%%
for i=1:popsize
    Chromosome(i,:) = 2*(ulb/tpl)*(rand(1,tpl)-0.5);  %initializing first generation
end
%%
for generation=1:maxgens    %LOOPS STEPS 2 to 5
%%
%STEP 2: Evaluation
F_obj = zeros(1);
for i=1:popsize
   F_obj(i) = OBJFUNC(Chromosome(i,:),tpl,prob_num);
end
%%
%STEP 3: Selection
%Selection by Tournament
for i=1:popsize/4   %Tournament process
    red = abs(floor(popsize*rand(1))) + 1;  %two random chromosomes to compete
    blu = abs(floor(popsize*rand(1))) + 1;
    if OBJFUNC(Chromosome(red,:),tpl,prob_num) > OBJFUNC(Chromosome(blu,:),tpl,prob_num)    %competition
        Chromosome(red,:) = Chromosome(blu,:);      %Blue Wins the tournament and overwrites Red
    else
        Chromosome(blu,:) = Chromosome(red,:);      %Red Wins the tournament and overwrites Blue
    end
end
%computing probabilities

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
%STEP 4: Cross-Over
for i=1:ctr
    cp(i) = floor(abs((tpl-1)*rand(1)))+1;   %crossover points
end
for i=1:ctr-1
   Chromosome(k(i),:) = CROSSOVER(PChromosome(i,:),PChromosome(i+1,:),cp(i),tpl);   %crossover ci and ci+1
end
Chromosome(k(ctr),:) = CROSSOVER(PChromosome(ctr,:),PChromosome(1,:),cp(ctr),tpl);    %crossover ck and c1
%%
%STEP 5: Mutation
%Per Chromosome mutation
mu = round(pm*popsize); %#ofchromosomestomutate = mutationrate*populationsize
for i=1:mu
    cngn = abs(floor(popsize*rand(1))) + 1; %random popsize number
    Chromosome(cngn,:) = 2*(ulb/tpl)*(rand(1,tpl)-0.5);   %mutation
end
%%
%Per Allele mutation
%totalgen = tpl*popsize;
%mu = round(pm*totalgen); %#ofmutations = mutationrate*tuples*populationsize
%for i=1:mu
%    cngn = round(totalgen*rand(1));    %random totalgen number
%    cn = floor(cngn/tpl) + 1;          %chromosome number
%    gn = mod(cngn,tpl) + 1;            %generation number
%    rmn = 2*(ulb/tpl)*(rand(1)-0.5);                 %random mutation number
%    Chromosome(cn,gn) = rmn;           %actual mutation event
%end
%%
%STEP 6: Post-Evaluation
for i=1:popsize
   F_obj(i) = OBJFUNC(Chromosome(i,:),tpl,prob_num);
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
fprintf('Fittest: %.16f\t\tRuntime: %.2f seconds\n',fittest,toc);
if fittest < tol %&& abs(DECODE(Chromosome(fi,:),tpl) - round(DECODE(Chromosome(fi,:),tpl))) < 0.00005
    break;
end
%Step 7: Repeat Generation Iteration
end
%%
%STEP 8: Solution (Best Chromosomes)
clc;
fprintf('===================Generations: %d\tPopulation: %d===================\n',maxgens,popsize);
fprintf('Final Generation\n');
Chromosome = SORT(Chromosome,popsize,tpl,prob_num);
disp(Chromosome)
fprintf('Generation Number: %d\nPopulation Size: %d\nCrossover Rate: %.2f\nMutation Rate: %.2f\n',maxgens,popsize,pc,pm);
fprintf('\nThe Fittest Chromosome is\n');
disp(Chromosome(fi,:)')
fprintf('Using Tournament Selection\n');
fprintf('With the Fitness Value of %.16f\nand Total Runtime of %f seconds\n',fittest,toc);
fprintf('The Root for Problem #%d is x = %.16f\n',prob_num,DECODE(Chromosome(fi,:),tpl));
fprintf('=====================================================================\n');
%%
end

function F_obj_ret = OBJFUNC(x,tpl,prob_num)
    t = DECODE(x,tpl);
    if prob_num == 1
        F_obj_ret = t - cos(t);
    elseif prob_num == 2
        F_obj_ret = exp(-t)*(t-2);
    else
        F_obj_ret = t^2 - t - 12;
    end
    F_obj_ret = abs(F_obj_ret);
end

function Root = DECODE(x,tpl)
    Root = 0;
    for i=1:tpl
        Root = Root + x(i);
    end
end

function Chromosome_ret = CROSSOVER(PChromosome1,PChromosome2,cp_rec,tpl)
    Chromosome_ret = PChromosome1;
    for i=cp_rec+1:tpl    %crossovers after the crossover point
        Chromosome_ret(i) = PChromosome2(i);
    end
end

function Chromosome_ret = SORT(Chromosome,popsize_rec,tpl,prob_num)
    %bubble sort
    for i=1:popsize_rec-1
       for j=1:popsize_rec-1
          if OBJFUNC(Chromosome(j,:),tpl,prob_num) > OBJFUNC(Chromosome(j+1,:),tpl,prob_num)
              %swap
              temp = Chromosome(j,:);
              Chromosome(j,:) = Chromosome(j+1,:);
              Chromosome(j+1,:) = temp;
          end
       end
    end
    Chromosome_ret = Chromosome;
end