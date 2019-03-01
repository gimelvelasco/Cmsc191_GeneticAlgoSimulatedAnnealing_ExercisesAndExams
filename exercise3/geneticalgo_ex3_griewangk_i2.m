function geneticalgo_ex3_griewangk_i2
%%
%VELASCO, Gimel David F.
%2012-58922
%Cmsc 191
%Genetic Algorithm
%Exercise 3 (1 of 4)
%%
%GRIEWANGK'S FUNCTION
%The 2-tuple Griewangk's Function is:
%           f(a,b) = 1 + ((a^2)+(b^2))/4000 - (cos(a/sqrt(1))*cos(b/sqrt(2)))
%where a,b is in [-600,600]
%Objective: The Fitness Function must be very close if not equal to 0
%%
clc;    %for ease of documentation
mg_array = [10 50 500];
ps_array = [10 50 500];
pc_array = [0.75 0.5 0.25];
for mg_index=1:3
for ps_index=1:3
tic;
fittest = 9999;    %dump value
while fittest > 0.5
%STEP 1: Initialization
maxgens = mg_array(mg_index);   %Maximum Generations %{10, 50, 500}
popsize = ps_array(ps_index);   %Population Size %{10, 50, 500}
pc = pc_array(ps_index);       %crossover rate %{0.75, 0.5, 0.25}
pm = 0.1;       %mutation rate
tpl = 2;        %i = 2 i.e. 2-Tuple Griewangk Function &and for ease of debugging
%initialize/preallocate zero 2d chromosome array of size popsize x tuple here
Chromosome = zeros(popsize+1,tpl);
%%
pnum = zeros(1);
for i=1:popsize
    Chromosome(i,:) = round(1200*(rand(1,tpl)-0.5));  %initializing first generation
    pnum(i) = i;
end
%%
for generation=1:maxgens    %LOOPS STEPS 2 to 5
%%
%STEP 2: Evaluation
F_obj = zeros(1);
for i=1:popsize
   F_obj(i) = OBJFUNC_GRIEWANGK(Chromosome(i,:));
end
%%
%STEP 3: Selection
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
    fprintf('No Parents for next generation was selected. Crossover Rate might be too small.\n');
    continue;
end
for i=1:ctr
    cp(i) = 1;   %crossover points %for i=2, the crossover point will only be at the center of the chromosome
end
for i=1:ctr-1
   Chromosome(k(i),:) = CROSSOVER(PChromosome(i,:),PChromosome(i+1,:),cp(i));   %crossover ci and ci+1
end
Chromosome(k(ctr),:) = CROSSOVER(PChromosome(ctr,:),PChromosome(1,:),cp(ctr));    %crossover ck and c1
%%

%STEP 5: Mutation
totalgen = tpl*popsize;
mu = round(pm*totalgen); %#ofmutations = mutationrate*tuples*populationsize
for i=1:mu
    cngn = round(totalgen*rand(1)); %random totalgen number
    cn = floor(cngn/tpl) + 1;          %chromosome number
    gn = mod(cngn,tpl) + 1;            %generation number
    rmn = round(1200*(rand(1)-0.5));     %random mutation number
    Chromosome(cn,gn) = rmn;     %actual mutation event
end
%%
%STEP 6: Repeat Steps 2 to 5 until maximum generations reached
end
%%
%STEP 7: Solution (Best Chromosomes)
%%
%Final Evaluation
for i=1:popsize
   F_obj(i) = OBJFUNC_GRIEWANGK(Chromosome(i,:));
end
%%
fittest = F_obj(1);
for i=1:popsize
    if fittest > F_obj(i)
       fittest = F_obj(i);
       cn = i;
    end
end
%%
%Step 8: Repeat Runs until the desired fitness value is met
end %while loop
pnum(popsize+1) = popsize+1;
F_obj(popsize+1) = 9999;
fprintf('===================Generations: %d\tPopulation: %d===================\n',maxgens,popsize);
%fprintf('Final Generation\n');
%Tf = table(pnum',Chromosome(:,1),Chromosome(:,2),F_obj(:),'VariableNames',{'Number','a','b','Fitness'});
%disp(Tf)
fprintf('Note: Dump Chromosome stays in the final row.\n');
fprintf('Generation Number: %d\nPopulation Size: %d\nCrossover Rate: %f\nMutation Rate: %f\n',maxgens,popsize,pc,pm);
fprintf('\nThe Fittest Chromosome is\n');
Tsol = table(Chromosome(cn,1),Chromosome(cn,2),fittest,'VariableNames',{'a','b','Fitness'});
disp(Tsol)
fprintf('With the Fitness Value of %f\nand Total Runtime of %f seconds\n',fittest,toc);
fprintf('=====================================================================\n');
%%
end
end
end

function F_obj_ret = OBJFUNC_GRIEWANGK(x)
    F_obj_ret = abs(1 + ((x(1)^2)+(x(2)^2))/4000 - (cos(x(1)/sqrt(1))*cos(x(2)/sqrt(2))));
    %Objective: The Fitness Function must be very close if not equal to 0
end

function Chromosome_ret = CROSSOVER(PChromosome1,PChromosome2,cp_rec)
    Chromosome_ret = PChromosome1;
    for i=cp_rec+1:2    %crossovers after the crossover point
        Chromosome_ret(i) = PChromosome2(i);
    end
end