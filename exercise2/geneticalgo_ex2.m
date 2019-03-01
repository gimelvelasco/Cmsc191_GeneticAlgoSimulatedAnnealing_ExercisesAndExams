function geneticalgo_ex2
%%
%VELASCO, Gimel David F.
%2012-58922
%Cmsc 191
%Genetic Algorithm
%Exercise 2
%%
%The Fitness Function is defined as:
%           f(x) = a + 2b + 3c + 4d + 5e + 6f - 240
%such that a,b,c,d,e,f is in [0,45]
%Note: The Fitness Function must be very close if not equal to 0
%%
%for g=1:50
fittest = 9999;    %dump value
while fittest ~= 0
clc;    %for ease of documentation
%STEP 1: Initialization
maxgens = 100;   %Maximum Generations
popsize = 50;   %Population Size %(50,10),(50,50),(100,10),(100,50)
pc = 0.25;       %crossover rate
pm = 0.1;       %mutation rate
%%
for i=1:popsize
    Chromosome(i,:) = round(45*rand(1,6));  %initializing first generation
    pnum(i) = i;
end
%%
for generation=1:maxgens    %LOOPS STEPS 2 to 5
%%
%STEP 2: Evaluation
for i=1:popsize
   F_obj(i) = OBJFUNC(Chromosome(i,:));
end
%%
%STEP 3: Selection
Total = 0;
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
    NewChromosome(i,:) = zeros(1,6);
for j=1:popsize-1
    if R(i) > C(j) && R(i) <= C(j+1)
        NewChromosome(i,:) = Chromosome(j,:);
    end
end
if NewChromosome(i,:) == zeros(1,6)
    NewChromosome(i,:) = Chromosome(1,:);
end
end
for i=1:popsize
    Chromosome(i,:) = NewChromosome(i,:);
end
%%
%STEP 4: Cross-Over
R(:) = rand(1,popsize);
ctr = 0;
for i=1:popsize
    if R(i) < pc
        %select parent
        ctr = ctr + 1;
        k(ctr) = i; %will save the positions of the parent chromosomes
        PChromosome(ctr,:) = Chromosome(i,:);
    end
end
if ctr == 0
    fprintf('No Parents for next generation was selected. Crossover Rate might be too small.\n End of Program.\n');
    return;
end
for i=1:ctr
    cp(i) = 1 + round(3*rand(1));   %crossover points
end
for i=1:ctr-1
   Chromosome(k(i),:) = CROSSOVER(PChromosome(i,:),PChromosome(i+1,:),cp(i));   %crossover ci and ci+1
end
Chromosome(k(ctr),:) = CROSSOVER(PChromosome(ctr,:),PChromosome(1,:),cp(ctr));    %crossover ck and c1
%%

%STEP 5: Mutation
totalgen = 6*popsize;
mu = round(pm*totalgen); %#ofmutations = mutationrate*tuples*populationsize
for i=1:mu
    cngn = round(totalgen*rand(1)); %random totalgen number
    cn = floor(cngn/6) + 1;          %chromosome number
    gn = mod(cngn,6) + 1;            %generation number
    rmn = round(45*rand(1));     %random mutation number
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
   F_obj(i) = OBJFUNC(Chromosome(i,:));
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
end
pnum(popsize+1) = popsize+1;
fprintf('Final Generation\n');
Tf = table(pnum',Chromosome(:,1),Chromosome(:,2),Chromosome(:,3),Chromosome(:,4),Chromosome(:,5),Chromosome(:,6),'VariableNames',{'Number','a','b','c','d','e','f'})
fprintf('The Fittest Chromosome is\n');
Tsol = table(Chromosome(cn,1),Chromosome(cn,2),Chromosome(cn,3),Chromosome(cn,4),Chromosome(cn,5),Chromosome(cn,6),'VariableNames',{'a','b','c','d','e','f'})
fprintf('With the Fitness Value of %d\n',fittest);
%Soln(g,:) = Chromosome(cn,:);
%%
%end
%fprintf('Solutions Found for the 6-tuple function\n');
%Tkaboom = table(Soln(:,1),Soln(:,2),Soln(:,3),Soln(:,4),Soln(:,5),Soln(:,6),'VariableNames',{'a','b','c','d','e','f'})
end

function F_obj_ret = OBJFUNC(x)
    F_obj_ret = abs((x(1) + 2*x(2) + 3*x(3) + 4*x(4) + 5*x(5) + 6*x(6)) - 240);
end

function Chromosome_ret = CROSSOVER(PChromosome1,PChromosome2,cp_rec)
    Chromosome_ret = PChromosome1;
    for i=cp_rec+1:6    %crossovers after the crossover point
        Chromosome_ret(i) = PChromosome2(i);
    end
end