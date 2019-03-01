function runGA(popsize,maxgens,pc,pm,tol,runs,tpl,ulb,test_func)

Chromosome = zeros(tpl);
for i=1:popsize
    Chromosome(i,:) = 2*ulb*(rand(1,tpl)-0.5);  %initializing first generation
end
for generation=1:maxgens   %LOOPS STEPS 2 to 5
%STEP 2: Selection
%Selection by Tournament
for i=1:popsize/4   %Tournament process
    red = abs(floor(popsize*rand(1))) + 1;  %two random chromosomes to compete
    blu = abs(floor(popsize*rand(1))) + 1;
    if evalFitness(Chromosome(red,:),tpl,test_func) > evalFitness(Chromosome(blu,:),tpl,test_func)    %competition
        Chromosome(red,:) = Chromosome(blu,:);      %Blue Wins the tournament and overwrites Red
    else
        Chromosome(blu,:) = Chromosome(red,:);      %Red Wins the tournament and overwrites Blue
    end
end
%%
%Pre-Evaluation
%%
%totalfitness = 0;
%F_obj = zeros(1);
%for i=1:popsize
%   F_obj(i) = evalFitness(Chromosome(i,:),tpl,test_func);
%    totalfitness = totalfitness + F_obj(i);
%end
%%
%Selection by Proportion
R = rand(1,popsize);
%R = zeros(1);
%for i=1:popsize
%    R(i) = F_obj(i)/totalfitness;
%end
k = zeros(1);
PChromosome = zeros(tpl);
cp = zeros(1);
ctr = 0;    %holds number of parents
for i=1:popsize
    if R(i) > pc
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
%STEP 4: Mutation
%Per Chromosome mutation
mu = round(pm*popsize); %#ofchromosomestomutate = mutationrate*populationsize
for i=1:mu
    cngn = abs(floor(popsize*rand(1))) + 1; %random popsize number
    q = evalFitness(Chromosome(cngn,:),tpl,test_func);
    if q < 1
        Chromosome(cngn,:) = 2*ulb*q*(rand(1,tpl)-0.5);   %mutation
    else
        Chromosome(cngn,:) = 2*ulb*(rand(1,tpl)-0.5);
    end
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
%STEP 5: Post-Evaluation
F_obj = zeros(1);
for i=1:popsize
   F_obj(i) = evalFitness(Chromosome(i,:),tpl,test_func);
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
fprintf('====================GENETIC ALGORITHM=============================\n');
%fprintf('Final Generation\n');
%Chromosome = SORT(Chromosome,popsize,tpl,test_func);
%disp(Chromosome)
fprintf('Generation Number: %d\nPopulation Size: %d\tCrossover Rate: %.2f\tMutation Rate: %.2f\n',maxgens,popsize,pc,pm);
%fprintf('\nThe Fittest Chromosome is\n');
%disp(Chromosome(fi,:)')
%fprintf('Using Tournament Selection\n');
fprintf('Fitness Function Value of \t\t\t %.16f\nTotal Runtime of \t\t\t\t %f seconds\n',fittest,toc);
fprintf('Global Minimum is at:\n');
%disp(Chromosome(fi,:))
for p=1:11
    fprintf('a%d = %f\n',p,Chromosome(fi,p));
end
for p=1:10
    fprintf('b%d = %f\n',p,Chromosome(fi,p+11));
end
fprintf('==================================================================\n');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end