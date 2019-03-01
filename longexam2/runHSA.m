function [a_sol, fittest] = runHSA(x,y,xmax,xmin,ulb,tpl,popsize,pc,cooling_ratio,num_neigh,initial_T,cooling_stop,test_func,runs,ci,ca,cb)
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
iteration(1) = sched;
while cooling_sched(sched) > cooling_stop     %iteration will stop if the cooling temperature reached less than 0.00000001
    T = cooling_sched(sched);               %sets the value of the temperature T
    %III.a. Do N/2 times
    for j=1:popsize/4
        %III.a.i. Select two parents at random
        red = abs(floor(popsize*rand(1))) + 1;  %two random chromosomes to compete
        blu = abs(floor(popsize*rand(1))) + 1;  %red and blu hold the index of the parents
        y_sim_red = getSimulated(Chromosome(red,:),x,tpl,test_func);
        y_sim_blu = getSimulated(Chromosome(blu,:),x,tpl,test_func);
        %III.a.ii. Generate two offsprings
        %%Recombination Operator (CROSSOVER)
        pc_trial = rand(1);
        if pc_trial < pc     %if trial made it in the crossover rate
            cp = floor(abs((tpl-1)*rand(1)))+1;   %random crossover point
            Child_Chromosome(1,:) = crossover(Chromosome(red,:),Chromosome(blu,:),cp,tpl);     %crossover red and blu
            Child_Chromosome(2,:) = crossover(Chromosome(blu,:),Chromosome(red,:),cp,tpl);     %they will have two children
            %%Neighborhood Operator (MUTATION)
            a_sol = Child_Chromosome(1,:);      %one child only will undergo mutation         
            for i=1:num_neigh
                a_tmp = 2*ulb*(rand(1,tpl)-0.5);    %selects a random neighbor for comparison. with a decreasing amount of randomness
                y_sim_sol = getSimulated(a_sol,x,tpl,test_func);
                y_sim_tmp = getSimulated(a_tmp,x,tpl,test_func);
                if evalFitness(y_sim_tmp,y,xmax,xmin) < evalFitness(y_sim_sol,y,xmax,xmin)  %if the neighbor is better, change the solution
                    a_sol = a_tmp;
                end
            end
            Child_Chromosome(1,:) = a_sol;           %will overwrite the child based on the neighborhood operator
            %%III.a.iii. Boltzman Trials
            y_sim_ch1 = getSimulated(Child_Chromosome(1,:),x,tpl,test_func);
            y_sim_ch2 = getSimulated(Child_Chromosome(2,:),x,tpl,test_func);
            
            E1 = evalFitness(y_sim_red,y,xmax,xmin) + evalFitness(y_sim_blu,y,xmax,xmin);
            E2 = evalFitness(y_sim_ch1,y,xmax,xmin) + evalFitness(y_sim_ch2,y,xmax,xmin);
            bp = boltzman(E1,E2,T);
            bp_trial = rand(1);
            if bp_trial >= bp
                %%III.a.iv. Overwrite Parents with the Trial Winner
                Chromosome(red,:) = Child_Chromosome(1,:);
                Chromosome(blu,:) = Child_Chromosome(2,:);
            end
            %%
        else                    %if the whole trial did not make it inside the crossover rate, it will have a tournament
            if evalFitness(y_sim_red,y,xmax,xmin) > evalFitness(y_sim_blu,y,xmax,xmin)    %competition
                Chromosome(red,:) = Chromosome(blu,:);      %Blue Wins the tournament and overwrites Red
            else
                Chromosome(blu,:) = Chromosome(red,:);      %Red Wins the tournament and overwrites Blue
            end
        end
    end
    fprintf('Run %d | Time: %.2f s | Temperature: %.16f\n',runs,toc,cooling_sched(sched));
    %III.b. Periodically Lower T
    cooling_sched(sched+1) = cooling_implementation(ci,ca,cb,cooling_ratio,sched,T);
    sched = sched+1;
    iteration(sched) = sched;
end
%%
%Post-Evaluation
y_sim_all = zeros(numel(y));
for i=1:popsize
    y_sim_all(i,:) = getSimulated(Chromosome(i,:),x,tpl,test_func);
end
F_obj = zeros(1);
for i=1:popsize
   F_obj(i) = evalFitness(y_sim_all(i,:),y,xmax,xmin);
end
%%
fittest = F_obj(1);
%fi = 1;
for i=1:popsize
    if fittest > F_obj(i)
       fittest = F_obj(i);
       %fi = i;
    end
end
%%
%SOLUTION
plot(iteration,cooling_sched)
xlabel('Iteration')
ylabel('Temperature')
%uncomment this section to see results
%%
%fprintf('==================HYBRID SIMULATED ANNEALING=======================\n');
%%
%fprintf('Population Size: %d\tCrossover Rate: %.2f\tMutation Rate: %.2f\n',popsize,pc,pm);
%fprintf('Minimum Energy:\t\t\t\t\t%.16f\nTotal Runtime:\t\t\t\t\t%f seconds\nFinal Temperature:\t\t\t\t%.16f\n',fittest,toc,cooling_sched(sched));
%fprintf('Global Minimum is at:\n');
%disp(Chromosome(fi,:))
%for p=1:11
%    fprintf('a%d = %f\n',p,Chromosome(fi,p));
%end
%for p=1:10
%    fprintf('b%d = %f\n',p,Chromosome(fi,p+11));
%end
%fprintf('===================================================================\n');
end