function [a_sol, energy] = runSA(x,y,xmax,xmin,cooling_ratio,tpl,test_func,initial_T,cooling_stop,ulb,num_neigh,runs,ci,ca,cb)
%%%%%%%%%%%%%%%%%%%%%%%%%SIMULATED ANNEALING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a_sol = 2*ulb*(rand(1,tpl)-0.5);               %initial guess of the solution a
cooling_sched = zeros(1);               %pre-allocation for speed
cooling_sched(1) = initial_T;                 %initializes the cooling schedule T0

tic;
sched = 1;                                  %index
iteration(1) = sched;
while cooling_sched(sched) > cooling_stop     %iteration will stop if the cooling temperature reached less than 0.00000001
    T = cooling_sched(sched);               %sets the value of the temperature T
    for j=1:num_neigh
        %r  = (cooling_ratio)^sched;             %is used so that the randomness of selecting a neighbor becomes narrower
        a_tmp = 2*ulb*(rand(1,tpl)-0.5);    %selects a random neighbor for comparison. with a decreasing amount of randomness
        y_sim_sol = getSimulated(a_sol,x,tpl,test_func);
        y_sim_tmp = getSimulated(a_tmp,x,tpl,test_func);
        if evalFitness(y_sim_tmp,y,xmax,xmin) < evalFitness(y_sim_sol,y,xmax,xmin)  %if the neighbor is better, change the solution
            a_sol = a_tmp;
        elseif evalFitness(y_sim_tmp,y,xmax,xmin) > evalFitness(y_sim_sol,y,xmax,xmin)  %if not, change the solution if it is lucky
            delta = evalFitness(y_sim_tmp,y,xmax,xmin) - evalFitness(y_sim_sol,y,xmax,xmin);
            p = P(delta,T);
            q = rand(1);
            if q <= p
                a_sol = a_tmp; 
            end
        end
    end
    fprintf('Run %d | Energy: %.16f | Time: %.2f s | Temperature: %.16f\n',runs,evalFitness(y_sim_sol,y,xmax,xmin),toc,cooling_sched(sched));
    cooling_sched(sched+1) = cooling_implementation(ci,ca,cb,cooling_ratio,sched,T);
    sched = sched+1;
    iteration(sched) = sched;

end

%SOLUTION
energy = evalFitness(y_sim_sol,y,xmax,xmin);
plot(iteration,cooling_sched)
xlabel('Iteration')
ylabel('Temperature')
%fprintf('====================LONG EXAM 1 PART 2============================\n');
%fprintf('==================SIMULATED ANNEALING=============================\n');
%fprintf('Energy: %.16f Time: %.2f seconds Temperature: %.16f\n',evalFitness(y_sim_sol,y,xmax,xmin),toc,cooling_sched(sched));
%fprintf('Approx Coefficients are:\n');
%disp(a_sol)
%for p=1:11
%    fprintf('a%d = %f\n',p,y_sim_sol(p));
%end
%for p=1:10
%    fprintf('b%d = %f\n',p,y_sim_sol(p+11));
%end
%fprintf('==================================================================\n');
end