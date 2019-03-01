function simanneal_ex5
%%
%VELASCO, Gimel David F.
%2012-58922
%Cmsc 191
%Simulated Annealing
%Exercise 5
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
%Objective: The Objective Function must be very close if not equal to 0 (Minimization)
%%
clc;    %for ease of documentation
tic;
%%%%%%%%%INITIALIZATION%%%%%%%%%%
prob_num = 3;                           %1 - Problem #1 | 2 - Problem #2 | o.w. - Problem #3
x_sol = 10*(rand(1)-0.5);               %initial guess of the solution x
cooling_sched = [100 50 40 30 20 10 5 2 1 0.5 0.25 0.1 0.001 0.0001 0.00001];    %initializes the cooling schedule
num_sched = numel(cooling_sched);       %gets the size of the cooling schedule
num_neigh = 10000;                      %initializes the size of the random neighbors
for i=1:num_sched
    T = cooling_sched(i);               %sets the value of the temperature T
    rand_neigh = 10*(rand(1,num_neigh)-0.5);     %random neighbors in [-1,1]
    for j=1:num_neigh
        x_tmp = rand_neigh(j);
        if E(x_tmp,prob_num) < E(x_sol,prob_num)
            x_sol = x_tmp;   
        elseif E(x_tmp,prob_num) > E(x_sol,prob_num)
            delta = E(x_tmp,prob_num) - E(x_sol,prob_num);
            p = P(delta,T);
            q = rand(1);
            if q <= p
                x_sol = x_tmp; 
            end
        end
    end
end

%SOLUTION
clc;
fprintf('=====================================================================\n\t\t\t\tSOLUTION\n\n');
fprintf('With the Objective Value of %.16f\nand Total Runtime of %f seconds\n',E(x_sol,prob_num),toc);
fprintf('The Root for Problem #%d is x = %.16f\n',prob_num,x_sol);
fprintf('=====================================================================\n');
%%
end

function F_obj_ret = E(x,prob_num)
    if prob_num == 1
        F_obj_ret = x - cos(x);
    elseif prob_num == 2
        F_obj_ret = exp(-x)*(x-2);
    else
        F_obj_ret = x^2 - x - 12;
    end
    F_obj_ret = abs(F_obj_ret);
end

function p_ret = P(delta,T)
    p_ret = exp(-(delta)/T);
end