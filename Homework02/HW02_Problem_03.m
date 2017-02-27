function [prey, predator, time] = HW02_Problem_03(alpha, beta, delta, gamma, prey_0, predator_0, dt, t_max, choice)
%function [prey, predator, time] = HW03_Problem_03(alpha, beta, delta, gamma, prey_0, predator_0, dt, t_max, choice)
%
%This function does 2 things.  First, it plots the number of prey and
%predator as a function of time.  Second, it plots predator vs prey as a
%phase plot
%
%Input: alpha: a constant
%       beta: a constant
%       delta: a constant
%       gamma: a constant
%       prey_0: initial number of prey present
%       predator_0: initial number of predators present
%       dt: time-step
%       t_max: the time interval over which the calculations will be made
%       choice: should be '1' for the Euler method, '2' for R-K method, any
%               other value results in a error
%
%Output: prey: number of prey present
%        predator: number of predators present
%        time: a vector containing times at which the number of prey and
%              predators was evaluated

%----------------------------error check-----------------------------------

if nargin ~= 9
    error('This function requires 9 inputs');
end

%---------------------error check------------------------------------------

if choice ~= 1 && choice ~= 2
    error('This function will use either Eulers method or the R-K method.  Please input "1" for Eulers method or "2" for the R-K method.');
end

%---------------------initial values---------------------------------------

prey(1) = prey_0;
predator(1) = predator_0;

time = [0:dt:t_max];

%--------------------calculate---------------------------------------------

for count = 1:length(time) - 1
    
    %Euler method
    if choice == 1
    prey(count + 1) = prey(count) + dt.*(alpha.*prey(count) - beta.*prey(count).*predator(count));
    predator(count + 1) = predator(count) + dt.*(delta.*prey(count).*predator(count) - gamma.*predator(count));
    
    %R-K method
    else
        prey_prime = prey(count) + .5.*dt;
        predator_prime = predator(count) + .5.*dt.*(delta.*prey(count).*predator(count) - gamma.*predator(count));
        
        prey(count + 1) = prey(count) + dt.*(alpha.*prey_prime - beta.*prey_prime.*predator_prime);
        predator(count + 1) = predator(count) + dt.*(delta.*prey_prime.*predator_prime - gamma.*predator_prime);
    end
    
end

%---------------------plot-------------------------------------------------

subplot(2,1,1); plot(time, prey, time, predator);
title('Prey & Predator vs. Time');
xlabel('Time (dimensionless)');
ylabel('Population of Prey & Predator');
legend('Prey','Predator');

subplot(2,1,2); plot(prey, predator, 'g');
title('Phase Plot, Predator vs. Prey');
xlabel('Prey');
ylabel('Predator');


