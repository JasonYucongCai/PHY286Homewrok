function [x] = HW3_problem1(x0, alpha, beta, n, N_iteration)
%function [x] = HW3_problem1(x0, alpha, beta, n, N_iteration)
%
%This function produces an iterative map to generate psuedo-random numbers
%using the function x(i + 1) = mod(alpha * x(i) * 2^n + beta, 2^n) / (2^n)
%Then a plot is created, plotting x vs. N_iteration.
%
%Inputs: x0: inital x value, should be between 0-1
%        alpha: constant, integer between 0-15
%        beta: constant: integer between 0-3
%        n: constant: should be of the form 2^m, where m is an integer >= 3
%        N_iteration: number of iterations to go through (e.g. how many x 
%                     values should be output)
%
%Output: x: a vector of length N_iteration that contains values between 0-1
%
%NOTE: good alpha values seem to be those that are odd, excluding 1
%      there seem to be a few exceptions to this based on x0
%      ex) x0 = .001, alpha = 15, beta = 1, n = 8 
%       
%      even alpha values seem to cause the fuction to converge upon a
%      particular value and are thus not good choices for generating random
%      numbers
%
%      good beta values seem to be any beta that is not zero
%      beta = 0 seems to produce the most cases where random numbers are
%      not generated

%-----------------------initial setup--------------------------------------

%establish length of 'x' vector
x = zeros(N_iteration, 1);

%establish initial 'x'
x(1) = mod(alpha * x0 * 2^n + beta, 2^n) / (2^n);

%-----------------------calculate------------------------------------------

%calculate next 'x' based on the previous one
for count = 1:N_iteration - 1
    x(count + 1) = mod(alpha * x(count) * 2^n + beta, 2^n) / (2^n);
end

%-----------------------plot-----------------------------------------------

%plot
plot(1:N_iteration, x, '-ok');
xlabel('Number of iterations, N');
ylabel('Output, x');
title(sprintf('Iterative map with x_{0} = %.3f, alpha = %d, beta = %d, n = %d', x0, alpha, beta, n));

%display
disp(' ');
fprintf('Good choices (with a few exceptions):');
fprintf('\n\t Odd Alpha > 1');
fprintf('\n\t Beta > 0');
disp(' ');
fprintf('Bad choices (with a few exceptions):');
fprintf('\n\t Even Alpha > 1');
fprintf('\n\t Beta = 0');
disp(' ');

