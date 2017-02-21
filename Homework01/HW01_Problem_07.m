function [simple_interest, compound_interest] = HW01_Problem_07(P, R, T);
%function [simple_interest, compound_interest] = HW01_Problem_07(P, R, T)
%
%This function will return the simple and compound interest based on the
%given inputs.  
%
%Input: P: the principle amount(s)
%Input: R: the annual interest rate
%Input: T: the total time of deposit
%
%Output: simple_interest: interest calculated based on the initial P,R,T
%Output: compound_interest: interest calculated based on each consecutively
%calculated P and R and T
%
%NOTE: the output displays better with a semicolon at the end of the
%command

%--------------------error check------------------------------

if nargin ~= 3
    error('You must enter three arguments.  They must be numeric.');
end

%--------------------simple interest------------------------------

simple_interest = P .* R .* T; % calculates simple interest and allows for use of vectors 

%--------------------coumpound interest------------------------------

P_final = P; % holder variable for P so that it can be used to calculate the interest later on

for i = 1:T % calculates the total money after T years
   P_final = (1+R).*P_final; 
end
compound_interest = P_final - P; %subtracts initial principle from the final principle to calculate the interest

%--------------------display results------------------------------

disp(' ');
fprintf('The simple interest is: %10.2f\n', simple_interest);
disp(' ');
fprintf('The compound interest is: %10.2f\n', compound_interest);





