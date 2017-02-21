function [nth_term] = HW01_Problem_04(n)
%function [nth_term] = HW01_Problem_04(n)
%
%This function computes a value of the closed form of the fibonacci
%sequence and compares it to the value of the regular fibonacci sequence
%
%Input: n: the term number of the fibonacci sequence desired
%
%Output: nth_term: the value of term n
%
%NOTE: the output displays better with a semicolon at the end of the
%command

%--------------------error check------------------------------

if nargin ~= 1
      error('This function requires one argument.');
end

if n <= 0
    error('This function should print at least one number.');
end

%--------------------calculate--------------------------------

nth_term = ((1+sqrt(5))^n - (1-sqrt(5))^n) / (2^n * sqrt(5)); % calculates the closed form of the fibonacci sequence
percent_error = abs((HW01_Problem_03(n) - nth_term) / HW01_Problem_03(n) * 100); % percent error

%--------------------display results------------------------------

disp(' ');
fprintf('The binet term is: %10.2f\n', nth_term);
disp(' ');
fprintf('The percent error between the closed form and the regular sequence is: %10.2f\n', percent_error);
