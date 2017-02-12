function [str_output] = HW01_Problem_02(str_input, n)
%function [str_output] = HW01_Problem_02(str_input, n)
%
%This function will display a string input on the same line 'n' times
%
%Input: str_input: a string that the user enters
%Input: n: the number of times the string is to be repeated
%
%Output: str_output: a vector that holds all of the string values

%--------------------error check------------------------------

if nargin ~= 2
      error('This function requires two arguments. Please use a string as the first input and an integer as the second!');
end

if n <= 0
    error('This function should print at least one string.');
end

%--------------------calculate--------------------------------

str_output = str_input; % initialize 

for i = 1:n-1
    str_output = [str_output ' ' str_input]; % create matrix of inputs with a space between each string
end 
     

