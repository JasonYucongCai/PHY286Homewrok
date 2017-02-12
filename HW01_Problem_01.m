function [str_output] = HW01_Problem_01(str_input, n)
%function [str_output] = HW01_Problem_01(str_input, n)
%
%This function will display a string input on a different line 'n' times
%
%Input: str_input: the string value input 
%Input: n: the number of times str_input will be displayed
%
%Output: str_output: the input to be displayed on 'n' lines

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
       str_output = [str_output; str_input]; % create matrix of inputs
  end 
     

  