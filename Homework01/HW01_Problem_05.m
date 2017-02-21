function [output] = HW01_Problem_05(input)
%function [output] = HW01_Problem_05(input)
%
%This function determines whether a number is prime or not
%
%Input: input: a number that the user enters
%
%Output: output: displays a message saying if the number entered is prime

%--------------------error check------------------------------

if nargin ~= 1
      error('This function requires one argument.');
end

%--------------------set up-----------------------------------

remainder = 1;  % variable that holds the remainder of (input / all numbers less than input and > 1)
IsNotPrime = false; % will be set to true if remainder ever equals 0

%--------------------calculate--------------------------------

if (input < 2) % 1 is not prime
    IsNotPrime = true;
else
    
for i = 2:input-1 % divides input by numbers from (2 to input - 1)
    remainder = rem(input, i); % store remainder value
    
    if (remainder == 0) % if input can be divided by a number with no remainder then it is not prime
    IsNotPrime = true; 
    end
end
end

%--------------------display results------------------------------

if (IsNotPrime == true) % display results with a fun easy to read message!!
    output = sprintf('The number %d is not a prime number!', input);
else
    output = sprintf('The number %d is a prime number!', input);
end
