function [output] = HW01_Problem_06(n)
%function [output] = HW01_Problem_06(n)
%
%This function lists the first 'n' prime numbers
%
%Input: n: determines how many primes to print
%
%Output: output: will display the first 'n' prime numbers

%--------------------error check------------------------------

if nargin ~= 1
      error('This function requires one argument.');
end

if n <= 0
    error('This function should print at least one number.');
end

%--------------------set up-----------------------------------

potential_prime = 3; %starts at 3 since since output starts at 2 (first prime) and thus for any n > 1 will thus start at 3
counter = 2; %will keep track of how many primes have been displayed, starts at 2 since for n > 1 two values will be added to output initially
output = 2; % starts at 2 since 2 is the first prime

%--------------------calculate--------------------------------

while counter <= n
    
    is_prime = true; %assume that the number is prime until proven otherwise
   
    for i = 2:potential_prime - 1 %check to see if it's prime
    if rem(potential_prime,i) == 0 %it is now proven that the number is not prime
        is_prime = false;
    end
    end
    
    if is_prime == true
    output = [output; potential_prime]; %potential_prime is in fact a prime number
    counter = counter + 1;
    end
    
    potential_prime = potential_prime + 1; %increment the number so as to continue to check higher numbers until the counter cond. is met
end
   

