function [nth_term] = HW01_Problem_03(n)
%function [nth_term] = HW01_Problem_03(n)
%
%This function returns the nth term of the fibonacci sequence
%
%Input: n: determines which term to display from the fibonacci sequence
%
%Output: nth_term: the term to be displayed from the fibonacci sequence

%--------------------error check------------------------------

if nargin ~= 1
      error('This function requires one argument.');
end

%--------------------set up-----------------------------------

if n == 1 || n == 2 % the first two terms are 1 and 1
    nth_term = 1;
else % n > 2
    
n1 = 1; % starts fibonacci sequence off with the first three terms
n2 = 1;
nth_term = n1 + n2; 

%--------------------calculate--------------------------------

for i = 4:n  % calculates terms of the fibonacci sequence
   n1 = n2;
   n2 = nth_term;
   nth_term = n1 + n2;
end
end