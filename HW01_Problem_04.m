function[]=HW01_Problem_04(n)
% The function accepts a single integer (n) as input and returns the nth 
% Fibonacci number using the closed form expression (called "Binet's 
% formula") in Wikipedia's page. Make the function print the %error
% between the closed-form computed value and the "true" nth Fibonacci 
% number in the series.
%
% Using the modern definition include the initial integer 0.
% The input must be itegers.
% the program will return the n't Fibonacci number.
%Yuocng Cai

%-------------------------------------------------------------------------
if nargin ~=1 %nargin -- new concept! see help nargin
  disp('Need exactly one inputs. See help');
    %error('Need exactly two inputs. See help');
  return;
end
%Input comfirm, Structure comfirm
%-------------------------------------------------------------------------

x=1/sqrt(5)*(   ((1+sqrt(5))/2)^n - ((1-sqrt(5))/2)^n  );
format short;
fprintf('the %i Fibonacci number is  %.f\n',n,x);
%-------------------------------------------------------------------------
