function[]=HW01_Problem_04(n)
% using the modern definition include the initial integer 0.
% the program will return the n't Fibonacci number.
%Yuocng Cai

%-------------------------------------------------------------------------
if nargin ~=2 %nargin -- new concept! see help nargin
  disp('Need exactly two inputs. See help');
    %error('Need exactly two inputs. See help');
  return;
end
%-------------------------------------------------------------------------

x=1/sqrt(5)*(   ((1+sqrt(5))/2)^n - ((1-sqrt(5))/2)^n  );
format short;
fprintf('the %i Fibonacci number is  %.f\n',n,x);