function[]=HW01_Problem_05(n)
%The function that accepts a single input integer and determines if the 
%input integer is prime.
%
%Single input integer and determins if the input integer is prime.
%The program will determine weather the integer is prime or not.
%Yuocng Cai
%-------------------------------------------------------------------------
if nargin ~=1 %nargin -- new concept! see help nargin
  disp('Need exactly one inputs. See help');
    %error('Need exactly two inputs. See help');
  return;
end
%Input comfirm, Structure comfirm
%-------------------------------------------------------------------------

if n<0
    n=(-n);
end
%Input comfirm, legal adjusting.
%-------------------------------------------------------------------------
A=factor(n);
%The size of the number.
%-------------------------------------------------------------------------
if size(A)==1
    fprintf('%i is a prime number.\n',n);
else
    fprintf('%i is not a prime number.\n',n);
end
%-------------------------------------------------------------------------
