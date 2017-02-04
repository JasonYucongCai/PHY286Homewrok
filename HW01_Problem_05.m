function[]=HW01_Problem_05(n)
%single input integer and determins if the input integer is prime.
%The program will determine weather the integer is prime or not.
%Yuocng Cai
if n<0
    n=(-n);
end

A=factor(n);

if size(A)==1
fprintf('%i is a prime number.\n',n);
else
    fprintf('%i is not a prime number.\n',n);
end
