function[]=HW01_Problem_06(n,firstn)
%the first argument as the integer that were to be valued, the second
%argument as the returns the determinted number of of primes. For example,
%HW01_Problem_06(256,2)
%first value as the input value and second value as to show the number of
%the first prime numbers.
%Yuocng Cai
%-------------------------------------------------------------------------
if nargin ~=2 %nargin -- new concept! see help nargin
  disp('Need exactly two inputs. See help');
    %error('Need exactly two inputs. See help');
  return;
end
%-------------------------------------------------------------------------
if n<0
    n=(-n);
end

A=factor(n);
N=size(A);
N=N(1,2);
%N is the total number of sum of powers of the primes.(How many numbers)
%fprintf('\nthe total number of numbers is %i\n\n',N);

bp=1;%the position to store the number.
co=0;%stand for count.

%set up.
%------------------------------------------------------------------------

bprime(1,bp)=A(1,1);% set up the initial 
nup(1,bp)=0;

for i=1:1:(N-1)
  co=co+1;%count the power of prime
   if  A(1,i)~=A(i+1)

        nup(1,bp)=co;% the power of prime;
        bp=bp+1;
        bprime(1,bp)=A(i+1);
        co=0;
    end
   
end
%find the list of all prime and (bp-1)'s power of the prime
%-------------------------------------------------------------------------


switch bp
    case 1
        nup(1,bp)=co+1;
   % case 2
    %    nup(1,2)=N-nup(1,1);
        
    otherwise
    nup(1,bp)=N;
    for i=1:1:(bp-1)
        nup(1,bp)=nup(1,bp)-nup(1,i);
    end

end

 
%to find the power of the last prime
%-------------------------------------------------------------------------
if bp<firstn
    finaln=bp;
else
    finaln=firstn;
end

fprintf('\nThe first %i prime factor(s) is(are):\n',finaln);
for i=1:1:finaln
    fprintf('the prime %i with power %i\n', bprime(1,i),nup(1,i))
end

    %print out
%-------------------------------------------------------------------------


if size(A)==1
fprintf('\n%i is a prime number.\n',n);
else
    fprintf('\n%i is not a prime number.\n',n);
end
