function[]=HW01_Problem_07(P,R,T)
% The function accepts three inputs: the initial principal (P), the annual
% interest rate (R) and total time of deposit (T). The function should 
% appropriately display the simple interest  (calculated as P*R*T) as well
% as the compound interest, for the time specified.
% Calculate the compound interest using a loop
%
% The function will display the simple and compound interest, the input
% should include three factors, in the order of principal, the annual
% interest rate and the total time of deposit.
% For example:  HW01_Problem_07(Principal,Rate,Time)
% The function assume the time is round down. For example, 2.9 years count 
% as 2 eyars.
%
% This code also satisfied the Extra Credit, which accept the row vectors
% as the input.
% Yucong Cai




%hello
%---------------------------------------------------------------------
if nargin ~=3 %nargin -- new concept! see help nargin
  disp('Need exactly three inputs. See help');
    %error('Need exactly two inputs. See help');
  return;
end
%input confirm, structure comfirm.
%-------------------------------------------------------------------------
N=size(P);
N=N(1,2);

for n=1:1:N
    if P(1,n)<0
    disp('principal, interest rate, and the interest must be bigger than 0');
return;
    end
end



%the number of input(Principal)
%-------------------------------------------------------------------------

if (T<0)||(R<0)
    disp('principal, interest rate, and the interest must be bigger than 0');
return;
end

if T<1;
disp('there is no interest earned');
return;
end    
%input confirm, legal comfirm.
%-------------------------------------------------------------------------


for n=1:1:N 
   
si(1,n)=P(1,n)*R*floor(T);
%disp(si(1,n));

fprintf('the simple interest of the first member is %.2f.\n',si(1,n));


end
%simple interest
%-------------------------------------------------------------------------
fprintf('\n\n');

for n=1:1:N 
    ci(1,n)=0;
    
    for a=1:1:floor(T)
        
    ci(1,n)=(ci(1,n)+ P(1,n))*R;

    end
        %disp(si(1,n));
    fprintf('the compound interest of the first member is %.2f.\n',ci(1,n));


end





%compound interest
%-------------------------------------------------------------------------
