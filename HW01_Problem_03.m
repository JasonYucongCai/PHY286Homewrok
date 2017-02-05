function[]=HW01_Problem_03(n)
% The function that takes a single integer (n) as input and returns the 
% nth Fibonacci number without using any the function.
%
% The function using the modern definition include the initial integer 0.
% the program will return the n't Fibonacci number. n must be an integer.
% Yuocng Cai
if nargin ~=1 %nargin -- new concept! see help nargin
  disp('Need exactly one inputs. See help');
    %error('Need exactly two inputs. See help');
  return;
end
%Input comfirm, Structure comfirm
%-------------------------------------------------------------------------
switch n
    case 0
        x=0;
    case 1;
        x=1;
    otherwise
        x=1;
        y=0;
%-------------------------------------------------------------------------      
        if n<0
            new=(-n);
        else
            new=n;
        end
%resolve the negative integer
%-------------------------------------------------------------------------
        for i=1:1:(new-1)
            %y=a{n-1}
            z=x;% store the an;
            x=x+y;% a{n-1}+an comput the an+1
            y=z; %y=an;
            
        
       % fprintf('at round %i, the value of x=%i y=%i\n ',i+1,x,y);
        end
% find the absolute value of the Fibonacci number.
%-------------------------------------------------------------------------
        if n<0
            switch 1
               
                case mod(n,2)==0
                    x=(-x);
                    disp('command 1 has been excued');
                otherwise
               %x= x;
               %disp('command 2 has been excued');
            end
       
        end
% find the sign of the Fibonacci number
%-------------------------------------------------------------------------        
end 

fprintf('the %i Fibonacci number is  %i\n',n,x);
% print out the number
%-------------------------------------------------------------------------
