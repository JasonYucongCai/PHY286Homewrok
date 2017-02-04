function[]=HW01_Problem_03(n)
% using the modern definition include the initial integer 0.
% the program will return the n't Fibonacci number.
%Yuocng Cai
switch n
    case 0
        x=0;
    case 1;
        x=1;
    otherwise
        x=1;
        y=0;
        
        if n<0
            new=(-n);
        else
            new=n;
        end
        
        for i=1:1:(new-1)
            %y=a{n-1}
            z=x;% store the an;
            x=x+y;% a{n-1}+an comput the an+1
            y=z; %y=an;
            
        
       % fprintf('at round %i, the value of x=%i y=%i\n ',i+1,x,y);
        end
         
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
       
        
end 

fprintf('the %i Fibonacci number is  %i\n',n,x);