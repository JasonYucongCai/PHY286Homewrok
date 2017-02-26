
function[]=HW02_Problem_03(prey,predator,T,dt,ER)

%check
%-------------------------------------------------------------------------
if nargin ~=5 %nargin -- new concept! see help nargin
  disp('Need exactly five inputs. See help');
    %error('Need exactly two inputs. See help');
  return;
end
%Input comfirm, Structure comfirm
%-------------------------------------------------------------------------
if prey<1 || predator<1 || dt>T || dt<0
 %   fprintf('Check your input, prey or predator must be bigger or equal\n');
    error('Check your input, prey or predator must be bigger or equal to 1, T>=dt>0.');
end



a=1/3;b=4/3;
v=1/3;u=4/3;

x(1)=prey;
y(1)=predator;
n=floor(T/dt);


switch ER
    case 'Euler'
        disp('Using Euler methods');
        
        for i=1:n-1 %where i was the index
        dx(i)=a*x(i)-b*y(i);
        dy(i)=v*x(i)-u*y(i);
        x(i+1)=x(i)+dx(i)*dt;
        y(i+1)=y(i)+dy(i)*dt;
        
        
        if y(i+1)<=0
            y(i+1)=0;
        end
        if x(i+1)<=0
            x(i+1)=0;
        end
        
        
        end
                
    case 'RK'
        disp('Using RK methods');
      
        for i=1:n-1 
        rkx=x(i)+(1/2)*(a*x(i)-b*y(i))*dt;
        rky=y(i)+(1/2)*(v*x(i)-u*y(i))*dt;
        x(i+1)=x(i)+(a*rkx-b*rky)*dt;
        y(i+1)=y(i)+(v*rkx-u*rky)*dt;
        
        if y(i+1)<=0
            y(i+1)=0;
        end
        if x(i+1)<=0
            x(i+1)=0;
        end
        
        end
        
    
    otherwise
        error('Check your input choice Euler or RK');

end
time=linspace(0,n*dt,n);

%Input comfirm, legal adjusting.
%-------------------------------------------------------------------------
figure
subplot(2,1,1); 
plot(time,x,'b.-',time,y,'r.-');
title('Prey & Predator vs. Time');
xlabel('Time (dimensionless)');
ylabel('Population of Prey & Predator');
legend('Prey','Predator');

subplot(2,1,2); 
plot(x, y);
title('Phase Plot, Predator vs. Prey');
xlabel('Prey');
ylabel('Predator');

end
%{
function[]=HW02_Problem_03(prey,predator,T,dt,ER)

%check
%-------------------------------------------------------------------------
if nargin ~=5 %nargin -- new concept! see help nargin
  disp('Need exactly five inputs. See help');
    %error('Need exactly two inputs. See help');
  return;
end
%Input comfirm, Structure comfirm
%-------------------------------------------------------------------------
if prey<1 || predator<1 || dt>T || dt<0
 %   fprintf('Check your input, prey or predator must be bigger or equal\n');
    error('Check your input, prey or predator must be bigger or equal to 1, T>=dt>0.');
end



a=1/3;b=4/3;
v=1/3;u=4/3;

x(1)=prey;
y(1)=predator;
n=floor(T/dt)


switch ER
    case 'Euler'
        disp('Using Euler methods');
        
        for i=1:n-1 %where i was the index
        dx(i)=a*x(i)-b*y(i);
        dy(i)=v*x(i)-u*y(i);
        x(i+1)=x(i)+dx(i)*dt;
        y(i+1)=y(i)+dy(i)*dt;
        end
                
    case 'RK'
        disp('Using RK methods');
      
        for i=1:n-1 
        rkx=x(i)+(1/2)*(a*x(i)-b*y(i))*dt;
        rky=y(i)+(1/2)*(v*x(i)-u*y(i))*dt;
        x(i+1)=x(i)+(a*rkx-b*rky)*dt;
        y(i+1)=y(i)+(v*rkx-u*rky)*dt;    
        end
        
    
    otherwise
        error('Check your input choice Euler or RK');

end
time=linspace(0,n*dt,n)

%Input comfirm, legal adjusting.
%-------------------------------------------------------------------------
figure
subplot(2,1,1); 
plot(time,x,'b.-',time,y,'r.-');
title('Prey & Predator vs. Time');
xlabel('Time (dimensionless)');
ylabel('Population of Prey & Predator');
legend('Prey','Predator');

subplot(2,1,2); 
plot(x, y);
title('Phase Plot, Predator vs. Prey');
xlabel('Prey');
ylabel('Predator');

end


function [prey, predator, time] = HW03_Problem_03(alpha, beta, delta, gamma, prey_0, predator_0, dt, t_max, choice)
%function [prey, predator, time] = HW03_Problem_03(alpha, beta, delta, gamma, prey_0, predator_0, dt, t_max, choice)
%
%This function does 2 things.  First, it plots the number of prey and
%predator as a function of time.  Second, it plots predator vs prey as a
%phase plot
%
%Input: alpha: a constant
%       beta: a constant
%       delta: a constant
%       gamma: a constant
%       prey_0: initial number of prey present
%       predator_0: initial number of predators present
%       dt: time-step
%       t_max: the time interval over which the calculations will be made
%       choice: should be '1' for the Euler method, '2' for R-K method, any
%               other value results in a error
%
%Output: prey: number of prey present
%        predator: number of predators present
%        time: a vector containing times at which the number of prey and
%              predators was evaluated

%---------------------error check------------------------------------------

if choice ~= 1 && choice ~= 2
    error('This function will use either Eulers method or the R-K method.  Please input "1" for Eulers method or "2" for the R-K method.');
end

%---------------------initial values---------------------------------------

prey(1) = prey_0;
predator(1) = predator_0;

time = [0:dt:t_max];

%--------------------calculate---------------------------------------------

for count = 1:length(time) - 1
    
    %Euler method
    if choice == 1
    prey(count + 1) = prey(count) + dt.*(alpha.*prey(count) - beta.*prey(count).*predator(count));
    predator(count + 1) = predator(count) + dt.*(delta.*prey(count).*predator(count) - gamma.*predator(count));
    
    %R-K method
    else
        prey_prime = prey(count) + .5.*dt;
        predator_prime = predator(count) + .5.*dt.*(delta.*prey(count).*predator(count) - gamma.*predator(count));
        
        prey(count + 1) = prey(count) + dt.*(alpha.*prey_prime - beta.*prey_prime.*predator_prime);
        predator(count + 1) = predator(count) + dt.*(delta.*prey_prime.*predator_prime - gamma.*predator_prime);
    end
    
end

%---------------------plot-------------------------------------------------

subplot(2,1,1); plot(time, prey, time, predator);
title('Prey & Predator vs. Time');
xlabel('Time (dimensionless)');
ylabel('Population of Prey & Predator');
legend('Prey','Predator');

subplot(2,1,2); plot(prey, predator);
title('Phase Plot, Predator vs. Prey');
xlabel('Prey');
ylabel('Predator');

%}
