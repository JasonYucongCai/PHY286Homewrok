function [] = HW02_Problem_02(v0, theta, dt, B, m,range_approx_tolorance)
%function [x_normal, y_normal, x_airdrag, y_airdrag, x_density, y_density, Py] = HW02_Problem_02(v0, theta, dt, B, m)
%
%This function displays the trajectory of an object in three ways.  First,
%without taking into account air resistance or density.  Second, by taking
%into account only air resistance.  Third, by taking into account both air
%resistance and density.
%
%Input: v0: the initial velocity of the object (m/s)
%       theta: the angle at which the object is to be fired (degrees)
%       dt: time-step
%       B: drag coefficient
%       m: mass of the object
%       range_approx_tolorance:constant where the range apporximation was
%       allowed
%
%Output: x_normal: the x position of the object fired without taking into
%                  account air resistance or density
%        y_normal: the y position of the object fired without taking into
%                  account air resistance or density
%        x_airdrag: the x position of the object fired taking into account
%                   only air resistance
%        y_airdrag: the y position of the object fired taking into account
%                   only air resistance
%        x_density: the x position of the object fired taking into account
%                   both air resistance and density
%        y_density: the y position of the object fired taking into account
%                   both air resistance and density
%        Py: calculates density as a function of height (y_density)
                  
%-----------------------initial values-------------------------------------
 
 
%-------------------------------------------------------------------------
if nargin ~=6 %nargin -- new concept! see help nargin
  disp('Need exactly six inputs. See help');
    %error('Need exactly two inputs. See help');
  return;
end
%Input comfirm, Structure comfirm
%-------------------------------------------------------------------------
if v0<0 || dt<0 || B<0 || m<0 || range_approx_tolorance<0
    error('check your input, v0,dt, B, m,range_approx_tolorance could not be negative.');
end
 
if theta<=0|| theta>=(pi/2)
    error('check your input, theta must between 0 to pi/2 and not including 0 and pi/2 .');
    
end
 
%Input comfirm, legal adjusting.
%-------------------------------------------------------------------------
 
 
 
 
%position
x_normal(1) = 0;
y_normal(1) = 0;
 
x_airdrag(1) = 0;
y_airdrag(1) = 0;
 
x_density(1) = 0;
 
y_density(1) = 0;
 
%velocity
vx_normal(1) = v0.*cos(theta);
vy_normal(1) = v0.*sin(theta);
 
vx_airdrag(1) = v0.*cos(theta);
vy_airdrag(1) = v0.*sin(theta);
 
vx_density(1) = v0.*cos(theta);
vy_density(1) = v0.*sin(theta);
 
%constants
g = 9.8;
y_0 = 10000; %approximate value from book, calculated as (kB*T/m*g)
 
%pressure
p_0 = 1.225; %approximate value of air density at sea level
 
%--------------------------compute-----------------------------------------
%(adjusted)
%first calculation
%-------------------------------------------------------------------------
 
%projectile with neither drag nor density
 
counter_n = 1;
 
while y_normal(counter_n) >= 0
    %velocity
    vx_normal(counter_n + 1) = vx_normal(counter_n) + dt.*0;
    vy_normal(counter_n + 1) = vy_normal(counter_n) + dt.*(-g);
    
    %position
    x_normal(counter_n + 1) = x_normal(counter_n) + dt.*(vx_normal(counter_n));
    y_normal(counter_n + 1) = y_normal(counter_n) + dt.*(vy_normal(counter_n));
   
    %increment counter_n    
    
    counter_n = counter_n + 1;
end%-------------------without out air drag
 
%projectile with drag only
counter_d = 1;
while y_airdrag(counter_d) >= 0
   
    %velocity
    vx_airdrag(counter_d + 1) = vx_airdrag(counter_d) + dt.*0 + (-B./m).*vx_airdrag(counter_d).*sqrt(vx_airdrag(counter_d).*vx_airdrag(counter_d) + vy_airdrag(counter_d).*vy_airdrag(counter_d)).*dt;
    vy_airdrag(counter_d + 1) = vy_airdrag(counter_d) + dt.*(-g) + (-B./m).*vy_airdrag(counter_d).*sqrt(vx_airdrag(counter_d).*vx_airdrag(counter_d) + vy_airdrag(counter_d).*vy_airdrag(counter_d)).*dt;
    
    %position
    x_airdrag(counter_d + 1) = x_airdrag(counter_d) + dt.*(vx_airdrag(counter_d));
    y_airdrag(counter_d + 1) = y_airdrag(counter_d) + dt.*(vy_airdrag(counter_d));
    
    %increment counter_d
    counter_d = counter_d + 1;
end%---------------------with air drag
 
%projectile with both drag and density
counter_dd = 1;
while y_density(counter_dd) >= 0
    %solve for density
    Py(counter_dd) = p_0.*exp(-y_density(counter_dd)./y_0);
    
    %velocity
    vx_density(counter_dd + 1) = vx_density(counter_dd) + dt.*0 + (-B./m).*(Py(counter_dd)./p_0).*vx_density(counter_dd).*sqrt(vx_density(counter_dd).*vx_density(counter_dd) + vy_density(counter_dd).*vy_density(counter_dd)).*dt;
    vy_density(counter_dd + 1) = vy_density(counter_dd) + dt.*(-g) + (-B./m).*(Py(counter_dd)./p_0).*vy_density(counter_dd).*sqrt(vx_density(counter_dd).*vx_density(counter_dd) + vy_density(counter_dd).*vy_density(counter_dd)).*dt;
    
    %position
    x_density(counter_dd + 1) = x_density(counter_dd) + dt.*(vx_density(counter_dd));
    y_density(counter_dd + 1) = y_density(counter_dd) + dt.*(vy_density(counter_dd));
    
    %increment counter_dd
    counter_dd = counter_dd + 1;
end%-------------------------with air drag and air density
 
%compute range of projectile
r_range_y_compare = - (y_normal(end - 1) / y_normal(end));
range_normal = (x_normal(end - 1) + r_range_y_compare * x_normal(end)) / (r_range_y_compare + 1);
 
r_airdrag = - (y_airdrag(end - 1) / y_airdrag(end));
range_airdrag = (x_airdrag(end - 1) + r_airdrag*x_airdrag(end)) / (r_airdrag + 1);
 
r_density = - (y_density(end - 1) / y_density(end));
range_density = (x_density(end - 1) + r_density*x_density(end)) / (r_density + 1);
 
 
%-------------------------------------------------------------------------
%finding the maxmium
%know that the funciton Range about the launch angle is continus and
%differentiable, also known that the the second order differential of the
%function exists and smaller than 0 all the time. though 0 an pi/2 was the
%limit point of the set, however, we know that it's equal to 0 and is not
%of the interst. We know that the maxmium exists, and is the strict local
%maxmium of the function.
%We use the fact that the series 1/2 +/- (1/2)^2 +/-(1/2)^3 ... was able 
%to converage to any point in R with a radius of 2.
%the out for loop with index h was used to increase the accuracy,
%effectively increaed the accuracy by 10e3 to 10e4 times, where h=0
%provided the initial references at pi/4;
%----------------------------------------------------------------------
%initial value set up
tic;
    theta2=pi/4;
    range_count = 1;
    %cleanup;
    vx_range=zeros();
    vy_range=zeros();
    x_range=zeros();
    y_range=zeros();
    vx_range(1) = v0.*cos(theta2);
    vy_range(1) = v0.*sin(theta2);
    y_range_maxium=0;
    x_range_maxium=-inf;
    theta_max2=theta2;
    x_range_minum=zeros();
    n=1;
    loop_count=0;
converage_converage_ratio=1/2;
%-----------------------------------------------------------------------
 
for h=0:2%round number, add the number to get more accuracy.
    x_range_minum(end)=0;
while (abs(x_range_maxium(end)-x_range_minum(end)))>=range_approx_tolorance  
    x_range_converage1=zeros();
    x_range_converage2=zeros();
    y_range_converage1=zeros();
    y_range_converage2=zeros();
    
    for k=1:1:2
        theta2=theta_max2+h*(((-1)^k)*((converage_converage_ratio)^(n))*(pi/4));        
    %cleanup;
    vx_range=zeros();
    vy_range=zeros();
    x_range=zeros();
    y_range=zeros();
    vx_range(1) = v0.*cos(theta2);
    vy_range(1) = v0.*sin(theta2);
    range_count = 1;
    %cleanup;
    
    while y_range(range_count) >= 0
     
    %solve for density
    Py(range_count) = p_0.*exp(-y_range(range_count)./y_0);
    
    %velocity
    vx_range(range_count + 1) = vx_range(range_count) + dt.*0 + (-B./m).*(Py(range_count)./p_0).*vx_range(range_count).*sqrt(vx_range(range_count).^2 +  vy_range(range_count).^2).*dt;
    vy_range(range_count + 1) =  vy_range(range_count) + dt.*(-g) + (-B./m).*(Py(range_count)./p_0).* vy_range(range_count).*sqrt(vx_range(range_count).^2 + vy_range(range_count).^2).*dt;
    
    %position
   x_range(range_count + 1) = x_range(range_count) + dt.*(vx_range(range_count));
   y_range(range_count + 1) = y_range(range_count) + dt.*( vy_range(range_count));
    
    %increment range_count
    range_count = range_count + 1;
end%-----while loop end--------------------with air drag and air density-------------
 loop_count=loop_count+1;
 
r_range_y_compare = - (y_range(end - 1) / y_range(end));
range_end= (x_range(end - 1) + r_range_y_compare*x_range(end)) / (r_range_y_compare + 1);
x_range(end)=range_end;
range_compare(k)=range_end;
 
if k==1    
x_range_converage1=x_range;
y_range_converage1=y_range;
elseif k==2
x_range_converage2=x_range;
y_range_converage2=y_range;
end
 
 end%----------for k loop end---------------------------------------
 
    if range_compare(1)>x_range_maxium(end) || range_compare(2)>x_range_maxium(end)
        if range_compare(1)>range_compare(2)
        theta_max2=theta_max2-(((converage_converage_ratio)^(n))*(pi/4));
            y_range_maxium=y_range_converage1;
            x_range_minum=zeros();
            x_range_minum=x_range_maxium;
            x_range_maxium=x_range_converage1;
              n=n+1; 
        elseif range_compare(1)<range_compare(2)
            theta_max2=theta_max2+(((converage_converage_ratio)^(n))*(pi/4));
            y_range_maxium=y_range_converage2;
            x_range_minum=zeros();
            x_range_minum=x_range_maxium;
            x_range_maxium=x_range_converage2;
            n=n+1; 
        elseif range_compare(1)== range_compare(2)
            if h==0
            h=1;
            end
        fprintf('caustion: range_compare(1)==range_compare(2)= %f,theta_max2=%f\n',range_compare(1),theta_max2);
         n=n+1;
        end
        
    else
        if range_compare(1)>range_compare(2)
            x_range_minum=x_range_converage1;
 
        else
            x_range_minum=x_range_converage2;
           
        end
    n=n+1;
    if n>200;
        disp('adjust the tolorance regin');
        break;
    end
    
    end
 
end%while loop end----------------while abs max -min
n=1;
end%for loop end------------------

fprintf('\nFinduing the maxmium range used %d loops.\nTime for finding the maxmium was shown below\n',loop_count);
toc;
 
%----------------------------plot------------------------------------------
figure
plot(x_normal,y_normal,'r-o', x_airdrag, y_airdrag,'m-x',x_density,y_density, 'y.-',x_range_maxium, y_range_maxium, 'b.-');
 
 
legend('Without Drag or density', 'With Drag only', 'With Drag and Density','max with air resistance and change of density');
xlim([0 inf]);
ylim([0 inf]);
xlabel('X (m)');
ylabel('Y (m)');
title('Projectile Motion');
  
%---------------------------display other info----------------------------- 
fprintf('\n The angle that maxmizes the range of the projectile that\n accounted for air drag and density is %.5f degree\n with a calculation tolorance of %f m.\n the final distance is %.3f m\n\n', theta_max2/pi*180,range_approx_tolorance,x_range_maxium(end));
fprintf('\nThe projectile that did not take into account drag or density \n traveled approximately %.3f meters.\n', range_normal); 
fprintf('The projectile that only took into account drag traveled approx\nimately %.3f meters.\n', range_airdrag); 
fprintf('The projectile that took into account both drag and density trav\neled approximately %.3f meters.\n', range_density);  
