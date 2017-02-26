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
r_normal = - (y_normal(end - 1) / y_normal(end));
range_normal = (x_normal(end - 1) + r_normal * x_normal(end)) / (r_normal + 1);

r_airdrag = - (y_airdrag(end - 1) / y_airdrag(end));
range_airdrag = (x_airdrag(end - 1) + r_airdrag*x_airdrag(end)) / (r_airdrag + 1);

r_density = - (y_density(end - 1) / y_density(end));
range_density = (x_density(end - 1) + r_density*x_density(end)) / (r_density + 1);

%-------------------------------------------------------------------------
%comput the maxmium angle for the normal------------------------up
%chose the reference angle
%known v0, theta, dt, B, m,range_approx_tolorance
%converage the variables
%{
theta2=theta+ (pi)/3;
theta_max2=theta;
if theta2>=pi/2 %creat a initial comparation
theta2=theta-(pi/2);
end

x_normal(counter_n)=range_normal;%fix the range msitakes
xmax=x_normal;
xmin=zeros();
yax2=y_normal;
n=200;

%-----------------------------------------------------------------------
k=1;
rangecom(1)=range_normal;
%}
%{

for i=1:n-1
    theta2=pi/2/n*i
    
    counter_n = 1;
    %cleanup;
    vx_normal2=zeros();
    vy_normal2=zeros();
    x_normal2=zeros();
    y_normal2=zeros();
    vx_normal2(1) = v0.*cos(theta2);
    vy_normal2(1) = v0.*sin(theta2);
    
    while y_normal2(counter_n) >= 0
    %velocity
    vx_normal2(counter_n + 1) = vx_normal2(counter_n) + dt.*0;
    vy_normal2(counter_n + 1) = vy_normal2(counter_n) + dt.*(-g);
    %position
    x_normal2(counter_n + 1) = x_normal2(counter_n) + dt.*(vx_normal2(counter_n));
    y_normal2(counter_n + 1) = y_normal2(counter_n) + dt.*(vy_normal2(counter_n));
    %increment counter_n
    counter_n = counter_n + 1;
    end%-------------------without out air drag
    

    r_normal2 = - (y_normal2(end - 1) / y_normal2(end));
    range2 = (x_normal2(end - 1) + r_normal2 * x_normal2(end)) / (r_normal2 + 1);
    
    if range2>xmax(end)
        x_normal2(end)=range2;
        yax2=y_normal2;
        xmin=xmax;
        xmax=x_normal2;
    end
       
    end

%}

theta2=pi/4;
    counter_dd = 1;
    %cleanup;
    vx_normal2=zeros();
    vy_normal2=zeros();
    x_normal2=zeros();
    y_normal2=zeros();
    vx_normal2(1) = v0.*cos(theta2);
    vy_normal2(1) = v0.*sin(theta2);
    
while vy_normal2(counter_dd) >= 0
    %solve for density
    Py(counter_dd) = p_0.*exp(-y_density(counter_dd)./y_0);
    
    %velocity
    vx_normal2(counter_dd + 1) = vx_normal2(counter_dd) + dt.*0 + (-B./m).*(Py(counter_dd)./p_0).*vx_normal2(counter_dd).*sqrt(vx_normal2(counter_dd).^2 +  vy_normal2(counter_dd).^2).*dt;
    vy_normal2(counter_dd + 1) =  vy_normal2(counter_dd) + dt.*(-g) + (-B./m).*(Py(counter_dd)./p_0).* vy_normal2(counter_dd).*sqrt(vx_normal2(counter_dd).^2 + vy_normal2(counter_dd).^2).*dt;
    
    %position
   x_normal2(counter_dd + 1) = x_normal2(counter_dd) + dt.*(vx_normal2(counter_dd));
   y_normal2(counter_dd + 1) = y_normal2(counter_dd) + dt.*( vy_normal2(counter_dd));
    
    %increment counter_dd
    counter_dd = counter_dd + 1;
end%-------------------------with air drag and air density

r_normal2 = - (y_normal2(end - 1) / y_normal2(end));
range2= (x_normal2(end - 1) + r_normal2*x_normal2(end)) / (r_normal2 + 1);
x_normal2(end)=range2;

yax2=y_normal2;
xmax=x_normal2;
theta_max2=theta2;
xmin=zeros();
n=1;
n=2;
ratio=1/2;

while (abs(xmax(end)-xmin(end)))>=range_approx_tolorance  
    
    
    tempar1=zeros();
    tempar2=zeros();
    y21=zeros();
    y22=zeros();
    
    for k=1:1:2
        theta2=theta_max2+(((-1)^k)*((ratio)^(n))*(pi/2));        
    %cleanup;
    vx_normal2=zeros();
    vy_normal2=zeros();
    x_normal2=zeros();
    y_normal2=zeros();
    vx_normal2(1) = v0.*cos(theta2);
    vy_normal2(1) = v0.*sin(theta2);
    counter_dd = 1;
    %cleanup;
    while y_normal2(counter_dd) >= 0
    %solve for density
    Py(counter_dd) = p_0.*exp(-y_normal2(counter_dd)./y_0);
    
    %velocity
    vx_normal2(counter_dd + 1) = vx_normal2(counter_dd) + dt.*0 + (-B./m).*(Py(counter_dd)./p_0).*vx_normal2(counter_dd).*sqrt(vx_normal2(counter_dd).^2 +  vy_normal2(counter_dd).^2).*dt;
    vy_normal2(counter_dd + 1) =  vy_normal2(counter_dd) + dt.*(-g) + (-B./m).*(Py(counter_dd)./p_0).* vy_normal2(counter_dd).*sqrt(vx_normal2(counter_dd).^2 + vy_normal2(counter_dd).^2).*dt;
    
    %position
   x_normal2(counter_dd + 1) = x_normal2(counter_dd) + dt.*(vx_normal2(counter_dd));
   y_normal2(counter_dd + 1) = y_normal2(counter_dd) + dt.*( vy_normal2(counter_dd));
    
    %increment counter_dd
    counter_dd = counter_dd + 1;
end%-------------------------with air drag and air density

r_normal = - (y_normal2(end - 1) / y_normal2(end));
range2= (x_normal2(end - 1) + r_normal*x_normal2(end)) / (r_normal + 1);
x_normal2(end)=range2;
rangecom(k)=range2;

if k==1    
tempar1=x_normal2;
y21=y_normal2;
elseif k==2
tempar2=x_normal2;
y22=y_normal2;
end

end
    
    
    
    if rangecom(1)>xmax(end) || rangecom(2)>xmax(end)
        if rangecom(1)>rangecom(2)
        theta_max2=theta_max2-(((ratio)^(n))*(pi/2));
            %yax2=zeros();
            yax2=y21;
            xmin=zeros();
            xmin=xmax;
            xmax=tempar1;
              n=n+1; 
        elseif rangecom(1)<rangecom(2)
            theta_max2=theta_max2+(((ratio)^(n))*(pi/2));
            %yax2=zeros();
            yax2=y22;
            xmin=zeros();
            xmin=xmax;
            xmax=tempar2;
            n=n+1; 
        elseif rangecom(1)== rangecom(2)
            error(' rangecom(1)==rangecom(2)');
        end
        
        
    else
        if rangecom(1)>rangecom(2)
            xmin=tempar1;

        else
            xmin=tempar2;
           
        end
        
    n=n+1;
    if n>200;
        disp('adjust the tolorance regin');
        break;
    end
    
    end
 %fprintf('value xmax %f xmin %f theta2 %f theta_max2 %f\n',xmax(end),xmin(end),theta2,theta_max2);
    
 %{   
%change angle-----------------------------------------------   

k=k+1;
rangecom(k)=range2;
if  abs(rangecom(k)-rangecom(k-1))<=range_approx_tolorance
    theta2=(theta2+theta_max2)/2;
end
if  (k>=30) && abs(rangecom(k)-rangecom(k-20))<=range_approx_tolorance
    theta2=(theta2+pi/2)/2;
end


   if abs(rangecom(k)-rangecom(k-1))>range_approx_tolorance
    if range2>xmax(end)%-----------------------------------------------       
        
        switch 1
        case theta_max2<theta2
            if (theta2*(range2/xmax(end)))>=(pi/2)
                theta_max2=theta2;
                theta2=(pi/2+theta)/2;
            else
                theta_max2=theta2;
                theta2=theta2*(range2/xmax(end));
            end
        case theta_max2>theta2
            theta2=theta2*(xmax(end)/range2);
            otherwise
    
            disp('angles are equal1');
            break;    
        end
    elseif range2==xmax(end)%----------------------------------------
        break;
    else%range2<max(end)---------------------------------------------
       %theta_max2 didn't change  
       theta2=(theta2+theta_max2)/2;
        
    end
   end
%change range----------------------------------------------------
    

%}
    %{
    k=k+1;
rangecom(k)=range2;

%{
if  (k>5)&&(xmax(end)-xmin(end))<=range_approx_tolorance
    if mod(k,5)==1
    theta2=(theta2+theta_max2)/2;
    elseif mod(k,16)==1
%    theta2=theta_max2-abs(theta_max2+theta2)/2;
 theta2=theta2/2;
    end
    
end
%}

if range2==xmax(end)&& (theta2~=theta_max2);
theta2=(theta2+theta_max2)/2;
end


 if range2>xmax(end)&&theta_max2>theta2;
 theta_max2=theta2;
 theta2=theta2*(xmax(end)/range2);
 
 elseif range2>xmax(end)&&theta_max2<theta2
 theta_max2=theta2;
 theta2=(pi/2+theta)/2;
 end
 
 if range2<xmax(end)&&theta_max2<theta2
    theta2=theta_max2*(range2/xmax(end));     
 elseif range2<xmax(end)&& theta_max2>theta2
  theta2=(pi/2+theta2)/2;
 end
 

disp(theta2)
%}
 
 
 %   
    %{
    
    x_normal2(end)=range2;
        
    if range2>xmax(end)
        yax2=zeros();
        yax2=y_normal2;
        xmin=zeros();
        xmin=xmax;
        xmax=x_normal2;
    elseif range2==xmax(end)
        disp('range didn not change');
        break;
    else%range2<max(end)
        xmin=x_normal2;
       %max didn't change        
      % xmin=x_nromal;
    end
    fprintf('round %d',k);
    fprintf('value xmax %f xmin %f theta2 %f theta_max2 %f\n',xmax(end),xmin(end),theta2,theta_max2);
   %}
   
   %change angle and range
%-----------------------------------------------------------------
end


%comput the maxmium angle for the normal-----------------------bot
%-------------------------------------------    

%----------------------------plot------------------------------------------
figure

%plot(x_normal,y_normal,'r.-', x_airdrag, y_airdrag,'b.-');
%legend('Without Drag or density', 'With Drag only');
%hold on;
%plot(x_density, y_density, 'g.-', )
plot(x_normal,y_normal,'r-o', x_airdrag, y_airdrag,'m-x',x_density,y_density, 'y.-',xmax, yax2, 'b.-');


legend('Without Drag or density', 'With Drag only', 'With Drag and Density','max with air resistance');
xlim([0 inf]);
ylim([0 inf]);
xlabel('X (m)');
ylabel('Y (m)');
title('Projectile Motion');
  
%---------------------------display other info----------------------------- 
fprintf('\n The maxminum angle the ball could be throuwn under the\n condition of drag and change of flew density \nis %.5f with a tolorance of %f\nthe final distance is %.3f', theta_max2,range_approx_tolorance,xmax(end));
fprintf('\nThe projectile that did not take into account drag or density \n traveled approximately %.3f meters.\n', range_normal); 
fprintf('The projectile that only took into account drag traveled approx\nimately %.3f meters.\n', range_airdrag); 
fprintf('The projectile that took into account both drag and density trav\neled approximately %.3f meters.\n', range_density);  
