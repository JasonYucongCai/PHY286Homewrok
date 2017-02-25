function [x_normal, y_normal, x_airdrag, y_airdrag, x_density, y_density, Py] = HW02_Problem_02(v0, theta, dt, B, m)
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

%position
x_normal(1) = 0;
y_normal(1) = 0;

x_airdrag(1) = 0;
y_airdrag(1) = 0;

x_density(1) = 0;
y_density(1) = 0;

%velocity
vx_normal(1) = v0.*cosd(theta);
vy_normal(1) = v0.*sind(theta);

vx_airdrag(1) = v0.*cosd(theta);
vy_airdrag(1) = v0.*sind(theta);

vx_density(1) = v0.*cosd(theta);
vy_density(1) = v0.*sind(theta);

%constants
g = 9.8;
y_0 = 10000; %approximate value from book, calculated as (kB*T/m*g)

%pressure
p_0 = 1.225; %approximate value of air density at sea level

Py(1) = p_0.*exp(-y_density(1)./y_0);

%--------------------------compute-----------------------------------------

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
end

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
end

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
end

%compute range of projectile
r_normal = - (y_normal(end - 1) / y_normal(end));
range_normal = (x_normal(end - 1) + r_normal*x_normal(end)) / (r_normal + 1);

r_airdrag = - (y_airdrag(end - 1) / y_airdrag(end));
range_airdrag = (x_airdrag(end - 1) + r_airdrag*x_airdrag(end)) / (r_airdrag + 1);

r_density = - (y_density(end - 1) / y_density(end));
range_density = (x_density(end - 1) + r_density*x_density(end)) / (r_density + 1);

%----------------------------plot------------------------------------------

plot(x_normal,y_normal,'r.-', x_airdrag, y_airdrag,'b.-', x_density, y_density, 'g.-');
legend('Without Drag or density', 'With Drag only', 'With Drag and Density');
xlim([0 inf]);
ylim([0 inf]);
xlabel('X (m)');
ylabel('Y (m)');
title('Projectile Motion');
  
%---------------------------display other info----------------------------- 

disp(' ');
fprintf('The projectile that did not take into account drag or density traveled approximately %.3f meters.', range_normal); 
disp(' ');
fprintf('The projectile that only took into account drag traveled approximately %.3f meters.', range_airdrag); 
disp(' ');
fprintf('The projectile that took into account both drag and density traveled approximately %.3f meters.', range_density);  
 