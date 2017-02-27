function [range, x, y, vx, vy] = range_AR_and_D(theta, v0, B, m, dt)
%
%function [range, x, y, vx, vy] = range_AR_and_D(theta, v0, B, m, dt)
%
%This function calculates the position, velocity, and range of a projectile
%that takes into account both air resistance and density.  It uses inputs
%descriped in HW02_Problem_02.

%--------------------------intialize-----------------------------------------

%position
x(1) = 0;
y(1) = 0;

%velocity
vx(1) = v0.*cosd(theta);
vy(1) = v0.*sind(theta);

%constants
g = 9.8;
y_0 = 10000; %approximate value from book, calculated as (kB*T/m*g)

%pressure
p_0 = 1.225; %approximate value of air density at sea level

Py(1) = p_0.*exp(-y(1)./y_0);

%--------------------------compute-----------------------------------------

%projectile with both drag and density
counter_dd = 1;
while y(counter_dd) >= 0
    %solve for density
    Py(counter_dd) = p_0.*exp(-y(counter_dd)./y_0);
    
    %velocity
    vx(counter_dd + 1) = vx(counter_dd) + dt.*0 + (-B./m).*(Py(counter_dd)./p_0).*vx(counter_dd).*sqrt(vx(counter_dd).*vx(counter_dd) + vy(counter_dd).*vy(counter_dd)).*dt;
    vy(counter_dd + 1) = vy(counter_dd) + dt.*(-g) + (-B./m).*(Py(counter_dd)./p_0).*vy(counter_dd).*sqrt(vx(counter_dd).*vx(counter_dd) + vy(counter_dd).*vy(counter_dd)).*dt;
    
    %position
    x(counter_dd + 1) = x(counter_dd) + dt.*(vx(counter_dd));
    y(counter_dd + 1) = y(counter_dd) + dt.*(vy(counter_dd));
    
    %increment counter_dd
    counter_dd = counter_dd + 1;
end

%range
r = - (y(end - 1) / y(end));
range = (x(end - 1) + r*x(end)) / (r + 1);
















