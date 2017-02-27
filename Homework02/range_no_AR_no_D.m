function [range, x, y, vx, vy] = range_no_AR_no_D(theta, v0, B, m, dt)
%
%function [range, x, y, vx, vy] = range_no_AR_no_D(theta, v0, B, m, dt)
%
%This function calculates the position, velocity, and range of a projectile
%that does not take into account air resistance or density.  It uses inputs
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

%--------------------------compute-----------------------------------------

%projectile with neither drag nor density
counter = 1;
while y(counter) >= 0
   
    %velocity
    vx(counter + 1) = vx(counter) + dt.*0;
    vy(counter + 1) = vy(counter) + dt.*(-g);
    
    %position
    x(counter + 1) = x(counter) + dt.*(vx(counter));
    y(counter + 1) = y(counter) + dt.*(vy(counter));
   
    %increment counter_n
    counter = counter + 1;
end

%range
r = - (y(end - 1) / y(end));
range = (x(end - 1) + r*x(end)) / (r + 1);
