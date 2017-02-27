function [range_density] = HW02_Problem_02(theta, v0, B, m, dt, initialguess, increment, tol)
%
%function [range_density] = HW02_Problem_02(theta, v0, B, m, dt, initialguess, increment, tol)
%
%This function makes use of 3 other functions that calculate the range,
%position, and velocity of a projectile with given inputs.  They are range_no_AR_no_D.m, range_only_AR.m, and range_AR_and_D.m.  With this
%information, this function plots the projectile motion of the 3 different
%projectiles, plots the range of the third projectile (takes into account
%density and air resistance) vs. theta (from 1 to 90 degrees), and
%calculates the angle that will yield the maximum range for the third
%projectile.  This last step is done using an optimization method found in
%appendix B of the textbook.
%
%SI units are assumed and theta is in degrees.
%
%Inputs: theta: angle of fire of the projectile
%        v0: initial velocity of the projectile
%        B: drag coefficient
%        m: mass
%        dt: time step
%        intialguess: angle to start with when finding the angle that will
%                     maximize the range of the third projectile (a good
%                     choice is 45)
%        increment: value that determines a bracket that holds a maximum
%                   (a good choice is .1)
%        tol: a tolerance zone for deciding how accurate the angle
%             calculated should be
%
%Outputs: range_density: the maximum range of the third projectile
          
%----------------------------error check-----------------------------------
if nargin ~= 8
    error('This function requires 8 inputs');
end

if theta <= 0 || v0 <= 0 || B < 0 || m < 0 || dt <= 0 || initialguess < 1 || initialguess > 90 || increment <= 0 || tol <= 0 
    error('There was an error.  Recall that you must choose an angle between 0 and 90 degrees, an object cannot have negative mass or velocity, dt should be greater than 0, the intial guess for checking the max angle should be between 0 and 90, the increment chosen should be greater than 0, and the tolerance zone should be greater than 0.');
end

%----------------------------calculate-------------------------------------

[range1,x1_vec,y1_vec,vx1_vec,vy1_vec] = range_no_AR_no_D(theta,v0, B, m, dt);
[range2,x2_vec,y2_vec,vx2_vec,vy2_vec] = range_only_AR(theta,v0, B, m, dt);
[range3,x3_vec,y3_vec,vx3_vec,vy3_vec] = range_AR_and_D(theta,v0, B, m, dt);

Range1 = range1;
Range2 = range2;
Range3 = range3;

%----------------------------plot------------------------------------------

%plot projectile motion of 3 cases
% 1) no drag or density taken into account
% 2) only drag taken into account
% 3) both drag and density taken into account
subplot(2,1,1); plot(x1_vec, y1_vec,'r.-', x2_vec, y2_vec, 'b.-', x3_vec, y3_vec, 'g.-');
legend('1) Without Drag or density', '2) With Drag only', '3) With Drag and Density');
xlim([0 inf]);
ylim([0 inf]);
xlabel('X (m)');
ylabel('Y (m)');
title('Projectile Motion');

%calculates the range of 3) (from above) for each theta (1:90)
theta = 1:90;
for i=1:length(theta)
   [r,x3_vec,y3_vec,vx3_vec,vy3_vec] = range_AR_and_D(theta(i), v0, B, m, dt);
   range3(i) = r;
end

subplot(2,1,2); plot(theta,range3);
legend('Possible ranges of 3)');
xlim([0 inf]);
ylim([0 inf]);
xlabel('Theta (degrees)');
ylabel('Range (m)');
title('Range vs. Theta');

%----------------------------find max angle for projectile with drag and density------------------------------------------

dx = increment; %increment
x(1) = initialguess; %intiail guess 
x(2) = x(1) + dx; %x value that we change and input and stuff

if range_AR_and_D(x(2), v0, B, m, dt) <= range_AR_and_D(x(1), v0, B, m, dt)  %interchange x(1) and x(2)
    temp = x(2);
    x(2) = x(1);
    x(1) = temp;
    dx = -dx;
end

count = 3;
while abs(2*dx) >= tol
    x(count) = x(count - 1) + dx; %increment x
    if range_AR_and_D(x(count-1), v0, B, m, dt) > range_AR_and_D(x(count-2),  v0, B, m, dt) && (range_AR_and_D(x(count-1), v0, B, m, dt) > range_AR_and_D(x(count), v0, B, m, dt)) %we found a max and it is bracketed
        if abs(2*dx) < tol
           range_density = x(count);
           
        else
           dx = -dx/2;
        end
        range_density = x(count);
    end
         
    if count >= 1000
        error('No maximum could be found!');
    end
    
    count = count + 1;
end

%----------------------------display------------------------------------------
disp(' ');
fprintf('The range of 1)is %.3f m', Range1);
disp(' ');
fprintf('The range of 2)is %.3f m', Range2);
disp(' ');
fprintf('The range of 3)is %.3f m', Range3);
disp(' ');
disp(' ');
fprintf('The angle that maximizes the range of 3) based on the given inputs is %.3f degrees', range_density);
disp(' ');