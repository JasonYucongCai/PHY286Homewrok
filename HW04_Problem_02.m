function [Yt] = HW04_Problem_02(x0, sigma, c, dx, length_x, t_max)
%
%function [Yt] = HW04_Problem_02(x0, sigma, c, dx, length_x, t_max)
%
%This functions represents a string with one end that oscillates according
%to A*sin(omega * t).  The other end of the string is fixed.
%
%Input: x0: used only if one wants the string to have some initial position
%           (e.g. the string starts pulled up, centered at x0
%       sigma: used only if one wants the string to have some initial
%              position; sigma is the width of the gaussian pluck
%       c: the speed of the wave
%       dx: time step
%       length_x: the length of the string
%       t_max: the time over which the wave is calculated
%
%Output: Yt: a matrix of the position and time of the wave
%
%GOOD INPUT ------- HW04_Problem_02(.4, 0, 500, .01, 1, 1);

%initial cond.-------------------------------------------------------------
%--------------------------------------------------------------------------

r=1;
dt = r*dx/c;

%establish vectors for 'x' and 't'
x = linspace(0, length_x, round(length_x/dx))';
t = linspace(0, t_max, round(t_max/dt));

%establish initial height of string (e.g. if a part is pulled up at some
%point
Yt = zeros(length(x), length(t));
y0 = set_initial_shape(x, x0, sigma);

%updates the Yt array to account for the initial position of the wave
Yt(:,1) = y0;


%calculate-----------------------------------------------------------------
%--------------------------------------------------------------------------

for tn=2:length(t)
    
    %establish 'current' position
    yc = Yt(:,tn-1);
    yc = setbc(yc, tn);
    
    %establish 'previous' position
    if tn == 2
       yp = yc;
    else
       yp = Yt(:,tn-2);
    end
    
    %make sure boundary conditions are still satisfied
    yp = setbc(yp, tn);
    
    %propagate the wave (e.g. calculate the 'new' position)
    yn = propagate_wave(yc, yp, r);
    
    %make sure boundary conditions are still satisfied
    yn = setbc(yn, tn);
   
    %update the vector
    Yt(:,tn) = yn;
   
    %plot
   if nargout == 0
       plot(x, yp, 'c', x, yc, 'k', x, yn, 'm');
       legend('Prev', 'Current', 'Next');
       ylim((max(y0) + 1).*[-3 3]);
       title(sprintf('Time step %d of %d', tn, length(t)));
       drawnow;
   end
   
end

% ===========================================================
%                       Propogate Wave
%============================================================
function yn = propagate_wave(yc, yp, r)

%establish the size of the 'new' position
yn = zeros(size(yc));

%solve the wave equation
for xn=2:length(yc)-1
    yn(xn) = 2*(1-r.*r).*yc(xn) - yp(xn) + r.*r.*(yc(xn+1)+yc(xn-1)); 
end

% ===========================================================
%                       Initial Shape
%============================================================
function y = set_initial_shape(x, x0, sigma)

if sigma == 0
    y = 0;
else
    y = exp(-(x-x0).^2/sigma);
end

% ===========================================================
%                       Boundary Conditions
%============================================================
function yout = setbc(yin, t)

A = .5;
omega = 1;

yout = yin;
yout(1) = A.*sin(omega.*(t-1)./10);
yout(end) = 0;