function [Nt, Nt_exact, time] = HW02_Problem_01(N_0, a, b, dt, t_end)
%function [Nt, Nt_exact, time] = HW02_Problem_01(N_0, a, b, dt, t_end)
%
%This function solves the equation N'(t) = aN(t) - bN^2(t) where N(t)
%represents a given population over some period of time. a and b are both
%constants
%
%Input: N_0: the initial population at time = 0
%       a: a constant
%       b: a constant
%       dt: time step (for this problem dt greatly affects the plot,
%           it should be quite small for most input values ~.001 or less)
%       t_end: the time span over which the population will be
%              calculated(again t_end greatly affects the plot, it should be around 1 for most input values)
%
%Output: Nt: the approximate solution to the differential equation given
%            above
%        Nt_exact: the exact solution to the differential equation given
%                  above
%        time: each time corresponds to a population value
%
%NOTES: for small Nt, b ~ 3, a ~ 10
%       for large Nt (around 1000), b ~ .01, a ~ 10
%       
%       when the ratio a/b = N_0 a horizontal line will result at y = N_0

%----------------------------error check-----------------------------------

if nargin ~= 5
    error('This function requires 5 inputs');
end

%-----------------initial values----------------

time = [0:dt:t_end]; 
Nt = zeros(size(time));
Nt(1) = N_0;

% solve for exact analytical solution
Nt_exact = (a.*N_0) ./ (b.*N_0 + (a - b.*N_0).*exp(-a.*time));

%-----------------perform calculation-----------

for count = 1:length(time) -1
    Nt(count + 1) = Nt(count) + dt.*(a.*Nt(count) - b.*Nt(count).^2);
end

%-----------------plot--------------------------

plot(time, Nt_exact, 'r*-', time, Nt, 'b.-');
xlabel('Time');
ylabel('Population N(t)');
legend('Exact Sol.', 'Approx. Sol.');
title({'Population vs. Time'; ['N_0 = ', num2str(N_0)]; ['a = ', num2str(a)]; ['b = ', num2str(b)]});



