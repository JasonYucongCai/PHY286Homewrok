function [] = HW01_Problem_08(T, Vmin, Vmax)
%function [] = HW01_Problem_08(T, Vmin, Vmax)
%
%This function will plot the wind chill factor based on a given temperature
%and wind speed.
%
%Input: T: the temperature (F)
%Input: Vmin: the minimum wind speed (mph)
%Input: Vmin: the maximum wind speed (mph)
%
%Output: n/a: the output for this function will be a plot

%--------------------error check------------------------------

if nargin ~= 3
    error('You must enter three arguments.  They must be numeric.');
end

if Vmin >= Vmax
    error('The minimum wind speed must be LESS THAN the maximum wind speed.');
end

if Vmin < 0
   error('The minimum wind speed must be a positive value.');  
end

%--------------------set up----------------------------------

wind_speeds = linspace(Vmin, Vmax, 100); % generate a vector of 100 wind speeds evenly spaced between the min and max value input 

for i = 1:100
    WCF(i) = 35.7 + .6.*T - (35.7.*(wind_speeds(i).^.16)) + .43.*T.*(wind_speeds(i).^.16);% calculates the wind chill factor for each wind speed
end

%--------------------plot------------------------------------

plot(wind_speeds, WCF, 'r-');
legend('Wind speed vs. WCF');
xlabel('Wind Speed (mph)');
ylabel('Wind Chill Factor (F)');
title('Wind speed vs. WCF');
