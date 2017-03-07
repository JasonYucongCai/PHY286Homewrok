function [theta, omega, t, energy, energy2] = simple_pendulum(theta0, omega0, L, dt, tmax)
% 
% function [theta, omega, t] = simple_pendulum(theta0, omega0, L)
% Calculate the time varying theta and omega (angular disp. and velocity) for a 
% simple pendulum in gravity
% Inputs: 
% 	theta0: initial angle in radians
% 	omega0: initial angular vel. in radians/s
% 	L: length of pendulum (in m)
% 	dt: time step (s)
% 	tmax: max time to compute (s)
% 
% Outputs: 
% 	theta: angular displacement in time
% 	omega: angular velocity in time
% 	t: time 
% 

t = [0:dt:tmax]; %create a vector time
g = 9.8; %assume we're on earth
mass = 1;

theta = zeros(size(t));
omega = zeros(size(t));
energy = zeros(size(t));
energy2 = zeros(size(t));

theta(1) = theta0;
omega(1) = omega0;
energy(1) = .5.*mass.*L.^2.*(omega(1).^2 + (g.*theta(1).^2)./L);
energy2(1) = .5.*mass.*L.^2.*(omega(1).^2 + (g.*theta(1).^2)./L);

for n1=2:length(t)
	omega(n1) = omega(n1-1) - g.*theta(n1-1).*dt./L;
	theta(n1) = theta(n1-1) + omega(n1).*dt;
    
    energy(n1) = .5.*mass.*L.^2.*(omega(n1).^2 + (g.*(theta(n1).^2))./L);
end

for n1=2:length(t)
	omega(n1) = omega(n1-1) - g.*theta(n1-1).*dt./L;
	theta(n1) = theta(n1-1) + omega(n1 -1 ).*dt;
    
    energy2(n1) = .5.*mass.*L.^2.*(omega(n1).^2 + (g.*(theta(n1).^2))./L);
end

figure(1); clf; %plot in figure window 1 and clear it 
subplot(3,1,1);
plot(t, theta, '-k'); xlabel('Time, t[s]'); ylabel('\theta [Rad]'); 
grid on
subplot(3,1,2);
plot(t, omega, '--k'); xlabel('Time, t[s]'); ylabel('\omega [Rad/s]');
grid on
subplot(3,1,3);
plot(t, energy, '--k', t, energy2, 'r'); xlabel('Time, t[s]'); ylabel('energy [J]');
grid on
