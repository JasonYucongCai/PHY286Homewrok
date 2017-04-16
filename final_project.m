function [n, m, h, I_k, I_na, I_leak, I_injected, V, time] = final_project(n0, m0, h0, I_excite, V0, dt, tmax)
%
%Note book gives n0 = 0, m0 = 0, h0 = 1 for t = 0
%based on other figures, V0 = 0 and I_excite could be a range of values(I
%use .01 based on figure 12.40)

%==========================================================================
%                           Initial Setup
%==========================================================================

%establish time vector
time = linspace(0, tmax, round(tmax/dt));

%establish n, m, h, I, V, alpha, beta vectors
n = zeros(size(time));
n(1) = n0;
alpha_n = zeros(size(time));
beta_n = zeros(size(time));

m = zeros(size(time));
m(1) = m0;
alpha_m = zeros(size(time));
beta_m = zeros(size(time));

h = zeros(size(time));
h(1) = h0;
alpha_h = zeros(size(time));
beta_h = zeros(size(time));

I_k = zeros(size(time));
I_na = zeros(size(time));
I_leak = zeros(size(time));
I_injected = zeros(size(time));

V = zeros(size(time));
V(1) = V0;

%constants
C_m = 1;        %uF / cm^2
G_k = .036;     %S / cm^2
V_k = -12;      %mV
G_na = .12;     %S / cm^2
V_na = 115;     %mV
G_leak = .3;    %mS / cm^2
V_leak = 10.6;  %mV


%==========================================================================
%                             Calculate
%==========================================================================

for count = 1:length(time) - 1
    
    %EULER METHOD
    
    %n related vectors
    alpha_n(count) = (.01.*(10-V(count))) ./ (-1 + (exp((1/10).*(10-V(count)))));
    beta_n(count) = .125.*exp(-V(count)./20);
    n(count + 1) = n(count) + dt.*(alpha_n(count).*(1-n(count)) - beta_n(count).*n(count));
    
    %m related vectors
    alpha_m(count) = (.1.*(25-V(count))) ./ (-1 + (exp((1/10).*(25-V(count)))));
    beta_m(count) = 4.*exp(-V(count)./18);
    m(count + 1) = m(count) + dt.*(alpha_m(count).*(1-m(count)) - beta_m(count).*m(count));   
    
    %h related vectors
    alpha_h(count) = .07.*exp(-V(count)./20);
    beta_h(count) = 1 ./ (1 + exp((1/10).*(30-V(count))));
    h(count + 1) = h(count) + dt.*(alpha_h(count).*(1-h(count)) - beta_h(count).*h(count));
    
    %I related vectors
    I_k(count) = G_k.*(n(count).^4).*(V(count) - V_k);
    I_na(count) = G_na.*(m(count).^3).*h(count).*(V(count)-V_na);
    I_leak(count) = G_leak.*(V(count) - V_leak);
    
    %injected current, only for a specific range of time
    if count >= min(find(time >= .2*tmax)) && count <= max(find(time <= .3*tmax))
    I_injected(count) = I_excite;
    end
    
    %V vector
    V(count + 1) = V(count) + dt.*((1/C_m).*(I_injected(count) - I_k(count) - I_na(count) - I_leak(count)));

end


%==========================================================================
%                             Plot
%==========================================================================

plot(time, V);









