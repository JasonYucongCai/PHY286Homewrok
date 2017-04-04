function [x, y, z, t] = HW3_problem2(x0y0z0, tmax, dt, r)
%
%function [x, y, z, t] = HW3_problem2(x0y0z0, tmax, dt, r)
%
%This function produces a bifurcation plot of the Lorenz model.
%It produces a plot of r vs. z at points only when x ~ 0.
%
%Inputs: x0y0z0: a vector containing the initial values of x, y, z
%        tmax: the maximum amount of time to calculate
%        dt: time step
%        r: a vector of the 'r' coefficient (25 -> 160)
%
%Outputs: x: vector of x values output by the Lorenz model
%         y: vector of y values output by the Lorenz model
%         z: vector of z values output by the Lorenz model
%         t: time vector
%
%NOTE -- good input values seem to be x0y0z0 = [1 0 0], tmax = 50, 
%                                     dt = .001, r = [25:.1:160]

%--------------------initial setup-----------------------------------------

%initial values
sigma = 10;
b = 8/3;

x(1) = x0y0z0(1);
y(1) = x0y0z0(2);
z(1) = x0y0z0(3);

t = 0:dt:tmax;

plot_r = [];
plot_z = [];

%--------------------calculate---------------------------------------------

for nn = 1:length(r)
    R = r(nn);

%calculate RK2
for count = 1:length(t) - 1
    
    %step 1 -- equations
    equation_1 = sigma.*(y(count)-x(count));
    equation_2 = -x(count).*z(count) + R.*x(count) - y(count);
    equation_3 = x(count).*y(count) - b.*z(count);
    
    t_p = t(count) + (1/2).*dt;

    %step 2 -- equations'
    x_prime = x(count) + (1/2).*dt.*equation_1;
    y_prime = y(count) + (1/2).*dt.*equation_2;
    z_prime = z(count) + (1/2).*dt.*equation_3;
    
    %step 3 -- RK equations
    RK_equation_1 = sigma.*(y_prime-x_prime);
    RK_equation_2 = -x_prime.*z_prime + R.*x_prime - y_prime;
    RK_equation_3 = x_prime.*y_prime - b.*z_prime;
    
    %step 4 -- calculate x, y, z
    x(count + 1) = x(count) + dt.*RK_equation_1;
    y(count + 1) = y(count) + dt.*RK_equation_2;
    z(count + 1) = z(count) + dt.*RK_equation_3;
    
end

%--------------------find z when x ~ 0-------------------------------------

%find values of z at time t to plot
[time, output] = getTime(x, y, z, t, tmax);

%get rid of the extra 0's in the vector
index_z = find(output > 0); 

%use only unique values of plot_t and plot_z
plot_r = [plot_r; ones(size(unique(output(index_z)))).*R];
plot_z = [plot_z; unique(output(index_z))];

end

%--------------------plot and display--------------------------------------

%plot
plot(plot_r, plot_z, '.k');
xlabel('r');
ylabel('z');
title('Bifurcation diagram of the Lorenz Model');


%==========================================================================
%function that returns 'in phase' values for t and z when x = 0
%==========================================================================
function [t_vector, z_vector] = getTime(x, y, z, t, tmax)

%use times after initial transient has decayed
index = find(t>tmax/1.5);

%create new x, y, z vectors without transients (POSSIBLE ERROR)
indexed_x = x(index);
indexed_y = y(index); %**
indexed_z = z(index); %**
indexed_t = t(index);

%create a time vector of length indexed_t to be later filled with some
%values of t that should be plotted
t_vector = zeros(length(indexed_t), 1);
z_vector = zeros(length(indexed_z), 1);

%create a vector of the 'sign' of x
sign_x_vector = sign(indexed_x);

%create a vector of the diff of x and add one spot at the end to keep the
%length consistent
diff_x_vector = [diff(sign_x_vector) 0];

for count = 1:length(indexed_x)

    %determines when x changes sign
    if abs(diff_x_vector(count)) == 2
        
        %stores the value of t, z into the plot_t, plot_z vectors
        t_vector(count) = indexed_t(count);
        z_vector(count) = indexed_z(count);
    end
    
end


