function [ output_args ] = HW03_Problem_022( x_initial, y_initial,z_initial,delta,t_max,delta_r)
%cloned form
%   Detailed explanation goes here

%a=10;
%b=8/3;

t_axis=[1:delta:t_max];
N=length(t_axis);


x(1)=x_initial;
y(1)=y_initial;
z(1)=z_initial;

plot_r = [];
plot_z = [];

for r_input=25:delta_r:160

for i=1:N-1
%----------------------------------RK mediate------------------------------
[fxm,~,~]=RK_function(x(i),y(i),z(i),r_input);
[~,fym,~]=RK_function(x(i),y(i),z(i),r_input);
[~,~,fzm]=RK_function(x(i),y(i),z(i),r_input);

dxm=x(i)+1/2*fxm*delta;
dym=y(i)+1/2*fym*delta;
dzm=z(i)+1/2*fzm*delta;


%----------------------------------RK mediate------------------------------
    
[fxrk,~,~]=RK_function(dxm,dym,dzm,r_input);
[~,fyrk,~]=RK_function(dxm,dym,dzm,r_input);
[~,~,fzrk]=RK_function(dxm,dym,dzm,r_input);

x(i+1)=x(i)+fxrk*delta;
y(i+1)=y(i)+fyrk*delta;
z(i+1)=z(i)+fzrk*delta;

end


%----------------------------------------------------------cloned code -------------------------------------------
index = find(t_axis>t_max/1.5);

indexed_x = x(index);
indexed_y = y(index); %**
indexed_z = z(index); %**
indexed_t = t_axis(index);
t_vector = zeros(length(indexed_t), 1);
z_vector = zeros(length(indexed_z), 1);

sign_x_vector = sign(indexed_x);
diff_x_vector = [diff(sign_x_vector) 0];




for count = 1:length(indexed_x)

    %determines when x changes sign
    if abs(diff_x_vector(count)) == 2
        
        %stores the value of t, z into the plot_t, plot_z vectors
        t_vector(count) = indexed_t(count);
        z_vector(count) = indexed_z(count);
    end
    
end

%[time, output] = getTime(x, y, z, t, tmax);
%function [t_vector, z_vector] = getTime(x, y, z, t, tmax)
%get rid of the extra 0's in the vector

index_z = find(z_vector > 0); 
%use only unique values of plot_t and plot_z
plot_r = [plot_r; ones(size(unique(z_vector(index_z)))).*r_input];
plot_z = [plot_z; unique(z_vector(index_z))];





end

figure (1)
plot(plot_r, plot_z, '.k');
xlabel('r');
ylabel('z');
title('Bifurcation diagram of the Lorenz Model');
%----------------------------------------------------------cloned code -------------------------------------------

function [dx,dy,dz] = RK_function(x,y,z,r) 
% evaluate the function at x and return its value 
a=10;
b=8/3;
dx=a*(y-x);
dy=-x*z+r*x-y;
dz=x*y-b*z;
end




end
