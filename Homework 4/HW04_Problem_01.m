function [Vxy, Ex] = HW04_Problem_01(V_in, V_out, length_inner, length_outer, num_grid_elements, TOL, alpha)
%
%function [Vxy, Ex] = HW04_Problem_01(V_in, V_out, length_inner, length_outer, num_grid_elements, TOL, alpha)
%
%This functions calculates the potential inside a box containing a
%parallel plate capacitor.  For the specific problem at hand, V_out = 0 and
%V_in = 1.  The potentials of the two plates in the box are V_in and -V_in.
%3 graphs are returned; the first is a 3-D graph of the potential, the
%second is a vector plot of the electric field generated, and the third is
%the graph of the equipotential surfaces.
%
%Input: V_in: the potential of the parellel plates inside the box, note
%             that the two plates will have the same magnitude of 
%             potential, but one will be negative and the other positive
%       V_out: the potential of the borders of the box, for this problem 0.
%       length_inner: the distance between the 2 plates of the capacitor
%       length_outer: the length of the box
%       num_grid_elements: how many 'regions' the box will be divided into
%       TOL: a value to determine when the 'convergence' has been reached
%       alpha: also determines when 'convergence' has been reached, should
%              be between 1 and 2 for best results
%
%Output: Vxy: a matrix containing the potential 
%        Ex: a matrix containing the electric field

%error---------------------------------------------------------------------
%--------------------------------------------------------------------------
if length_inner >= length_outer
    error('The outer length must be greater than the inner length!!');
end

%initial setup-------------------------------------------------------------
%--------------------------------------------------------------------------
x = linspace(-length_outer/2, length_outer/2, num_grid_elements);
y = linspace(-length_outer/2, length_outer/2, num_grid_elements);

x_index = find(x >= -length_inner/2 & x <= length_inner/2);
y_index = find(y >= -length_inner/2 & y <=length_inner/2);

%assume square shape
length_x = num_grid_elements; 
length_y = num_grid_elements;

%set V0 = 0 for all x and y
Vxy = zeros(num_grid_elements, num_grid_elements); 


%boundary conditions-------------------------------------------------------
%--------------------------------------------------------------------------
Vxy = set_boundary_cond(Vxy, V_in, V_out, length_x, x_index, y_index);

converged = 0;
iter = 1;

%calculate-----------------------------------------------------------------
%--------------------------------------------------------------------------
while converged == 0
    %storing 'old' Vxy to calculate the tolerance necessary for convergence
    Vxy_old = Vxy;
    
    %calculate
	for nx = 2:length_x-1
		for ny=2:length_y-1
			Vxy(ny,nx) = (Vxy(ny-1,nx) + Vxy(ny+1,nx) + Vxy(ny,nx+1) + Vxy(ny,nx-1) )./4; 
		end
    end
    
    %make sure that the boundary conditions are still satisfied
	Vxy = set_boundary_cond(Vxy, V_in, V_out, length_x, x_index, y_index);
    
    %calculate dV using SOR
	deltaV = mean(mean(abs(Vxy - Vxy_old)./alpha));

    %determine if convergence has happened
	if deltaV <= TOL
		converged = 1;
    end
    
    %increment number of iterations
	iter = iter + 1;
end

%calculate Electric Field
[Ex, Ey] = calc_electric_field(Vxy, x, y);

%set up graphing vectors
x_vector = linspace(-5,5,length(Vxy));
y_vector = linspace(-5,5,length(Vxy));


%display and plot----------------------------------------------------------
%--------------------------------------------------------------------------
disp(' ');
fprintf('Took %d steps to converge (alpha=%0.3f, Ngrid = %d, TOL = %f) \n', iter, alpha, num_grid_elements, TOL);
disp(' ');

if nargout == 0
	figure(1);
	contourf(x, y, Vxy, num_grid_elements); colorbar
	title('Equipotential for 2d parallel plate');
	xlabel('x'); ylabel('y');

	figure(2);
	quiver(x, y, Ex, Ey);
	set(gca, 'Fontsize', 15);
	xlabel('x'); ylabel('y');
	title('Corresponding E field');
    
    figure(3);
    mesh(x_vector, y_vector, Vxy);
    title('3-D Graph of the Potential');
    xlabel('x'); ylabel('y'); zlabel('Potential [V]');

end



% ===========================================================
%                       Boundary Conditions
%============================================================
function [Vout] = set_boundary_cond(Vxy, V_in, V_out, length_x, x_index, y_index)

%establishes the plates within the 'box' to have a V of V_in and -V_in
Vxy(1 + round(length_x / 4) :  length_x - round(length_x / 4), 1 + round(length_x / 4)) = V_in;
Vxy(1 + round(length_x / 4) :  length_x - round(length_x / 4), length_x - round(length_x / 4)) = -V_in;

%establishes the 'sides' of the box to be v_out
Vxy(1,:) = V_out; 
Vxy(length_x, :) = V_out; 
Vxy(:,1) = V_out; 
Vxy(:,length_x) = V_out; 

%returns the updated Vxy
Vout = Vxy; 

% ===========================================================
%                 Calculate Electric Field
%============================================================
function [Ex, Ey] = calc_electric_field(Vxy, x, y)

Ly = size(Vxy, 1);
Lx = size(Vxy, 2);
dx = mean(diff(x));
dy = mean(diff(y));

Ex = zeros(size(Vxy));
Ey = zeros(size(Vxy));

for ny=2:Ly-1
	for nx=2:Lx-1
		Ex(ny,nx) = -(Vxy(ny,nx+1)-Vxy(ny,nx-1))./(2*dx);
		Ey(ny,nx) = -(Vxy(ny+1,nx)-Vxy(ny-1,nx))./(2*dy);
	end
end
