%% ELEC 7310 - Advanced Electrodynamics I - Project 1
% Jonathan Gunasekaran
% 2026/09/19

%% 2D Numerical Simulation of Laplace's Equation
% spatial step delta (meters)
delta = 0.01;

% bounds of simulation space (meters)
lower_limit = 0;
upper_limit = 10;

% create 2D grid
x = lower_limit:delta:upper_limit;
z = lower_limit:delta:upper_limit;
[X,Z] = meshgrid(x,z);  % x is column, z is row

% place plates at center of simulation space
% assume the area of each sense plate to be 10cm x 10cm
center_sim_space = (lower_limit + upper_limit) / 2;
plate_width = 0.1; % meters
plates_x_ind = find( abs(x-center_sim_space) <= plate_width ); % any points 1m away from center in space along x-axis
plates_z_ind = floor(length(z)/2);

% set distance between shutter and sense plates
d = 0.01; % meters
shutter_z_ind = plates_z_ind + d./delta;

% compute time-stepped shutter grid indices
shutter_x_ind = find( x-center_sim_space <= 0 & x-center_sim_space > -plate_width ); 
tstep_ind = 1;
shutter_x_ind = plates_x_ind(tstep_ind:tstep_ind+plate_width/delta);



% set initial boundary condition values for edges of simulation space
V = zeros(size(X));
E = 100;  % static E field of 100 V/m

% x is column, z is row
V(1,:) = 0;                  % bottom of sim space is grounded
V(end,:) = upper_limit * E;  % top of sim space is held at fixed potential
V(:,1) = z.*E;               % sides of sim space are held at fixed gradient values
V(:,end) = z.*E;


% numerical iteration until voltage converges
prev_V = zeros(size(X));
convergence_bound = 0.1;

iteration = 1;
while max(max(V - prev_V)) > convergence_bound
    if mod(iteration, 100) == 0
        fprintf("Currently processing iteration %d\n", iteration);
    end

    prev_V = V;

    for z_ind = 2:size(V,2)-1
        V(2:end-1, z_ind) = 1/4 * (V(1:end-2,z_ind) + V(3:end,z_ind) + ...
                                    V(2:end-1,z_ind-1) + V(2:end-1,z_ind+1)); 
    end

    % ground shutter plate positions
    V(shutter_z_ind, shutter_x_ind) = 0;

    iteration = iteration + 1;
end

% get flux density from voltage


% integrate over the area of the plates to find the total charge






%% Analytical Solution for Amplifier Output

