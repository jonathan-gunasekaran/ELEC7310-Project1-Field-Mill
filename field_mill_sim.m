%% ELEC 7310 - Advanced Electrodynamics I - Project 1
% Jonathan Gunasekaran
% 2026/09/19

%% 2D Numerical Simulation of Laplace's Equation
% spatial step delta (meters)
delta = 0.01;

% bounds of simulation space (meters)
lower_limit = 0;
upper_limit = 5;

% create 2D grid
x = lower_limit:delta:upper_limit;
z = lower_limit:delta:upper_limit;
[X,Z] = meshgrid(x,z);  % x is column, z is row


% determine grid indices for sense plates
% place plates at center of simulation space
% assume the area of each sense plate to be 10cm x 10cm
center_sim_space = (lower_limit + upper_limit) / 2;
plate_width = 0.10; % meters
plates_x_ind = find( abs(x-center_sim_space) <= plate_width ); % any points 1m away from center in space along x-axis
plates_z_ind = floor(length(z)/2);
% plates_z_ind = 100;


% determine grid indices for shutter 
% set distance between shutter and sense plates
d = 0.02; % meters
shutter_z_ind = plates_z_ind + d./delta;

% compute time-stepped shutter grid indices
shutter_x_ind = find( x-center_sim_space <= 0 & x-center_sim_space > -plate_width ); 
tstep_ind = 1;
shutter_x_ind = plates_x_ind(tstep_ind:tstep_ind+plate_width/delta);


% define initial condition
E = 100;  % static E field of 100 V/m

% assume initial voltage is zero everywhere
V = zeros(size(X));

% assume initial voltage has gradient that is uniform in z across entire simulation space
% V = repmat((z.*E).', size(z,1), size(x,2));


% enforce boundary condition values at edges of simulation space
% x is column, z is row
V(1,:) = 0;                  % bottom of sim space is grounded
V(end,:) = upper_limit * E;  % top of sim space is held at fixed potential
V(:,1) = z.*E;               % sides of sim space are held at fixed gradient values
V(:,end) = z.*E;



% numerical iteration until voltage converges
prev_V = zeros(size(X));
convergence_bound = 0.01;

tic
iteration = 1;
% make sure that the magnitude of the voltage different between consecutive 
% iterations is within the convergence bound criteria for all simulation
% points
while max(max(abs(V - prev_V))) > convergence_bound
    if mod(iteration, 100) == 0
        fprintf("Currently processing iteration %d\n", iteration);
    end

    prev_V = V;

    % Laplace update equation
    for z_ind = 2:size(V,1)-1
        V(z_ind,2:end-1) = 1/4 * (V(z_ind,1:end-2) + V(z_ind,3:end) + ...
                                    V(z_ind-1,2:end-1) + V(z_ind+1,2:end-1)); 
    end

    % ground shutter plate positions (force set to zero voltage)
    V(shutter_z_ind, shutter_x_ind) = 0;

    iteration = iteration + 1;
end
toc

figure; clf;
% imagesc(x,z,V);
% scatter3(X(:),Z(:),V(:),[],V(:),"filled");
surf(X,Z,V);
hold on;
plot3(x(shutter_x_ind), repmat(z(shutter_z_ind), size(shutter_x_ind)), repmat(-1, size(shutter_x_ind)), 'r', 'LineWidth',3);
plot3(x(plates_x_ind), repmat(z(plates_z_ind), size(plates_x_ind)), repmat(-1, size(plates_x_ind)), 'g', 'LineWidth',3);
shading interp;
colorbar;
view(0,90);

% get flux density from voltage


% integrate over the area of the plates to find the total charge






%% Analytical Solution for Amplifier Output



%% Comparison
