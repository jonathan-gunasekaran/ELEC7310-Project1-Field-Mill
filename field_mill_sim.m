%% ELEC 7310 - Advanced Electrodynamics I - Project 1
% Jonathan Gunasekaran
% 2026/09/19

% Constants
epsilon0 = 8.854e-12; % permittivity of free-space
E = 100;              % static E field of 100 V/m

% numerical iteration until voltage converges
convergence_bound = 0.01;


%% 2D Numerical Simulation of Laplace's Equation
% spatial step delta (meters) --> if constant velocity, essentially equal to v
delta = 0.01;

% bounds of simulation space (meters)
lower_limit = 0;
upper_limit = 4;

% create 2D grid
x = lower_limit:delta:upper_limit;
z = lower_limit:delta:upper_limit;
[X,Z] = meshgrid(x,z);  % x is column, z is row


% determine grid indices for sense plates
% place plates at center of simulation space
% assume the area of each sense plate to be 10cm x 10cm
center_sim_space = (lower_limit + upper_limit) / 2;
plate_width = 0.40; % meters

float_point_tol = 1e-10;
plates_x_inds = find( abs(x-center_sim_space) <= plate_width + float_point_tol); % any points 1m away from center in space along x-axis
plates_z_ind = floor(length(z)/2);
% plates_z_ind = 100;

% get grid indices for each individual plate
plate_a_x_inds = plates_x_inds(1:floor(end/2)+1);
plate_b_x_inds = plates_x_inds(floor(end/2)+1:end);

% determine grid indices for shutter 
% set distance between shutter and sense plates
d = 0.02; % meters
shutter_z_ind = plates_z_ind + d./delta;

% determine grid indices for bottom gnd plane
d_gnd = -0.01;
gnd_x_inds = plates_x_inds;
gnd_z_ind = plates_z_ind + d_gnd./delta;

% shutter_x_ind = find( x-center_sim_space <= 0 & x-center_sim_space > -plate_width ); 


%--------------------------------------------------------------------------
% create initial voltage grid (set to zero everywhere)
V = zeros(size(X));

% tvec = 0:plate_width/delta;
% movement = zeros(size(tvec));
tvec = 0:plate_width/delta*2; % back and forth
movement = [ones(size(plate_a_x_inds)), -ones(size(plate_b_x_inds))];  % 1 -> right; -1 -> left 

prev_shutter_x_start_ind = 1; % starting point of shutter

% Q(t) array for each plate
Qa = zeros(size(tvec));
Qb = zeros(size(tvec));


% create figure for plot of voltage solution to Laplace's equation at each time step
figure; clf;

tic
% compute time-stepped shutter grid indices
for tstep_ind = 1:length(tvec)
    % reset previous iteration voltage grid for each successive timestep so
    % a new solution is obtained at each timestep
    prev_V = zeros(size(V));

    % define initial condition
    % at t=0, assume initial voltage has gradient that is uniform in z across entire simulation space
    % otherwise, just start with solution from previous time step (should minimize convergence time)
    if tstep_ind == 1
        V = repmat((z.*E).', size(z,1), size(x,2));
    end
    % V = repmat((z.*E).', size(z,1), size(x,2));
    
    % enforce boundary condition values at edges of simulation space
    % x is column, z is row
    V(1,:) = 0;                  % bottom of sim space is grounded
    V(end,:) = upper_limit * E;  % top of sim space is held at fixed potential
    V(:,1) = z.*E;               % sides of sim space are held at fixed gradient values
    V(:,end) = z.*E;
    
    % determine time-stepped grid indices for shutter
    % assuming linear motion from one side to the next
    if movement(tstep_ind) == 1
        shutter_x_start_ind = prev_shutter_x_start_ind + 1;
    else
        shutter_x_start_ind = prev_shutter_x_start_ind - 1;
    end
    shutter_x_inds = plates_x_inds(shutter_x_start_ind:shutter_x_start_ind+plate_width/delta);

    % tic
    iteration = 1;
    % make sure that the magnitude of the voltage different between consecutive 
    % iterations is within the convergence bound criteria for all simulation
    % points
    while max(max(abs(V - prev_V))) > convergence_bound
        if mod(iteration, 500) == 0
            fprintf("Currently processing iteration %d for timestep %d/%d\n", iteration, tstep_ind, length(tvec));
        end
    
        prev_V = V;
    
        % Laplace update equation
        for z_ind = 2:size(V,1)-1
            V(z_ind,2:end-1) = 1/4 * (V(z_ind,1:end-2) + V(z_ind,3:end) + ...
                                        V(z_ind-1,2:end-1) + V(z_ind+1,2:end-1)); 
        end
    
        % ground shutter plate positions (force set to zero voltage)
        V(shutter_z_ind, shutter_x_inds) = 0;
        V(gnd_z_ind, gnd_x_inds) = 0;
    
        iteration = iteration + 1;
    end
    % toc

    % imagesc(x,z,V);
    % scatter3(X(:),Z(:),V(:),[],V(:),"filled");
    surf(X,Z,V);
    % contour(X,Z,V);
    hold on;
    plot3(x(shutter_x_inds), repmat(z(shutter_z_ind), size(shutter_x_inds)), repmat(500, size(shutter_x_inds)), 'r', 'LineWidth',2);
    plot3(x(plate_a_x_inds), repmat(z(plates_z_ind), size(plate_a_x_inds)), repmat(500, size(plate_a_x_inds)), 'g', 'LineWidth',2);
    plot3(x(plate_b_x_inds), repmat(z(plates_z_ind), size(plate_b_x_inds)), repmat(500, size(plate_b_x_inds)), 'm', 'LineWidth',2);    
    plot3(x(gnd_x_inds), repmat(z(gnd_z_ind), size(gnd_x_inds)), repmat(500, size(gnd_x_inds)), 'k', 'LineWidth',2);
    drawnow;
    shading interp;
    colormap winter;
    colorbar;
    view(0,90);
    pause(0.001)
    hold off;

    % get flux density from voltage at sense plate grid indices
    % evaluated at boundary between grid indices
    % Da = -epsilon0.*( V(plates_z_ind+1,plate_a_x_inds) - V(plates_z_ind,plate_a_x_inds) );
    % Db = -epsilon0.*( V(plates_z_ind+1,plate_b_x_inds) - V(plates_z_ind,plate_b_x_inds) );
    Da = delta*epsilon0.*( 2*V(plates_z_ind,plate_a_x_inds) - V(plates_z_ind+1,plate_a_x_inds) - V(plates_z_ind-1,plate_a_x_inds) ) + ...
         delta*epsilon0.*( 2*V(plates_z_ind,plate_a_x_inds) - V(plates_z_ind,plate_a_x_inds+1) - V(plates_z_ind,plate_a_x_inds-1) );
    Db = delta*epsilon0.*( 2*V(plates_z_ind,plate_b_x_inds) - V(plates_z_ind+1,plate_b_x_inds) - V(plates_z_ind-1,plate_b_x_inds) ) + ...
         delta*epsilon0.*( 2*V(plates_z_ind,plate_b_x_inds) - V(plates_z_ind,plate_b_x_inds+1) - V(plates_z_ind,plate_b_x_inds-1) );

    % Da = delta.*epsilon0.*( V(plates_z_ind+1,plate_a_x_inds) - V(plates_z_ind-1,plate_a_x_inds) ) + ...
    %      delta.*epsilon0.*( V(plates_z_ind,plate_a_x_inds+1) - V(plates_z_ind,plate_a_x_inds-1) );
    % Db = delta.*epsilon0.*( V(plates_z_ind+1,plate_b_x_inds) - V(plates_z_ind-1,plate_b_x_inds) ) + ...
    %      delta.*epsilon0.*( V(plates_z_ind,plate_b_x_inds+1) - V(plates_z_ind,plate_b_x_inds-1) );


    [Ex, Ez] = gradient(-V, delta);
    Dx = epsilon0*Ex;
    Dz = epsilon0*Ez;
    Da = Dz(plates_z_ind, plate_a_x_inds);
    Db = Dz(plates_z_ind, plate_b_x_inds);

    % integrate over the area of the plates to find the total charge
    Qa(tstep_ind) = sum(Da);
    Qb(tstep_ind) = sum(Db);

end
toc

I_diff = diff(Qb) - diff(Qa);

% time vector is arbitrary
gain = 1e10;
Vout_numerical = gain.*I_diff;

figure; clf;
plot(tvec(1:end-1), Vout_numerical);
hold on;


%% Analytical Solution for Amplifier Output
vel = delta;

Vout_analytical = 2*vel*epsilon0*E * gain;
Vout_analytical = movement .* Vout_analytical;
plot(tvec, repmat(Vout_analytical,size(tvec)));

xlabel("Time Step (t)")


%% Comparison
