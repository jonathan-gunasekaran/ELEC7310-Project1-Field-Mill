%% ELEC 7310 - Advanced Electrodynamics I - Project 1
% Jonathan Gunasekaran
% 2026/09/19

% TODO:
% - plot E as quiver
% - show diff in E between analytical and numerical along sense plates axis
% - show variation as d increases, maybe also as A decreases if time allows
% - truncate plot space to just show field around sense plates
% - [DONE] initial condition ground the sense plates

clear all;

% Constants
epsilon0 = 8.854e-12; % permittivity of free-space
E = 100;              % static E field of 100 V/m

% Simulation Options
% numerical iteration until voltage converges to within preset tolerance bound
convergence_bound = 0.01;
sim_full_cycle = 1;  % if set, move shutter right and left over full cycle; otherwise, only half cycle


%% 2D Numerical Simulation of Laplace's Equation
% spatial step delta (meters)
delta_xyz = 0.01;
% time step (seconds) - somewhat arbitrary since velocity of shutter is assumed constant
delta_t = 1;

% bounds of simulation space (meters)
lower_spatial_limit = 0;
upper_spatial_limit = 4;

% create 2D grid
x = lower_spatial_limit:delta_xyz:upper_spatial_limit;
z = lower_spatial_limit:delta_xyz:upper_spatial_limit;
[X,Z] = meshgrid(x,z);  % x is column, z is row


% determine grid indices for sense plates
% place plates at center of simulation space
% assume the area of each sense plate to be square (width x width)
center_sim_space = (lower_spatial_limit + upper_spatial_limit) / 2;
plate_width = 0.20; % meters

% any points along x-axis that are within "plate_width" distance away from center grid index
float_point_tol_indexing = 1e-10;
sense_x_inds = find( abs(x-center_sim_space) <= plate_width + float_point_tol_indexing);
sense_z_ind = floor(length(z)/2);

% get grid indices for each individual plate
sense_a_x_inds = sense_x_inds(1:floor(end/2)+1);
sense_b_x_inds = sense_x_inds(floor(end/2)+1:end);

% determine grid indices for shutter 
% set distance between shutter and sense plates
dz_shutter = 0.02; % meters
shutter_z_ind = sense_z_ind + dz_shutter./delta_xyz;

% determine grid indices for bottom gnd plane
dz_gnd = -0.01;
gnd_x_inds = sense_x_inds;
gnd_z_ind = sense_z_ind + dz_gnd./delta_xyz;


%--------------------------------------------------------------------------
% create initial voltage grid (set to zero everywhere)
V = zeros(size(X));

% define time stamps and shutter movement pattern for simulation
if sim_full_cycle
    tvec = 0:plate_width/delta_xyz*2; % back and forth
    
    % NOTE: movement array represents motion direction between successive
    % time stamps given by tvec, so len(movement) = len(tvec) - 1
    movement = [diff(sense_a_x_inds), -diff(sense_a_x_inds)];  % 1 -> right; -1 -> left 
else
    tvec = 0:plate_width/delta_xyz;
    movement = diff(sense_a_x_inds);
end

shutter_x_start_ind = 1; % starting point of shutter

% Q(t) array for each plate
Qa = zeros(size(tvec));
Qb = zeros(size(tvec));


% create figure for movie plot of voltage solution to Laplace's equation at each time step
field_mill_movie_fig = figure; clf;

tic
% compute time-stepped shutter grid indices
for tstep_ind = 1:length(tvec)
    % reset previous iteration voltage grid for each successive timestep so
    % a new solution is obtained at each timestep
    prev_V = zeros(size(V));

    % define initial condition
    % assume initial solution is uniform gradient for each shutter position
    % (iterate until Laplace's equation converges to within tolerance criteria)
    V = repmat((z.*E).', size(z,1), size(x,2));

    % at t=0, assume initial voltage has gradient that is uniform in z
    % across entire simulation space; otherwise, just start with solution
    % from previous time step (should minimize convergence time) 
    % -->  DOES NOT WORK FOR WHATEVER REASON
    % if tstep_ind == 1
    %     V = repmat((z.*E).', size(z,1), size(x,2));
    % end
    
    % enforce boundary condition values at edges of simulation space
    % x is column, z is row
    V(1,:) = 0;                  % bottom of sim space is grounded
    V(end,:) = upper_spatial_limit * E;  % top of sim space is held at fixed potential
    V(:,1) = z.*E;               % sides of sim space are held at fixed gradient values
    V(:,end) = z.*E;
    
    % ground sense plates through pull-down at start of evaluation of
    % Laplace's equation (not force set throughout iteration)
    V(sense_z_ind,sense_a_x_inds) = 0;
    V(sense_z_ind,sense_b_x_inds) = 0;
    
    % determine time-stepped grid indices for shutter
    % disp(shutter_x_start_ind);
    shutter_x_inds = sense_x_inds(shutter_x_start_ind:shutter_x_start_ind+plate_width/delta_xyz);

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


    [Ex, Ez] = gradient(-V, delta_xyz);
    % Dx = epsilon0*Ex;
    % Dz = epsilon0*Ez;
    % Da = Dz(plates_z_ind, plate_a_x_inds);
    % Db = Dz(plates_z_ind, plate_b_x_inds);

    % plot movie of scalar voltage as shutter moves back and forth
    figure(field_mill_movie_fig);
    axis equal;
    % imagesc(x,z,V);
    % scatter3(X(:),Z(:),V(:),[],V(:),"filled");
    % surf(X,Z,V);
    V_level_diff = 20;
    equiV_levels = lower_spatial_limit:V_level_diff:upper_spatial_limit*E;
    contour(X,Z,V,equiV_levels);
    % quiver(X,Z,Ex,Ez);
    hold on;
    plot3(x(shutter_x_inds), repmat(z(shutter_z_ind), size(shutter_x_inds)), repmat(500, size(shutter_x_inds)), 'r', 'LineWidth',2);
    plot3(x(sense_a_x_inds), repmat(z(sense_z_ind), size(sense_a_x_inds)), repmat(500, size(sense_a_x_inds)), 'y', 'LineWidth',2);
    plot3(x(sense_b_x_inds), repmat(z(sense_z_ind), size(sense_b_x_inds)), repmat(500, size(sense_b_x_inds)), 'm', 'LineWidth',2);    
    plot3(x(gnd_x_inds), repmat(z(gnd_z_ind), size(gnd_x_inds)), repmat(500, size(gnd_x_inds)), 'k', 'LineWidth',2);
    drawnow;
    shading interp;
    colormap winter;
    colorbar;
    view(0,90);
    pause(0.001)
    hold off;

    % get flux density from voltage at sense plate grid indices
    % negative gradient of voltage evaluated at boundary between grid indices
    % testing with only top face of sense plates
    % Da = -epsilon0.*( V(plates_z_ind+1,plate_a_x_inds) - V(plates_z_ind,plate_a_x_inds) );
    % Db = -epsilon0.*( V(plates_z_ind+1,plate_b_x_inds) - V(plates_z_ind,plate_b_x_inds) );
    
    % evaluate gradient for all boundaries of grid cell (closed surface integral)
    Da = epsilon0.*( 2*V(sense_z_ind,sense_a_x_inds) - V(sense_z_ind+1,sense_a_x_inds) - V(sense_z_ind-1,sense_a_x_inds) ) + ...
         epsilon0.*( 2*V(sense_z_ind,sense_a_x_inds) - V(sense_z_ind,sense_a_x_inds+1) - V(sense_z_ind,sense_a_x_inds-1) );
    Db = epsilon0.*( 2*V(sense_z_ind,sense_b_x_inds) - V(sense_z_ind+1,sense_b_x_inds) - V(sense_z_ind-1,sense_b_x_inds) ) + ...
         epsilon0.*( 2*V(sense_z_ind,sense_b_x_inds) - V(sense_z_ind,sense_b_x_inds+1) - V(sense_z_ind,sense_b_x_inds-1) );
    
    % scale factor for plate area since only doing simulation in 2D
    Da = delta_xyz.*Da;
    Db = delta_xyz.*Db;

    % integrate over the area of the plates to find the total charge
    Qa(tstep_ind) = sum(Da);
    Qb(tstep_ind) = sum(Db);


    % update start index of shutter for next iteration
    % assuming linear motion from one side to the next
    if tstep_ind ~= length(tvec)
        if movement(tstep_ind) == 1
            shutter_x_start_ind = shutter_x_start_ind + 1;
        else
            shutter_x_start_ind = shutter_x_start_ind - 1;
        end
    end
end
toc

I_diff = diff(Qb) - diff(Qa);

% time vector is arbitrary (uniform time step)
% gain is arbitrary
gain = 1e10;
Vout_numerical = gain.*I_diff;


% create plot for charge on the two plates as a function of time
q_fig = figure; clf;
xlabel("Time Step (t)");
ylabel("Relative Q (C)");
plot(tvec, Qa*gain, 'DisplayName', 'relative Qa');
hold on;
plot(tvec, Qb*gain, 'DisplayName', 'relative Qb');

% create plot for output voltage from transimpedance amplifier
vout_fig = figure; clf;
xlabel("Time Step (t)")
ylabel("Relative V_{out}")

plot(tvec(1:end-1), Vout_numerical, 'DisplayName','V_{numerical}', 'LineWidth',2);
hold on;



%% Analytical Solution for Amplifier Output
vel = delta_xyz./delta_t;

Vout_analytical = 2*vel*epsilon0*E * gain;
Vout_analytical = Vout_analytical * delta_xyz; % account for 2D simulation - depth factor for area
Vout_analytical = movement .* Vout_analytical;

% add analytical solution to voltage plot
figure(vout_fig);
plot(tvec(1:end-1), Vout_analytical, 'DisplayName','V_{analytical}', 'LineWidth',2);
hold on;


%% Comparison

residualEz = Ez - E;

% residual_E_quiver
