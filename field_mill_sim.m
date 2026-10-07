%% ELEC 7310 - Advanced Electrodynamics I - Project 1
% Jonathan Gunasekaran
% 2026/09/19

% TODO:
% - [DONE] plot E as quiver
% - show diff in E between analytical and numerical along sense plates axis
% - show variation as d increases, maybe also as A decreases if time allows (time does not allow)
% - [DONE] truncate plot space to just show field around sense plates

clear all;

% Constants
epsilon0 = 8.854e-12; % permittivity of free-space
E = 100;              % static E field of 100 V/m

% Simulation Options
% numerical iteration until voltage converges to within preset tolerance bound
convergence_bound = 0.05;
sim_full_cycle = 1;  % if set, move shutter right and left over full cycle; otherwise, only half cycle

% plotting toggles
zoom_in_movie = 1;
plot_charge = 1;
plot_current = 1;
plot_voltage_out = 1;
plot_analytical = 0;
plot_comparison = 1;

%% 2D Numerical Simulation of Laplace's Equation
% spatial step delta (meters)
delta_xyz = 0.01;
% time step (seconds) - somewhat arbitrary since velocity of shutter is assumed constant
delta_t = 1;

% bounds of simulation space (meters)
lower_spatial_limit = -2;
upper_spatial_limit = 2;

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
sense_z_ind = floor(length(z)/2)+1;

% get grid indices for each individual plate
sense_a_x_inds = sense_x_inds(1:floor(end/2)+1);
sense_b_x_inds = sense_x_inds(floor(end/2)+1:end);

% determine grid indices for shutter 
% set distance between shutter and sense plates
dz_shutter = 0.02; % meters
shutter_z_ind = sense_z_ind + dz_shutter./delta_xyz;

% determine grid indices for bottom gnd plane
dz_gnd = -0.02;
gnd_x_inds = sense_x_inds;
gnd_z_ind = sense_z_ind + dz_gnd./delta_xyz;

% shift up so gnd plane is at center z index
gnd_z_ind = gnd_z_ind - dz_gnd./delta_xyz;
shutter_z_ind = shutter_z_ind - dz_gnd./delta_xyz;
sense_z_ind = sense_z_ind - dz_gnd./delta_xyz;


%-----------------------------------------------------------------------
% create initial voltage grid (set to zero everywhere)
V = zeros(size(X));

T = 2*plate_width/delta_xyz*delta_t; % period of full cycle

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
Qa_numerical = zeros(size(tvec));
Qb_numerical = zeros(size(tvec));


% store E-field at start and middle of shutter movement for comparison
% against analytical assumptions



% create figure for movie plot of voltage solution to Laplace's equation at each time step
field_mill_movie_fig = figure; clf;
if zoom_in_movie
    plot_size = 1.5*plate_width;
    movie_x_inds = find( abs(x-center_sim_space) <= plot_size + float_point_tol_indexing);
    movie_z_inds = find( abs(z-center_sim_space) <= plot_size + float_point_tol_indexing);
    movie_x_bounds = x([movie_x_inds(1), movie_x_inds(end)]);
    movie_z_bounds = z([movie_z_inds(1), movie_z_inds(end)]);
end

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
    V(1,:) = lower_spatial_limit * E;   % bottom of sim space is held at fixed potential
    V(end,:) = upper_spatial_limit * E; % top of sim space is held at fixed potential
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
    Dx = epsilon0*Ex;
    Dz = epsilon0*Ez;
    Da = Dz(sense_z_ind, sense_a_x_inds);
    Db = Dz(sense_z_ind, sense_b_x_inds);
    Da_dot_ds = Da * delta_xyz.^2;
    Db_dot_ds = Db * delta_xyz.^2;


    %-------------------------------------------------------------------
    % plot movie of scalar voltage as shutter moves back and forth
    figure(field_mill_movie_fig);
    axis equal;
    % surf(X,Z,V);
    V_level_diff = 20; % voltage difference between adjacent equipotential lines
    equiV_levels = lower_spatial_limit*E:V_level_diff:upper_spatial_limit*E;
    contour(X,Z,V,equiV_levels);
    hold on;
    quiver(X,Z,Ex,Ez);
    plot3(x(shutter_x_inds), repmat(z(shutter_z_ind), size(shutter_x_inds)), repmat(500, size(shutter_x_inds)), 'k', 'LineWidth',2);
    plot3(x(sense_a_x_inds), repmat(z(sense_z_ind), size(sense_a_x_inds)), repmat(500, size(sense_a_x_inds)), 'r', 'LineWidth',2);
    plot3(x(sense_b_x_inds), repmat(z(sense_z_ind), size(sense_b_x_inds)), repmat(500, size(sense_b_x_inds)), 'g', 'LineWidth',2);    
    plot3(x(gnd_x_inds), repmat(z(gnd_z_ind), size(gnd_x_inds)), repmat(500, size(gnd_x_inds)), 'k', 'LineWidth',2);
    
    if zoom_in_movie
        xlim(movie_x_bounds);
        ylim(movie_z_bounds);
    end
    xlabel("x (m)");
    ylabel("z (m)");
    title("Potential Evolution over Time for 2D Cross Section")
    drawnow;
    shading interp;
    colormap winter;
    colorbar;
    view(0,90);
    pause(0.001)
    hold off;
    %-------------------------------------------------------------------


    % % get flux density from voltage at sense plate grid indices - 
    % % evaluate negative gradient of voltage at all boundaries (between 
    % % grid indices) of each grid cell --> closed surface integral
    % Da_dot_ds = delta_xyz.*epsilon0.*( 2*V(sense_z_ind,sense_a_x_inds) - V(sense_z_ind+1,sense_a_x_inds) - V(sense_z_ind-1,sense_a_x_inds) ) + ...
    %             delta_xyz.*epsilon0.*( 2*V(sense_z_ind,sense_a_x_inds) - V(sense_z_ind,sense_a_x_inds+1) - V(sense_z_ind,sense_a_x_inds-1) );
    % Db_dot_ds = delta_xyz.*epsilon0.*( 2*V(sense_z_ind,sense_b_x_inds) - V(sense_z_ind+1,sense_b_x_inds) - V(sense_z_ind-1,sense_b_x_inds) ) + ...
    %             delta_xyz.*epsilon0.*( 2*V(sense_z_ind,sense_b_x_inds) - V(sense_z_ind,sense_b_x_inds+1) - V(sense_z_ind,sense_b_x_inds-1) );
    
    % scale factor for plate area since only doing simulation in 2D
    %  - depth integral becomes const multiplication factor (assuming 2D
    %    cross sections at any depth will have the same developed charge 
    %    distribution) 
    Da_dot_ds = (plate_width/delta_xyz+1).*Da_dot_ds;
    Db_dot_ds = (plate_width/delta_xyz+1).*Db_dot_ds;

    % integrate over the area of the plates to find the total charge
    Qa_numerical(tstep_ind) = sum(Da_dot_ds);
    Qb_numerical(tstep_ind) = sum(Db_dot_ds);


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

I_diff_numerical = -(diff(Qb_numerical) - diff(Qa_numerical))./delta_t;

% time vector is arbitrary (uniform time step)
% gain is arbitrary
gain = 1e13;
Vout_numerical = gain.*I_diff_numerical;


% create plot for charge on the two plates as a function of time
q_C2nC_factor = 1e9;
if plot_charge
    q_numerical_fig = figure; clf;
    plot(tvec, Qa_numerical*q_C2nC_factor, 'DisplayName', 'Qa', LineWidth=2);
    hold on;
    plot(tvec, Qb_numerical*q_C2nC_factor, 'DisplayName', 'Qb', LineWidth=2);
    
    xlabel("Time Step (s)");
    ylabel("Total Charge (nC)");
    title("Numerical Charge on Sense Plates vs. Time");
    legend;
end

% create plot for output voltage from transimpedance amplifier
if plot_voltage_out
    vout_numerical_fig = figure; clf;
    plot(tvec(1:end-1), Vout_numerical, 'DisplayName','Numerical', 'LineWidth',2);
    hold on;
    
    xlabel("Time Step (s)");
    ylabel("Output Voltage (V)");
    title("Output Voltage vs. Time");
end

%% Analytical Solution for Amplifier Output
vel = delta_xyz./delta_t;

% initial steady-state charge values with shutter over plate A
Qa_analytical = zeros(size(tvec));
Qb_analytical = ones(size(tvec)) * -epsilon0*E * plate_width^2; % account for 2D simulation - depth factor for area

for tind = 2:length(tvec)
    curr_motion = movement(tind-1);
    if curr_motion == 1 % to the right: negative charge flows into plate A, flows out of plate B
        Qa_analytical(tind) = Qa_analytical(tind-1) + (-epsilon0*E * plate_width * vel*delta_t);
        Qb_analytical(tind) = Qb_analytical(tind-1) - (-epsilon0*E * plate_width * vel*delta_t);
    else                % to the left: negative charge flows into plate B, flows out of plate A
        Qa_analytical(tind) = Qa_analytical(tind-1) - (-epsilon0*E * plate_width * vel*delta_t);
        Qb_analytical(tind) = Qb_analytical(tind-1) + (-epsilon0*E * plate_width * vel*delta_t);
    end
end

I_diff_analytical = -(diff(Qb_analytical) - diff(Qa_analytical))./delta_t;
Vout_analytical = I_diff_analytical * gain;


% Plotting
if plot_analytical
    if plot_charge
        q_analytical_fig = figure; clf;
        plot(tvec, Qa_analytical*q_C2nC_factor, 'DisplayName', 'Qa', LineWidth=2);
        hold on;
        plot(tvec, Qb_analytical*q_C2nC_factor, 'DisplayName', 'Qb', LineWidth=2);
        
        xlabel("Time Step (s)");
        ylabel("Total Charge (nC)");
        title("Analytical Charge on Sense Plates vs. Time");
        legend('Location','best');
    end
    
    if plot_current
        i_analytical_fig = figure; clf;
        plot(tvec(1:end-1), -diff(Qa_analytical*q_C2nC_factor)./delta_t, 'DisplayName', 'Ia', LineWidth=2);
        hold on;
        plot(tvec(1:end-1), -diff(Qb_analytical*q_C2nC_factor)./delta_t, 'DisplayName', 'Ib', LineWidth=2);
        
        xlabel("Time Step (s)");
        ylabel("Output Current (nA)");
        title("Analytical Output Current from Sense Plates vs. Time");
        legend('Location','best');
    end
    
    % add analytical solution to voltage plot
    if plot_voltage_out
        vout_analytical_fig = figure; clf;
        plot(tvec(1:end-1), Vout_analytical, 'DisplayName','Analytical', 'LineWidth',2);
        hold on;
        
        % Vout_analytical2 = -4 * epsilon0*E * plate_width^2 / T;
        % Vout_analytical2 = Vout_analytical2 * movement;
        % Vout_analytical2 = Vout_analytical2 * gain;

        xlabel("Time Step (s)");
        ylabel("Output Voltage (V)");
        title("Analytical Output Voltage vs. Time");
    end
end


if plot_voltage_out
    figure(vout_numerical_fig);
    plot(tvec(1:end-1), Vout_analytical, '--', 'DisplayName','Analytical', 'LineWidth',2);
    legend('Location','best');
end

%% Comparison

residualEz = Ez - E;

% residual_E_quiver
