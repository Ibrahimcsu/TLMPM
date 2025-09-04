% TL-MPM 2D Analysis Code for Excitation Function Comparison
% ======================================================
%
% Description:
% ------------
%   This script models the transverse vibration response of a beam subjected 
%   to prescribed excitation functions (e.g., sine, cosine, square-wave, chirp). 
%   The excitation is applied at boundary nodes or distributed along the beam, 
%   and the dynamic response (displacement, velocity, acceleration) is computed 
%   at selected points such as the beam tip or midspan.
%
% Features:
% ---------
% - Linear elastic material model
% - Explicit time integration
% - Standard bilinear Shape Functions
% - Total Lagrangian formulation
%
% Usage:
% ------
% 1. Define material properties (E, density, nu)
% 2. Create particle mesh and Eulerian grid
% 3. Set initial conditions and boundary conditions
% 4. Run simulation with appropriate time step
%
% Particle Mesh Creation:
% ----------------------
% If useing the Particle mesh generator: 
% xp = array of particle coordinates [npart × 2]
% Create particles within a rectangular domain with:
% number_of_particlesx_per_cell = particles per cell in x-direction
% number_of_particlesy_per_cell = particles per cell in y-direction
% small_width = cell size
% You may upload or define any particle mesh in the following format.
% Input: xp (N_p × 2) matrix
%        where N_p = number of particles
%        and the two columns represent the x- and y-coordinates of each particle.
%
% Eulerian Mesh Creation:
% ----------------------
% LOC = nodal coordinates [nnodes × 2]
% Create grid with:
% nnodesx = number of nodes in x-direction
% nnodesy = number of nodes in y-direction
% dxn, dyn = grid spacing in x and y directions
%
% Author: Ibrahim Bouzaid, Ph.D. Student
% Institution: Department of Civil and Environmental Engineering,
%              Colorado State University, Fort Collins, CO, USA
% Advisor: Prof. Paul H. Heyliger
%
% Authorship and Contributions:
%   This code is derived from original TL-MPM implementations 
%   provided by Prof. Paul H. Heyliger (Fortran and MATLAB versions).
%   The present work by is prepared by Ibrahim Bouzaid, introduces modifications 
%   to the mapping procedure and shape function evaluations, incorporates a fully 
%   vectorized structure to improve efficiency, and extends the formulation 
%   in line with the Total Lagrangian MPM as described in the published literature.
%
% Reference:
% -----------
% de Vaucorbeil, A., Nguyen, V. P., & Hutchinson, C. R. (2020). 
% A total-Lagrangian material point method for solid mechanics problems 
% involving large deformations. Computer Methods in Applied Mechanics 
% and Engineering, 360, 112783.
% -----------
% Version: 1.0
% Date: August 2025
%==========================================================================
%==========================================================================
clear all; close all; clc;
format long

% List of excitation types to test
excitation_types = {'cosine', 'sine', 'Switching'};

% Initialize results structure
results = struct();

% Common parameters
E = 2e8;
density = 2800;
nu = 0;
c11 = E/(1-nu^2);
c22 = c11;
c12 = c11*nu;
c21 = c12;
c66 = E/(2*(1+nu));
C = [c11 c12 0; c12 c22 0; 0 0 c66];

LLx = 15; LLy = 1;
dxn = 1; dyn = 1;
le = [dxn, dyn];
X = 0:dxn:LLx; Y = 0:dyn:LLy;
[x1, y1] = meshgrid(Y, X);
y = reshape(x1, length(Y)*length(X), 1);
x = reshape(y1, length(X)*length(Y), 1);
LOC = [x, y];
nnodesx = length(X); nnodesy = length(Y);
nnodes = length(X)*length(Y);
nem = (LLx/dxn)*(LLy/dyn);

% Particles setup
width = 15; height = 1;
number_of_particlesx_per_cell = 2;
number_of_particlesy_per_cell = 2;
T_num_particles_x = round(width*number_of_particlesx_per_cell/dxn);
T_num_particles_y = round(height*number_of_particlesy_per_cell/dxn);
small_width = width / T_num_particles_x;
small_height = height / T_num_particles_y;
centers = zeros(T_num_particles_x * T_num_particles_y, 2);
index = 1;
for i = 1:T_num_particles_x
    for j = 1:T_num_particles_y
        x_center = (i - 0.5) * small_width;
        y_center = (j - 0.5) * small_height;
        centers(index, :) = [x_center, y_center];
        index = index + 1;
    end
end
xp(:,1) = centers(:,1) + 0;
xp(:,2) = centers(:,2) + 0;

npart = length(xp(:,1));
vpart = zeros(npart, 2);

for i = 1:npart
    r1(i,:) = [(small_width/number_of_particlesx_per_cell)/2 0];
    r2(i,:) = [0 (small_width/number_of_particlesy_per_cell)/2];
end
r10 = r1; r20 = r2;

volpart = zeros(1, npart);
for i = 1:npart
    volpart(i) = 4*abs(r1(i,1)*r2(i,2)-r1(i,2)*r2(i,1));
end
partmass = density * volpart;
volpart0 = volpart;

% Time parameters
dt = 0.001;
T = 20; % Reduced time for testing
Tstep = round(T/dt);
time = (0:Tstep-1)*dt;

% Excitation parameters
g = 20000; f0 = 0.5; phi = 0; Nb = 8; Tb = Nb / f0;
t0 = 0.1; Aimp = 5000;

% Precompute shape functions
sff = zeros(nnodes, npart);
dsf = zeros(2, nnodes, npart);
for p = 1:npart
    for i = 1:nnodes
        [sf, dsfx, dsfy] = shapef(xp(p,1), xp(p,2), LOC(i,1), LOC(i,2), dxn, dyn);
        sff(i,p) = sf;
        dsf(1,i,p) = dsfx;
        dsf(2,i,p) = dsfy;
        if sff(i,p) == 0
            dsf(1,i,p) = 0;
            dsf(2,i,p) = 0;
        end
    end
end
nodmass = partmass * sff';
dsfxx = squeeze(dsf(1,:,:));
dsfyy = squeeze(dsf(2,:,:));

% Main loop over excitation types
for exc_idx = 1:length(excitation_types)
    excitation_type = excitation_types{exc_idx};
    fprintf('Running simulation for: %s\n', excitation_type);

    % Initialize variables for this simulation
    stress = zeros(2,2,npart);
    strain = zeros(2,2,npart);
    defgrad = repmat(eye(2), 1, 1, npart);
    velgrad = zeros(2,2,npart);
    vpart = zeros(npart,2);

    ebc = 1;
    if ebc == 1
        ebcn = 1:nnodesx:nnodes;
    else
        ebcn = [];
    end

    % Preallocate arrays for storing results
    Ken = zeros(1, Tstep);
    Sen = zeros(1, Tstep);
    sigma=zeros(2,2,npart);
    Ken(Tstep)=0;Sen(Tstep)=0;Pen(Tstep)=0;strs(2,2,Tstep,npart)=0;strn(2,2,Tstep,npart)=0;defgrad=eye(2).*ones(2,2,npart);sf = zeros(nnodes,npart);dsf=zeros(2,nnodes,npart);forcenodetotal=zeros(nnodes,2);fint=zeros(Tstep,nnodes,2);fext=zeros(Tstep,nnodes,2);ftot=zeros(Tstep,nnodes,2);trac=zeros(Tstep,nnodes,2);ptraction=zeros(npart,2);traction=zeros(nnodes,2);ptr=zeros(npart,2,Tstep);sigma=zeros(2,2,npart);nomo(Tstep,nnodes,2)=0;noms(Tstep,nnodes,1)=0;vn(Tstep,nnodes,2)=0;pm(Tstep,npart)=0;Fdot=zeros(2,2,npart);

    tip_displacement = zeros(1, Tstep);
    excitation_signal = zeros(1, Tstep);

    % Switching function initialization
    initial_switch_interval = 150;
    interval_increment = 50;
    current_switch_interval = initial_switch_interval;
    current_value = -1;
    switch_counter = 0;

    % Main time loop
    for nstep = 1:Tstep
        if mod(nstep, 100) == 0
            fprintf('  Step %d/%d/%d\n', exc_idx, nstep, Tstep);
        end

        % Reset variables
        nodemo = zeros(nnodes, 2);
        vnode = zeros(nnodes, 2);
        forcenodee = zeros(nnodes, 2);
        forcenodei = zeros(nnodes, 2);
        forcenodetotal = zeros(nnodes, 2);
        strs11(:,:)=stress(1,1,:);
        strs22(:,:)=stress(2,2,:);
        strs12(:,:)=stress(1,2,:);


        t = (nstep-1)*dt;

        % Apply excitation
        ptraction = zeros(npart, 2);
        t = (nstep-1)*dt;   % current simulation time
        switch excitation_type
            case 'cosine'
                ptraction(1:2:8, 2) = g * cos(2*pi*f0*t);

            case 'sine'
                ptraction(1:2:8, 2) = g * sin(2*pi*f0*t);

            case 'Switching'
                % Switching function update
                switching_function(nstep) = current_value;
                switch_counter = switch_counter + 1;
                if switch_counter >= current_switch_interval
                    current_value = -current_value;
                    switch_counter = 0;
                    current_switch_interval = current_switch_interval + interval_increment;
                end
                ptraction(1:2:8,2) = current_value * g;
            otherwise
                ptraction(1:2:8, 2) = 0; % default: no load
        end


        excitation_signal(nstep) = ptraction(1,2);

        frate=[0, 0];
        nodemo=(partmass'.*(vpart))'*sff';
        forcenodee=(frate'*nodmass)';
        forcenodeix=((-volpart'.*strs11)'*dsfxx'+(-volpart'.*strs12)'*dsfyy');
        forcenodeiy=((-volpart'.*strs12)'*dsfxx'+(-volpart'.*strs22)'*dsfyy');
        forcenodei=[forcenodeix ;forcenodeiy]';
        traction=((volpart'.*ptraction)'*sff')';

        % Impulse Momentum at the nodes
        forcenodetotal=forcenodee+forcenodei+traction;
        nodemo=(nodemo+dt*forcenodetotal')';

        % Update the Particle Velocities
        anode=((1e10*forcenodetotal)'./(1e10*nodmass))';
        vnode=((1e10*nodemo)'./(1e10*nodmass))';


        if ebc==1
            vnode(ebcn,:)=0;
            anode(ebcn,:)=0;
        end

        vnode(isnan(vnode))=0;
        vnode(isinf(vnode))=0;
        anode(isnan(anode))=0;
        anode(isinf(anode))=0;

        vpart=vpart+sff'*dt*anode;
        xp=xp+sff'*dt*vnode;

        nodemo=((partmass'.*(vpart))'*sff')';
        % Get updated velocities at the nodes
        vnode =( nodemo'./(nodmass))';

        % Impose the EBCs again
        if ebc==1
            vnode(ebcn,:)=0;
        end

        vnode(isnan(vnode))=0;
        vnode(isinf(vnode))=0;

        % Calculate stresses

        Fdot(1, 1,:)=vnode(:,1)'*dsfxx(:,:);
        Fdot(1, 2,:)=vnode(:,1)'*dsfyy(:,:);
        Fdot(2, 1,:)=vnode(:,2)'*dsfxx(:,:);
        Fdot(2, 2,:)=vnode(:,2)'*dsfyy(:,:);

        defgrad=defgrad+dt.*Fdot;
        velgrad=pagemtimes(Fdot,permute(defgrad, [2, 1, 3]));
        defrate=((velgrad + permute(velgrad, [2, 1, 3]))/2)*dt;
        J = arrayfun(@(k) det(defgrad(:,:,k)), 1:npart);
        volpart=volpart0.*J;
        strain=[defrate(1,1,:) defrate(2,2,:) 2*defrate(1,2,:)];
        stresss=pagemtimes(C,permute(strain,[2 1 3]));
        sigma=sigma+[stresss(1,1,:) stresss(3,1,:);stresss(3,1,:) stresss(2,1,:)];
        for n=1:npart
            stress(:,:,n)=sigma(:,:,n)*inv(defgrad(:,:,n)');
        end

        % Store data
        ptr(:,:,nstep)=ptraction(:,:);
        xxp(nstep,:,:)=xp(:,:);

        for n=1:npart
            S=[strs(1,1,nstep,n),strs(2,2,nstep,n),strs(1,2,nstep,n)];
            Ken(nstep) = Ken(nstep) + 0.5*partmass(n)*(vpart(n,1).^2+vpart(n,2).^2);
            Sen(nstep) = Sen(nstep) + (0.5) * (volpart(n)*(S*S')/E);
        end
        t=t+dt;

        % Store results
        tip_displacement(nstep) = mean(xp(end-10:end, 2)); % Average tip displacement

        if ~exist('all_xxp', 'var')
            all_xxp = zeros(length(excitation_types), Tstep, npart, 2);
            all_strs = zeros(length(excitation_types), Tstep, 2, 2, npart);
        end
        all_xxp(exc_idx, nstep, :, :) = xp;
        all_strs(exc_idx, nstep, :, :, :) = stress;

    end

    % Store results for this excitation type
    results.(excitation_type).time = time;
    results.(excitation_type).excitation = excitation_signal;
    results.(excitation_type).displacement = tip_displacement;
    results.(excitation_type).Ken = Ken;
    results.(excitation_type).Sen = Sen;

    fprintf('Completed: %s\n', excitation_type);
end

% Save results
save('excitation_comparison_results.mat', 'results');

%% Analysis and Plotting
close all;clc
fprintf('Generating analysis plots...\n');

% Define equations for each excitation type
equations = {
    'cosine: $f(t) = g \cdot \cos(2\pi f_0 t)$',...
    'sine: $f(t) = g \cdot \sin(2\pi f_0 t)$',...
    'Switching: $f(t) = \pm g$ (switching sign)',...
};

% Figure 1: Excitation Signals Comparison
figure('Position', [100, 100, 1200, 1000], 'Color', 'w');
BBB = ceil(length(excitation_types)/2);
ax_all = gobjects(length(excitation_types),1); % store subplot handles

for i = 1:length(excitation_types)
    ax = subplot(BBB, 2, i);
    ax_all(i) = ax; % store axis handle
    
    exc_type = excitation_types{i};
    time_data = results.(exc_type).time;
    excitation_data = results.(exc_type).excitation;
    
    plot(time_data, excitation_data, 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title([exc_type ' Excitation'], 'Interpreter', 'none');
    grid on;
    set(gca, 'FontSize', 10);
    xlim([0, min(20, T)]);
    
    % Add equation as text
    eq_idx = find(strcmp(exc_type, {'cosine', 'sine', 'Switching'}));
    if ~isempty(eq_idx) && eq_idx <= length(equations)
        text(0.98, 0.95, equations{eq_idx}, 'Units', 'normalized', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', ...
            'BackgroundColor', [1, 1, 0.8], 'FontSize', 9, 'Interpreter', 'latex');
    end
end

% --- Align Y labels across subplots ---
yLabelPositions = arrayfun(@(a) a.YLabel.Position(1), ax_all);
minPos = min(yLabelPositions); % farthest left position
for i = 1:length(ax_all)
    ax_all(i).YLabel.Position(1) = minPos;
end

saveas(gcf, 'excitation_signals_comparison.png');




% Figure 2: Displacement Response Comparison
figure('Position', [100, 100, 1200, 1000], 'Color', 'w');

colors = lines(length(excitation_types));
nRows  = BBB; 
nCols  = 2; 
ax_all = gobjects(length(excitation_types),1);

for i = 1:length(excitation_types)
    ax = subplot(nRows, nCols, i);
    ax_all(i) = ax;

    exc_type = excitation_types{i};
    plot(results.(exc_type).time, ...
         results.(exc_type).displacement / results.(exc_type).displacement(1,1), ...
         'LineWidth', 1.5, 'DisplayName', exc_type, 'Color', colors(i,:));

    xlim([0, min(20, T)]);
    title(exc_type, 'FontSize', 16, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 12);

    % Format tick labels to 3 decimals max
    xtickformat('%.3f');
    ytickformat('%.3f');

    % Show y-labels only in the left column
    if mod(i-1, nCols) == 0
        ylabel('Displacement (m)');
    else
        ylabel('');
        ax.YAxis.Visible = 'off';
    end

    % Show x-labels only in the bottom row
    if i > (nRows-1)*nCols
        xlabel('Time (s)');
    else
        xlabel('');
        ax.XAxis.Visible = 'off';
    end
end

% Align y-labels in left column
yLabelPositions = arrayfun(@(a) a.YLabel.Position(1), ax_all(mod((1:length(ax_all))-1,nCols)==0));
minPos = min(yLabelPositions);
for i = find(mod((1:length(ax_all))-1,nCols)==0)
    ax_all(i).YLabel.Position(1) = minPos;
end


saveas(gcf, 'displacement_response_comparison.png');



%% 1. Natural Frequency Analysis using FFT
figure('Position', [100, 100, 1200, 800], 'Color', 'w');
% sgtitle('Natural Frequency Analysis using FFT', 'FontSize', 16, 'FontWeight', 'bold');

natural_frequencies = zeros(length(excitation_types), 3); % Store first 3 natural frequencies
ax_all = gobjects(length(excitation_types),1); % store subplot handles

for i = 1:length(excitation_types)
    ax = subplot(BBB, 2, i);
    ax_all(i) = ax; % store axis handle
    exc_type = excitation_types{i};
    time = results.(exc_type).time;
    displacement = results.(exc_type).displacement;
    
    % Remove any NaN values
    valid_idx = ~isnan(displacement);
    time = time(valid_idx);
    displacement = displacement(valid_idx);
    
    % Detrend the signal
    displacement_detrended = detrend(displacement);
    
    % FFT analysis
    fs = 1/(time(2)-time(1)); % Sampling frequency
    N = length(displacement_detrended);
    f = (0:N-1)*(fs/N);
    
    Y = fft(displacement_detrended);
    P2 = abs(Y/N);
    P1 = P2(1:floor(N/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f_plot = f(1:floor(N/2)+1);
    
    % Find peaks in frequency domain (natural frequencies)
    [pks, locs] = findpeaks(P1, f_plot, 'MinPeakHeight', max(P1)*0.2, 'MinPeakDistance', 0.1);
    
    % Store first 3 natural frequencies
    num_freqs = min(3, length(pks));
    natural_frequencies(i, 1:num_freqs) = locs(1:num_freqs);
    
    % Plot frequency spectrum
    subplot(BBB,2, i);
    plot(f_plot, P1, 'LineWidth', 1.5, 'Color', colors(i,:));
    hold on;
    plot(locs, pks, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    

       % Show y-labels only in the left column
    if mod(i-1, nCols) == 0
        ylabel('Amplitude');
    else
        ylabel('');
        ax.YAxis.Visible = 'off';
    end

    % Show x-labels only in the bottom row
    if i > (nRows-1)*nCols
        xlabel('Frequency (Hz)');
    else
        xlabel('');
        ax.XAxis.Visible = 'off';
    end
    title(sprintf(exc_type, locs(1)), 'Interpreter', 'none');
    xlim([0 10]);
    grid on;
    
    % Annotate peaks
    for j = 1:min(3, length(pks))
        text(locs(j)+0.1, pks(j), sprintf('%.2f Hz', locs(j)), ...
            'VerticalAlignment', 'bottom', 'FontSize', 8);
    end
end
% --- Align Y labels across subplots ---
yLabelPositions = arrayfun(@(a) a.YLabel.Position(1), ax_all);
minPos = min(yLabelPositions); % farthest left position
for i = 1:length(ax_all)
    ax_all(i).YLabel.Position(1) = minPos;
end

saveas(gcf, 'natural_frequencies_analysis.png');

%% 2. Natural Frequency Comparison Table
figure('Position', [100, 100, 800, 300], 'Color', 'w');
natural_freq_table = cell(length(excitation_types), 4);
for i = 1:length(excitation_types)
    exc_type = excitation_types{i};
    freqs = natural_frequencies(i, :);
    natural_freq_table{i, 1} = exc_type;
    natural_freq_table{i, 2} = sprintf('%.3f', freqs(1));
    natural_freq_table{i, 3} = sprintf('%.3f', freqs(2));
    natural_freq_table{i, 4} = sprintf('%.3f', freqs(3));
end

uitable('Data', natural_freq_table(:, 2:4), ...
    'ColumnName', {'1st Nat Freq (Hz)', '2nd Nat Freq (Hz)', '3rd Nat Freq (Hz)'}, ...
    'RowName', excitation_types, ...
    'Position', [20, 20, 760, 250], ...
    'FontSize', 10);
title('Natural Frequencies Comparison');




%% Helper function for spider plot (if not available in your MATLAB version)
function spider_plot = spider_plot(P, options)
    % Simplified spider plot implementation
    figure;
    axes = gca;
    hold on;
    
    num_vars = size(P, 1);
    num_cases = size(P, 2);
    
    % Create spider web
    theta = linspace(0, 2*pi, num_vars+1);
    max_val = max(P(:)) * 1.2;
    
    for i = 1:num_vars
        plot([0 cos(theta(i))*max_val], [0 sin(theta(i))*max_val], 'k--', 'LineWidth', 0.5);
    end
    
    % Plot each case
    colors = lines(num_cases);
    for j = 1:num_cases
        rho = [P(:,j); P(1,j)];
        polarplot(theta, rho, 'Color', colors(j,:), 'LineWidth', 2);
    end
    
    spider_plot = gcf;
end

% Helper function
function [sf, dsfx, dsfy] = shapef(x, y, xn, yn, dxn, dyn)
sf = max(0, (dxn - abs(x-xn))/dxn) .* max(0, (dyn - abs(y-yn))/dyn);
dsfx = max(0, 1 - abs(y-yn)) .* -sign(x-xn) .* (y >= (yn-1)) .* (y <= (yn+1));
dsfy = max(0, 1 - abs(x-xn)) .* -sign(y-yn) .* (x >= (xn-1)) .* (x <= (xn+1));
end