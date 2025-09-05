% TL-MPM (Total Lagrangian Material Point Method) 2D Code
% ======================================================
%
% Description:
% ------------
% This code implements the Total Lagrangian Material Point Method (TL-MPM) 
% for 2D solid mechanics simulations. The method combines Lagrangian particle
% tracking with Eulerian grid-based calculations for large deformation problems.
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
%   The present work is prepared by Ibrahim Bouzaid, introduces modifications 
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

% MATERIAL PROPERTIES
% ====================
E = 2e7;          % Young's modulus (Pa)
materialDensity = 12000;       % Material density (kg/m³)
poissonsRatio = 0.3;          % Poisson's ratio

% Elasticity matrix for plane stress
c11 = E/(1-poissonsRatio^2);
c22 = c11;
c12 = c11*poissonsRatio;
c21 = c12;
c66 = E/(2*(1+poissonsRatio));

elasticityMatrix = [c11, c12, 0;
                    c12, c22, 0;
                    0,   0,   c66];

% INITIALIZATION
% ==============

% Eulerian Mesh Definition
% Domain dimensions
domainLengthX = 12;    % Length of domain in x-direction
domainLengthY = 3;     % Length of domain in y-direction
gridSpacingX = 1;      % Grid spacing in x-direction
gridSpacingY = 1;      % Grid spacing in y-direction

% Create nodal coordinates
X_coordinates = 0:gridSpacingX:domainLengthX;
Y_coordinates = 0:gridSpacingY:domainLengthY;
[Y_mesh, X_mesh] = meshgrid(Y_coordinates, X_coordinates);

% Reshape into column vectors
nodalY = reshape(Y_mesh, length(Y_coordinates)*length(X_coordinates), 1);
nodalX = reshape(X_mesh, length(X_coordinates)*length(Y_coordinates), 1);
nodalCoordinates = [nodalX, nodalY];

% Mesh properties
numberOfNodesX = length(X_coordinates);
numberOfNodesY = length(Y_coordinates);
numberOfNodes = length(X_coordinates) * length(Y_coordinates);
numberOfElements = (domainLengthX/gridSpacingX) * (domainLengthY/gridSpacingY);

% Particle Mesh Definition
% Solid body dimensions
Width = 10;        % Width of the solid body
Height = 1;       % Height of the solid body

% Particles per cell
particlesPerCellX = 2; % Number of particles per cell in x-direction
particlesPerCellY = 2; % Number of particles per cell in y-direction

% Calculate total number of particles
totalParticlesX = round(Width * particlesPerCellX / gridSpacingX);
totalParticlesY = round(Height * particlesPerCellY / gridSpacingX);

% Calculate particle spacing
particleSpacingX = Width / totalParticlesX;
particleSpacingY = Height / totalParticlesY;

% Initialize particle centers array
particleCenters = zeros(totalParticlesX * totalParticlesY, 2);

% Generate particle centers
particleIndex = 1;
for i = 1:totalParticlesX
    for j = 1:totalParticlesY
        % Calculate center coordinates for each particle
        centerX = (i - 0.5) * particleSpacingX;
        centerY = (j - 0.5) * particleSpacingY;
        particleCenters(particleIndex, :) = [centerX, centerY];
        particleIndex = particleIndex + 1;
    end
end

% Set initial particle positions (with vertical offset)
particlePositions(:, 1) = particleCenters(:, 1) + 0;   % X positions
particlePositions(:, 2) = particleCenters(:, 2) + 1;   % Y positions (offset by 1 unit)
numberOfParticles=length(particlePositions);

% Initialize particle geometry vectors (for calculating volume)
for iParticle = 1:numberOfParticles
    geometryVector1(iParticle,:) = [(Width/particlesPerCellX)/2, 0];
    geometryVector2(iParticle,:) = [0, (Height/particlesPerCellY)/2];
    deformationGradientSub{iParticle} = [1, 0; 0, 1]; % Initial deformation gradient
end

% Store initial geometry vectors for reference
geometryVector10 = geometryVector1;
geometryVector20 = geometryVector2;

% Calculate initial particle volumes
particleVolumes = zeros(1, numberOfParticles);
for iParticle = 1:numberOfParticles
    % Volume = 4 * |r1 × r2| (area of parallelogram)
    particleVolumes(iParticle) = 4 * abs(geometryVector1(iParticle,1) * geometryVector2(iParticle,2) - ...
                                     geometryVector1(iParticle,2) * geometryVector2(iParticle,1));
end

% Calculate particle masses
particleMasses = materialDensity * particleVolumes;

% VARIABLE INITIALIZATION
% ========================
initialParticleVolumes = particleVolumes;
cauchyStress = ones(2, 2, numberOfParticles);      % Cauchy stress tensor
strain = zeros(2, 2, numberOfParticles);           % Strain tensor
deformationGradient = eye(2) .* ones(2, 2, numberOfParticles); % Deformation gradient
velocityGradient = zeros(2, 2, numberOfParticles); % Velocity gradient
particleVelocities = zeros(numberOfParticles, 2);  % Particle velocities
sf = zeros(numberOfParticles, numberOfNodes);      % Shape functions
totalNodalForces = zeros(numberOfNodes, 2);        % Total nodal forces
secondPiolaKirchhoffStress = zeros(2, 2, numberOfParticles);   % second Piola Kirchhoff Stress
sigma=zeros(2,2,numberOfParticles);

% BOUNDARY CONDITIONS
% ===================
% Essential boundary conditions (fixed nodes)
essentialBoundaryConditionFlag = 1;
essentialBoundaryNodes = 1:numberOfNodesX:numberOfNodes;

% MESH VISUALIZATION AND SETUP ANALYSIS
% ======================================
% This section creates a visual of the initial Eulerian mesh
% and Lagrangian particle configuration with proper labeling
% and annotation for analysis and documentation purposes.

figure('Name', 'TL-MPM Mesh Configuration', 'Units', 'normalized', ...
       'Position', [0.1, 0.1, 0.8, 0.8], 'Color', 'white');

% Create main plot area
hold on;
grid on;
box on;
axis equal;

% Plot Eulerian mesh nodes
meshNodes = scatter(nodalCoordinates(:, 1), nodalCoordinates(:, 2), ...
                   80, 'o', 'MarkerEdgeColor', [0.2, 0.2, 0.6], ...
                   'MarkerFaceColor', [0.7, 0.7, 0.9], ...
                   'LineWidth', 1.5);

% Plot Lagrangian particles
particlePoints = scatter(particlePositions(:, 1), particlePositions(:, 2), ...
                        100, '.', 'MarkerEdgeColor', [0.6, 0.2, 0.2], ...
                        'MarkerFaceColor', [0.9, 0.7, 0.7], ...
                        'LineWidth', 1.5);

% Highlight essential boundary condition nodes
boundaryNodes = scatter(nodalCoordinates(essentialBoundaryNodes, 1), ...
                       nodalCoordinates(essentialBoundaryNodes, 2), ...
                       120, '^', 'MarkerEdgeColor', [0.8, 0.1, 0.1], ...
                       'MarkerFaceColor', [1.0, 0.6, 0.6], ...
                       'LineWidth', 2.5);

% Add labels for selected nodes and particles (avoid clutter)
labelSpacing = 0.1; % Distance offset for labels
maxLabels = min(8, numberOfNodes); % Limit number of labels to avoid clutter

% Label every nth node for clarity
for i = 1:ceil(numberOfNodes/maxLabels):numberOfNodes
    text(nodalCoordinates(i, 1) + labelSpacing, nodalCoordinates(i, 2) + labelSpacing, ...
         sprintf('N%d', i), 'FontSize', 8, 'FontWeight', 'bold', ...
         'Color', [0.2, 0.2, 0.6], 'BackgroundColor', [1, 1, 1, 0.7]);
end

% Label selected particles
maxParticleLabels = min(5, numberOfParticles);
for i = 1:ceil(numberOfParticles/maxParticleLabels):numberOfParticles
    text(particlePositions(i, 1) + labelSpacing, particlePositions(i, 2) + labelSpacing, ...
         sprintf('P%d', i), 'FontSize', 8, 'FontWeight', 'bold', ...
         'Color', [0.6, 0.2, 0.2], 'BackgroundColor', [1, 1, 1, 0.7]);
end

% Highlight boundary condition nodes with special labels
for i = 1:length(essentialBoundaryNodes)
    nodeIdx = essentialBoundaryNodes(i);
    text(nodalCoordinates(nodeIdx, 1) + labelSpacing, nodalCoordinates(nodeIdx, 2) - labelSpacing, ...
         sprintf('BC%d', i), 'FontSize', 9, 'FontWeight', 'bold', ...
         'Color', [0.8, 0.1, 0.1], 'BackgroundColor', [1, 1, 1, 0.8]);
end

% Set plot properties
xlim([-1, domainLengthX + 1]);
ylim([-1, domainLengthY + 1]);
xlabel('X Position (m)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Y Position (m)', 'FontSize', 12, 'FontWeight', 'bold');
title('TL-MPM Initial Mesh Configuration', 'FontSize', 14, 'FontWeight', 'bold');

% Add legend
legend([meshNodes, particlePoints, boundaryNodes], ...
       {'Eulerian Nodes', 'Lagrangian Particles', 'Boundary Nodes'}, ...
       'Location', 'northeastoutside', 'FontSize', 10);

% Add grid information annotation
gridInfoText = {sprintf('Eulerian Mesh: %dx%d nodes', numberOfNodesX, numberOfNodesY), ...
               sprintf('Grid spacing: %.2f x %.2f m', gridSpacingX, gridSpacingY), ...
               sprintf('Domain: %.1f x %.1f m', domainLengthX, domainLengthY), ...
               sprintf('Particles: %d', numberOfParticles), ...
               sprintf('Particles/cell: %dx%d', particlesPerCellX, particlesPerCellY), ...
               sprintf('Boundary nodes: %d', length(essentialBoundaryNodes))};

annotation('textbox', [0.4, 0, 0.25, 0.2], 'String', gridInfoText, ...
           'FitBoxToText', 'on', 'BackgroundColor', [1, 1, 1, 0.8], ...
           'EdgeColor', [0.5, 0.5, 0.5], 'FontSize', 9);

% SAVE MESH VISUALIZATION
% ========================
% saveas(gcf, 'TL_MPM_Mesh_Configuration.png');
fprintf('Mesh configuration plot saved as: TL_MPM_Mesh_Configuration.png\n');

% Display mesh information in command window
fprintf('\n=== MESH CONFIGURATION SUMMARY ===\n');
fprintf('Eulerian Mesh:\n');
fprintf('  - Domain size: %.1f x %.1f m\n', domainLengthX, domainLengthY);
fprintf('  - Grid spacing: %.2f x %.2f m\n', gridSpacingX, gridSpacingY);
fprintf('  - Number of nodes: %d (%dx%d)\n', numberOfNodes, numberOfNodesX, numberOfNodesY);
fprintf('  - Number of elements: %d\n', numberOfElements);
fprintf('\nLagrangian Particles:\n');
fprintf('  - Number of particles: %d\n', numberOfParticles);
fprintf('  - Particles per cell: %dx%d\n', particlesPerCellX, particlesPerCellY);
fprintf('  - Particle spacing: %.3f x %.3f m\n', particleSpacingX, particleSpacingY);
fprintf('  - Solid dimensions: %.1f x %.1f m\n', Width, Height);
fprintf('\nBoundary Conditions:\n');
fprintf('  - Fixed boundary nodes: %d\n', length(essentialBoundaryNodes));
fprintf('  - Boundary location: Left side (x=0)\n');

% TIME INTEGRATION SETUP
% ======================
currentTime = 0.0;
timeStepCounter = 0;

% CFL condition: Δt <= C * h_min / c
%   h_min : smallest grid spacing
%   c     : wave speed (sqrt(E/ρ) for elastic material)
%   C     : Courant number (≤ 1, often 0.4–0.8 for stability)

h_min   = min(gridSpacingX,gridSpacingY);       % replace with your grid spacing
c       = sqrt(E/materialDensity);              % material wave speed
C       = 0.1;                                  % Courant safety factor
dt = C * h_min / c;                             % stable time step

totalTime = 50;                                  % Choose the total time
totalTimeSteps = round(totalTime/dt);           % Number of time steps

% Initialize traction variables
particleTraction = zeros(numberOfParticles, 2);
nodalTraction = zeros(numberOfNodes, 2);
% Particle Positions over the total time
XP_History = zeros(numberOfParticles,2,totalTimeSteps);

% SHAPE FUNCTION CALCULATION
% ==========================
[sf, dsfx, dsfy] = shapef(particlePositions(:,1),particlePositions(:,2), nodalCoordinates(:,1),nodalCoordinates(:,2), gridSpacingX, gridSpacingY);
dsfx(isnan(dsfx(:,:)))=0;
dsfy(isnan(dsfy(:,:)))=0;
dsf(1,:,:)=dsfx;
dsf(2,:,:)=dsfy;
dsfxx(:,:)=dsf(1,:,:);
dsfyy(:,:)=dsf(2,:,:);
% Calculate nodal mass (mass mapping from particles to nodes)
% ===========================================================
nodalMass = particleMasses * sf';

% LOADING PROFILE PARAMETERS
% ==========================
initialSwitchInterval = 50;      % Initial number of time steps between load reversals
switchIntervalIncrement = 50;    % Increment to increase reversal interval each cycle
gravitationalAcceleration = 9.81; % Gravitational acceleration (m/s²)
tractionMagnitude = 10000;           % Magnitude of applied traction (Pa)
switchingFunction = zeros(1, totalTimeSteps);  % Load direction function
currentLoadDirection = -1;                     % Initial load direction (-1: downward, +1: upward)
cyclesCompleted = 0;                           % Number of completed load cycles
stepsSinceLastSwitch = 0;                      % Counter since last load reversal
currentSwitchInterval = initialSwitchInterval; % Current interval between reversals
NeumannBCs=1:2:8;                              % Nodes subject to Neumann boundary conditions


% MAIN SIMULATION LOOP
% ====================
tic
for timeStep = 1:totalTimeSteps
    fprintf('Time step: %d/%d\n', timeStep, totalTimeSteps);
    
    % Reset Eulerian mesh variables
    nodalMomentum(:) = 0;
    nodalVelocities(:) = 0;
    velocityGradient(:) = 0;
    externalNodalForces(:) = 0;
    internalNodalForces(:) = 0;
    totalNodalForces(:) = 0;
    
    % Extract stress components for vectorized operations
    stressComponent11(:,:) = cauchyStress(1,1,:);
    stressComponent22(:,:) = cauchyStress(2,2,:);
    stressComponent12(:,:) = cauchyStress(1,2,:);

    % DYNAMIC LOADING PROFILE UPDATE
    % ===============================
    % Update the switching function for frequency-modulated loading
    stepsSinceLastSwitch = stepsSinceLastSwitch + 1;
    if stepsSinceLastSwitch >= currentSwitchInterval
        % Reverse load direction
        currentLoadDirection = -currentLoadDirection;
        stepsSinceLastSwitch = 0;
        cyclesCompleted = cyclesCompleted + 0.5; % Half-cycle completed
        currentSwitchInterval = currentSwitchInterval + switchIntervalIncrement;
    end
    particleTraction(NeumannBCs,2) = currentLoadDirection * tractionMagnitude;
    ptr(timeStep,:,:)=particleTraction;


    % FORCE CALCULATION
    % =================
    % External force rate
    externalForceRate = [0, 0]; % Gravity or other body forces
    
    % Calculate nodal momentum
    nodalMomentum = (particleMasses' .* particleVelocities)' * sf';

    % External nodal forces
    externalNodalForces = (externalForceRate' * nodalMass)';
    
    % Internal nodal forces (stress divergence)
    internalNodalForcesX = ((-particleVolumes' .* stressComponent11)' * dsfx' + ...
                           (-particleVolumes' .* stressComponent12)' * dsfy');
    internalNodalForcesY = ((-particleVolumes' .* stressComponent12)' * dsfx' + ...
                           (-particleVolumes' .* stressComponent22)' * dsfy');
    internalNodalForces = [internalNodalForcesX; internalNodalForcesY]';
    
    % Traction forces
    nodalTraction = ((particleVolumes' .* particleTraction)' * sf')';
    
    % Total nodal forces
    totalNodalForces = externalNodalForces + internalNodalForces + nodalTraction;
    
    % Update nodal momentum (impulse-momentum theorem)
    nodalMomentum = (nodalMomentum + dt * totalNodalForces')';
    
    % UPDATE PARTICLE VELOCITIES AND POSITIONS
    % ========================================
    % Nodal acceleration
    nodalAcceleration = ((totalNodalForces)' ./ nodalMass)';
    
    % Nodal velocities
    nodalVelocities = ((nodalMomentum)' ./ nodalMass)';
    
    % Apply essential boundary conditions (fixed nodes)
    nodalVelocities(essentialBoundaryNodes,:) = 0;
    nodalAcceleration(essentialBoundaryNodes,:) = 0;
    
    % Handle numerical issues
    nodalVelocities(isnan(nodalVelocities)) = 0;
    nodalVelocities(isinf(nodalVelocities)) = 0;
    nodalAcceleration(isnan(nodalAcceleration)) = 0;
    nodalAcceleration(isinf(nodalAcceleration)) = 0;
    
    % Update particle velocities using FLIP method
    particleVelocities = particleVelocities + sf' * dt * nodalAcceleration;
    
    % Update particle positions
    particlePositions = particlePositions + sf' * dt * nodalVelocities;
    
    % UPDATE NODAL QUANTITIES
    % ========================
    % Recalculate nodal mass and momentum after particle update
    nodalMomentum = ((particleMasses' .* particleVelocities)' * sf')';
    
    % Update nodal velocities
    nodalVelocities = (nodalMomentum' ./ nodalMass)';
    
    % Reapply essential boundary conditions
    nodalVelocities(essentialBoundaryNodes,:) = 0;
    
    % Handle numerical issues
    nodalVelocities(isnan(nodalVelocities)) = 0;

    % STRESS UPDATE AND DEFORMATION CALCULATION
    % =========================================
    for n = 1:numberOfParticles
        Fdot = zeros(2,2);
        for i = find(any(dsf(:,:,n) ~= 0, 1))
            Fdot = Fdot + nodalVelocities(i,:)' * dsf(:,i,n)';
        end
        deformationGradient(:,:,n) = deformationGradient(:,:,n) + dt * Fdot;
        velocityGradient = Fdot * deformationGradient(:,:,n)';
        strain = 0.5 * (velocityGradient + velocityGradient') * dt;
        J = det(deformationGradient(:,:,n));
        particleVolumes(n) = initialParticleVolumes(n) * J;
        sigma(:,:,n) = sigma(:,:,n) + [c11*strain(1,1) + c12*strain(2,2), c66*2*strain(1,2);
                                        c66*2*strain(1,2), c22*strain(2,2) + c21*strain(1,1)];
        secondPiolaKirchhoffStress(:,:,n) = sigma(:,:,n) * inv(deformationGradient(:,:,n))';
        strs(timeStep,n,:,:) = secondPiolaKirchhoffStress(:,:,n);

    end

    % Save Particle positions for figures
    XP_History(:,:,timeStep)=particlePositions(:, :);

    % Update time
    currentTime = currentTime + timeStep;
end
toc


% OUTPUT SIMULATION INFORMATION
fprintf('Simulation completed:\n');
fprintf('Number of particles: %d\n', numberOfParticles);
fprintf('Number of nodes: %d\n', numberOfNodes);



%% VISUALIZATION AND VIDEO CREATION
% ================================

% Create figure for visualization
figure
set(gcf, 'units', 'inches', 'position', [3, 1, 8, 8]);
set(gcf, 'Color', 'white'); % White background for better visualization

% Initialize video writer
video = VideoWriter('TL_MPM_2D_Simulation.avi');
video.FrameRate = 15; % Set frame rate
open(video);

% Determine visualization steps (show every 10th step)
visualizationStep = 200;

fprintf('Creating visualization video...\n');

% Create storage for particle positions and stresses over time
% (Assuming these were stored during the simulation)
% If not stored, you would need to modify the main loop to save this data
% particle positions over time [timeSteps × npart × 2]


for timeIdx = 1:visualizationStep:totalTimeSteps
    % Clear previous frame
    clf;
    
    % Get current particle positions and stresses
    % Plot particles colored by stress
    scatter(XP_History(:, 1,timeIdx), XP_History(:, 2,timeIdx),'.',LineWidth=3);
    
    % Add grid and labels
    grid on;
    hold on;
    
    % Plot Eulerian mesh nodes for reference
    scatter(nodalCoordinates(:, 1), nodalCoordinates(:, 2), 20, 'k', 'filled', 'MarkerFaceAlpha', 0.3);
    
    title(sprintf('TL-MPM Simulation - Time: %.3f s (Step %d/%d)', ...
          timeIdx * dt, timeIdx, totalTimeSteps));
    xlabel('X Position (m)');
    ylabel('Y Position (m)');
    
    % Set axis limits
    axis equal;
    xlim([-1, domainLengthX + 1]);
    ylim([-1, domainLengthY + 1]);
    
    % Add simulation info text
    simulationInfo = {sprintf('Time: %.3f s', timeIdx * dt), ...
                     sprintf('Step: %d/%d', timeIdx, totalTimeSteps), ...
                     sprintf('Particles: %d', numberOfParticles), ...
                     sprintf('Nodes: %d', numberOfNodes)};
    annotation('textbox', [0.75, 0.75, 0.2, 0.15], 'String', simulationInfo, ...
               'FitBoxToText', 'on', 'BackgroundColor', 'white');
    
    % Capture frame for video
    frame = getframe(gcf);
    writeVideo(video, frame);
end

% Close video file
close(video);
fprintf('Video saved as: TL_MPM_2D_Simulation.avi\n');

%% ADDITIONAL POST-PROCESSING
% ==========================

figure('Position', [100, 100, 1200, 800]);

% Final configuration plot
scatter(XP_History(:, 1,end), XP_History(:, 2,end));
colorbar;
colormap('jet');
title('Final Configuration - Stress Distribution');
xlabel('X Position (m)');
ylabel('Y Position (m)');
axis equal;
grid on;


%%
plot(1:timeStep, ptr(:,1,2),'-')

function [sff, dsfx, dsfy] = shapef(x, y, xn, yn, dxn, dyn)  
% function [shapeFunctions, shapeFunctionDerivativesX, shapeFunctionDerivativesY] = ...
%          calculateShapeFunctions(particleX, particleY, nodeX, nodeY, gridSpacingX, gridSpacingY)
%     - Calculate the shape functions and derivatives
%     - Inputs: particle coordinates, node coordinates, grid spacing
%     - Outputs: shape functions and their derivatives
sff = max(0, (dxn - abs(xn - x.'))/dxn).* max(0, (dyn - abs(yn - y.'))/dyn);
% Calculate derivatives with respect to x
dsfx = max(0, 1 - abs(y.'-yn)).*-sign(x.'-xn).*double(y.' >= abs(yn-1)).*double(y.' <= yn+1).*sff;
% Calculate derivatives with respect to y
dsfy = max(0, 1 - abs(x.'-xn)).*-sign(y.'-yn).*double(x.' >= abs(xn-1)).*double(x.' <= xn+1).*sff;
% Normalize derivatives
dsfx=dsfx./sff;
dsfy=dsfy./sff;
end








