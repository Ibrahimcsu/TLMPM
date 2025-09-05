% BSMPM-2D: B-Spline Material Point Method for 2D Solid Mechanics
% ===============================================================
%
% Description:
% ------------
%   This script implements the B-Spline Material Point Method (BSMPM) for 
%   simulating solid mechanics problems in 2D. The method uses B-spline 
%   basis functions for improved accuracy and stability compared to standard 
%   MPM approaches, particularly for problems involving large deformations.
%
% Features:
% ---------
% - B-spline basis functions for enhanced accuracy
% - Linear elastic material model
% - Explicit time integration
%
% Key Variables:
% x, y        - Nodal coordinates
% xp          - Material point positions
% partmass    - Particle masses
% volpart     - Particle volumes
% defgrad     - Deformation gradient (2x2xMP)
% stress      - Cauchy stress at material points
% vpart       - Particle velocities
% velgrad     - Velocity gradient
% nodmass     - Nodal mass
% vnode       - Nodal velocity
% forcenodee  - External nodal forces
% forcenodei  - Internal nodal forces
% nodemo      - Nodal momentum
% sf          - Shape functions
% dsf         - Shape function derivatives
%
% Usage:
% ------
% 1. Define material properties (E, density, nu)
% 2. Set B-spline parameters (degree, domain size)
% 3. Create particle mesh and Eulerian grid
% 4. Set initial conditions and boundary conditions
% 5. Run simulation with appropriate time step
%
% Particle Mesh Creation:
% ----------------------
% xp = array of particle coordinates [npart × 2]
% Create particles within a rectangular domain with specified
% particles per cell in x and y directions
%
% Eulerian Mesh Creation:
% ----------------------
% LOC = nodal coordinates [nnodes × 2]
% Create grid with specified number of segments in x and y directions
%
% Authorship and Contributions:
%   This implementation builds upon the original BSMPM methodology with
%   modifications to the mapping procedure and shape function evaluations,
%   incorporating a fully vectorized structure to improve computational
%   efficiency.
%
% Reference:
% -----------
% Steffen, M., Wallstedt, P. C., Guilkey, J. E., Kirby, R. M., & Berzins, M. (2008).
% "Examination and comparison of different spatial discretization methods."
% Journal of Computational Physics, 227(24), 10059-10086.
%
% Version: 1.0
% Date: August 2025
%==========================================================================
%==========================================================================

clear all; close all; clc;
format long

% MATERIAL PROPERTIES AND SIMULATION PARAMETERS
E = 2e9;               % Young's modulus (Pa)
density = 12000;       % Material density (kg/m³)
nu = 0.3;              % Poisson's ratio

% Elasticity tensor for plane stress
c11 = E/(1-nu^2);
c22 = c11;
c12 = c11*nu;
c21 = c12;
c66 = E/(2*(1+nu));
C = [c11 c12 0;
     c12 c22 0;
     0   0   c66];

% B-spline basis function parameters
degree_x = 2;          % Degree of B-spline in x-direction
degree_y = 2;          % Degree of B-spline in y-direction
space_x = 6;           % Domain size in x-direction
space_y = 6;           % Domain size in y-direction

% Time parameters
t=0;
T = 0.3;               % Total simulation time
dt = 0.001;            % Time step size
Tstep = round(T/dt);   % Number of time steps

% External loading parameters
g = 30000;             % Gravity acceleration
initial_switch_interval = 200;  % Switching interval for oscillatory load
interval_increment = 0;         % Increment to increase switch interval

% B-SPLINE MESH GENERATION
% Generate B-spline basis functions and mesh
num_segments_x = 6;
num_segments_y = 6;
dxn = space_x/num_segments_x;
dyn = space_y/num_segments_y;

% Calculate number of basis functions
[nfunx, nfuny] = nBasis(space_x, num_segments_x, degree_x, space_y, num_segments_y, degree_y);
nnodes = nfunx * nfuny;

% Generate Eulerian mesh for visualization
LLx = space_x;
LLy = space_y;
X = 0:dxn:LLx;
Y = 0:dyn:LLy;
[x1, y1] = meshgrid(Y, X);
y = reshape(x1, length(Y)*length(X), 1);
x = reshape(y1, length(X)*length(Y), 1);

% MATERIAL POINTS INITIALIZATION
% Create a grid of material points
width = 6;
height = 1;
number_of_particlesx_per_cell = 5;
number_of_particlesy_per_cell = 5;

T_num_particles_x = round(width*number_of_particlesx_per_cell/dxn);
T_num_particles_y = round(height*number_of_particlesy_per_cell/dxn);
small_width = width / T_num_particles_x;
small_height = height / T_num_particles_y;

% Initialize particle centers
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

% Position particles with vertical offset
xp = zeros(size(centers));
xp(:, 1) = centers(:, 1) + 0;
xp(:, 2) = centers(:, 2) + 2;
npart = length(xp(:, 1)); % Number of particles

% Initialize particle volumes and masses
volpart = small_width * small_height * ones(1, npart);
partmass = density * volpart;
volpart0 = volpart;     % Store initial volumes
partmass0 = partmass;   % Store initial masses
traction=zeros(nnodes,2);

% Initialize particle velocities
vpart = zeros(npart, 2);

% INITIALIZE FIELD VARIABLES
% Stress and deformation variables
stress = zeros(2, 2, npart);     % Cauchy stress tensor
strain = zeros(2, 2, npart);     % Strain tensor
defgrad = repmat(eye(2), 1, 1, npart); % Deformation gradient
velgrad = zeros(2, 2, npart);    % Velocity gradient

% B-spline shape functions and derivatives
sfff = zeros(nfunx, nfuny, npart);
dsfxxx = zeros(nfunx, nfuny, npart);
dsfyyy = zeros(nfunx, nfuny, npart);

% Nodal variables
nodemo = zeros(nnodes, 2);       % Nodal momentum
nodmass = zeros(nnodes);      % Nodal mass
vnode = zeros(nnodes, 2);        % Nodal velocity
forcenodee = zeros(nnodes, 2);   % External nodal forces
forcenodei = zeros(nnodes, 2);   % Internal nodal forces
forcenodetotal = zeros(nnodes, 2); % Total nodal forces

% Boundary conditions
ebc = 1; % Essential boundary condition flag
if ebc == 1
    ebcn = 1:nfunx:length(x); % Fixed nodes
else
    ebcn = [];
end

% Energy tracking
Ken = zeros(Tstep, 1); % Kinetic energy
Sen = zeros(Tstep, 1); % Strain energy

% Data storage for analysis
strs = zeros(2, 2, Tstep, npart); % Stress history
xxp = zeros(Tstep, npart, 2);     % Position history
vp = zeros(Tstep, npart, 2);      % Velocity history
nomo = zeros(Tstep, nnodes, 2);   % Nodal momentum history
noms = zeros(Tstep, nnodes, 1);   % Nodal mass history
vn = zeros(Tstep, nnodes, 2);     % Nodal velocity history
fint = zeros(Tstep, nnodes, 2);   % Internal forces
fext = zeros(Tstep, nnodes, 2);   % External forces
trac = zeros(Tstep, nnodes, 2);   % Traction forces
velg = zeros(2, 2, Tstep, npart); % Velocity gradient history

% Particle traction
ptraction = zeros(npart, 2);
ptr = zeros(npart, 2, Tstep);     % Traction history

% Switching function for oscillatory loading
switching_function = zeros(1, Tstep);
current_switch_interval = initial_switch_interval;
current_value = -1;
switch_counter = 0;

% VISUALIZE INITIAL CONFIGURATION
figure
set(gcf, 'units', 'inches', 'position', [3 1 8 8])
hold on
plot(xp(:, 1), xp(:, 2), 'sq', 'Markersize', 5, 'LineWidth', 5)
if ebc == 1
    plot(x(ebcn), y(ebcn), '*', 'Markersize', 5, 'LineWidth', 5)
end
plot(x(:,:), y(:,:), '.', 'Markersize', 10, 'LineWidth', 3)
title('Initial Configuration: Material Points and Eulerian Mesh')
xlabel('X position')
ylabel('Y position')
legend('Material Points', 'Fixed Nodes', 'Eulerian Nodes', 'Location', 'best')
grid on

% MAIN SIMULATION LOOP
fprintf('Starting BSMPM simulation with %d time steps...\n', Tstep);
tic

for nstep = 1:Tstep
    fprintf('Time step: %d/%d\n', nstep, Tstep);
    t = t + dt;
    
    % Reset nodal quantities
    nodemo(:) = 0;
    nodmass(:) = 0;
    vnode(:) = 0;
    velgrad(:, :, :) = 0;
    forcenodee(:) = 0;
    forcenodei(:) = 0;
    forcenodetotal(:) = 0;
    
    % Compute B-spline shape functions and derivatives for all particles
    for p = 1:npart
        [sf, dsfx, dsfy] = BSMPM2d(space_x, space_y, num_segments_x, num_segments_y, ...
                                  degree_x, degree_y, xp(p, 1), xp(p, 2));
        sfff(:, :, p) = sf;
        dsfxxx(:, :, p) = dsfx;
        dsfyyy(:, :, p) = dsfy;
    end
    
    % Reshape for efficient computation
    sff = reshape(sfff, nnodes, npart);
    dsfxx = reshape(dsfxxx, nnodes, npart);
    dsfyy = reshape(dsfyyy, nnodes, npart);
    dsf = cat(3, dsfxx, dsfyy);
    dsf = permute(dsf, [3, 1, 2]);
    
    
    % Extract stress components for efficient computation
    strs11(:,:)=stress(1,1,:);
    strs22(:,:)=stress(2,2,:);
    strs12(:,:)=stress(1,2,:);
    dsfxx(:,:)=dsf(1,:,:);
    dsfyy(:,:)=dsf(2,:,:);

    
    % MAPPING: PARTICLE TO GRID (P2G)
    % Map mass and momentum to nodes
    nodmass = partmass * sff';
    nodemo = (partmass' .* vpart)' * sff';
    
    % Compute external forces (gravity)
    frate = [0, -10]; % Force rate
    forcenodee = (frate' * nodmass)';
    
    % Compute internal forces (stress divergence)
    forcenodeix = ((-volpart' .* strs11)' * dsfxx' + (-volpart' .* strs12)' * dsfyy');
    forcenodeiy = ((-volpart' .* strs12)' * dsfxx' + (-volpart' .* strs22)' * dsfyy');
    forcenodei = [forcenodeix; forcenodeiy]';
    
    % Compute traction forces
    traction = ((volpart' .* ptraction)' * sff')';
    
    % Total nodal forces
    forcenodetotal = forcenodee + forcenodei + traction;
    
    % Update nodal momentum (impulse-momentum theorem)
    nodemo = (nodemo + dt * forcenodetotal')';
    
    % Apply essential boundary conditions
    if ebc == 1
        nodemo(ebcn, :) = 0;
        forcenodetotal(ebcn, :) = 0;
    end
    
    % GRID UPDATE
    % Compute nodal acceleration and velocity
    anode = forcenodetotal ./ nodmass';
    vnode = nodemo ./ nodmass';
    
    % Apply essential boundary conditions
    if ebc == 1
        vnode(ebcn, :) = 0;
        anode(ebcn, :) = 0;
    end
    
    % Handle numerical issues
    vnode(isnan(vnode)) = 0;
    vnode(isinf(vnode)) = 0;
    anode(isnan(anode)) = 0;
    anode(isinf(anode)) = 0;
    
    % MAPPING: GRID TO PARTICLE (G2P)
    % Update particle velocities and positions
    vpart = vpart + sff' *dt* anode;
    xp = xp + sff' * dt*vnode;
    
    % Recompute nodal mass and momentum
    nodmass = partmass * sff';
    nodemo = (partmass' .* vpart)' * sff';
    
    % Update nodal velocities
    vnode = nodemo' ./ nodmass';
    
    % Apply essential boundary conditions
    if ebc == 1
        vnode(ebcn, :) = 0;
    end
    
    % Handle numerical issues
    vnode(isnan(vnode)) = 0;
    vnode(isinf(vnode)) = 0;
    
    % UPDATE DEFORMATION AND STRESS
    % Compute velocity gradient at particles
    velgrad(1, 1, :)=vnode(:,1)'*dsfxx;
    velgrad(1, 2, :)=vnode(:,1)'*dsfyy;
    velgrad(2, 1, :)=vnode(:,2)'*dsfxx;
    velgrad(2, 2, :)=vnode(:,2)'*dsfyy;
    defrate=((velgrad + permute(velgrad, [2, 1, 3]))/2)*dt;
    defgrad=(repmat(eye(2), 1, 1, npart)+velgrad*dt).*defgrad;
    J = arrayfun(@(k) det(defgrad(:,:,k)), 1:npart);
    volpart=volpart0.*J;
    strain=[defrate(1,1,:) defrate(2,2,:) 2*defrate(1,2,:)];
    stresss=pagemtimes(C,permute(strain,[2 1 3]));
    stress=stress+[stresss(1,1,:) stresss(3,1,:);stresss(3,1,:) stresss(2,1,:)];
    
    % DATA STORAGE AND ENERGY COMPUTATION
    % Store current state
    xxp(nstep, :, :) = xp;
    strs(:, :, nstep, :) = stress;
    
    % Compute energies
    for n = 1:npart
        S = [stress(1, 1, n), stress(2, 2, n), stress(1, 2, n)];
        Ken(nstep) = Ken(nstep) + 0.5 * partmass(n) * (vpart(n, 1)^2 + vpart(n, 2)^2);
        Sen(nstep) = Sen(nstep) + 0.5 * volpart(n) * (S * S') / E;
    end
end

toc
fprintf('BSMPM simulation completed successfully.\n');


%% Plot the particle's positions @ 50
scatter(xxp(50,:,1),xxp(50,:,2))


% B-SPLINE BASIS FUNCTION CALCULATION FUNCTIONS

function [nBasis_x, nBasis_y] = nBasis(space_x, num_segments_x, degree_x, space_y, num_segments_y, degree_y)
    % Calculate number of B-spline basis functions in each direction
    
    % Define the size of each segment
    segment_sizex = space_x / num_segments_x;
    segment_sizey = space_y / num_segments_y;
    
    % Create the internal knots based on the number of segments
    internal_knotsx = 0:segment_sizex:space_x;
    internal_knotsy = 0:segment_sizey:space_y;
    
    % Create the full knot vector with clamped ends
    knotsx = [zeros(1, degree_x), internal_knotsx, space_x * ones(1, degree_x)];
    knotsy = [zeros(1, degree_y), internal_knotsy, space_y * ones(1, degree_y)];
    
    % Calculate number of basis functions
    nBasis_x = length(knotsx) - degree_x - 1;
    nBasis_y = length(knotsy) - degree_y - 1;
end

function [sf, dx_values, dy_values] = BSMPM2d(space_x, space_y, num_segments_x, num_segments_y, degree_x, degree_y, X, Y)
    % Calculate 2D B-spline basis functions and their derivatives at point (X,Y)
    
    % Generate knot vectors
    [knots_x, knots_y] = generate_2D_knots(space_x, space_y, num_segments_x, num_segments_y, degree_x, degree_y);
    
    % Number of basis functions in each direction
    nBasis_x = length(knots_x) - degree_x - 1;
    nBasis_y = length(knots_y) - degree_y - 1;
    
    % Initialize matrices to store the 2D B-spline basis functions and derivatives
    sf = zeros(nBasis_x, nBasis_y);
    dx_values = zeros(nBasis_x, nBasis_y);
    dy_values = zeros(nBasis_x, nBasis_y);
    
    % Calculate and store the 2D B-spline basis functions and their derivatives
    for i = 1:nBasis_x
        for j = 1:nBasis_y
            % B-spline basis function
            sf(i, j) = bspline_basis(i, degree_x, knots_x, X) .* bspline_basis(j, degree_y, knots_y, Y);
            
            % Derivative in x direction
            dx_values(i, j) = bspline_derivative(i, degree_x, knots_x, X) .* bspline_basis(j, degree_y, knots_y, Y);
            
            % Derivative in y direction
            dy_values(i, j) = bspline_basis(i, degree_x, knots_x, X) .* bspline_derivative(j, degree_y, knots_y, Y);
        end
    end
end

function [knots_x, knots_y] = generate_2D_knots(space_x, space_y, num_segments_x, num_segments_y, degree_x, degree_y)
    % Generate knot vectors for both directions
    knots_x = generate_knot_vector(space_x, num_segments_x, degree_x);
    knots_y = generate_knot_vector(space_y, num_segments_y, degree_y);
end

function knots = generate_knot_vector(space, num_segments, degree)
    % Generate a clamped knot vector for B-splines
    
    % Define the size of each segment
    segment_size = space / num_segments;
    
    % Create the internal knots based on the number of segments
    internal_knots = 0:segment_size:space;
    
    % Create the full knot vector with clamped ends
    knots = [zeros(1, degree), internal_knots, space * ones(1, degree)];
end

function N = bspline_basis(i, p, knots, x)
    % Calculate B-spline basis function of degree p at point x
    
    if p == 0
        % Zeroth-degree B-spline (piecewise constant)
        N = double(knots(i) <= x & x < knots(i + 1));
    else
        % Recursive definition
        N1 = 0;
        if knots(i + p) ~= knots(i)
            N1 = ((x - knots(i)) / (knots(i + p) - knots(i))) .* bspline_basis(i, p - 1, knots, x);
        end
        
        N2 = 0;
        if knots(i + p + 1) ~= knots(i + 1)
            N2 = ((knots(i + p + 1) - x) / (knots(i + p + 1) - knots(i + 1))) .* bspline_basis(i + 1, p - 1, knots, x);
        end
        
        N = N1 + N2;
    end
end

function dN = bspline_derivative(i, p, knots, x)
    % Calculate derivative of B-spline basis function of degree p at point x
    
    if p == 0
        % Zeroth-degree B-spline derivative is zero (piecewise constant)
        dN = zeros(size(x));
    else
        % Derivative calculation
        dN1 = 0;
        if knots(i + p) ~= knots(i)
            dN1 = (p / (knots(i + p) - knots(i))) * bspline_basis(i, p - 1, knots, x);
        end
        
        dN2 = 0;
        if knots(i + p + 1) ~= knots(i + 1)
            dN2 = (p / (knots(i + p + 1) - knots(i + 1))) * bspline_basis(i + 1, p - 1, knots, x);
        end
        
        dN = dN1 - dN2;
    end
end