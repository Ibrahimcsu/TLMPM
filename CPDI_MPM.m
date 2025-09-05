% CPDI-MPM 2D Implementation for Solid Mechanics Simulation
% ======================================================
%
% Description:
% ------------
% This code implements the Convected Particle Domain Interpolation (CPDI)
% Material Point Method (MPM) for 2D solid mechanics problems.
%
% Features:
% ---------
% - Linear elastic material model
% - Explicit time integration
% - Convected Particle Domain Interpolation formulation
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
%   in line with the CPDI MPM as described in the published literature.
%
% Reference:
% -----------
% Sadeghirad, A., Brannon, R. M., & Burghardt, J. (2011). A convected particle domain 
% interpolation technique to extend applicability of the material point method for 
% problems involving massive deformations. International Journal for numerical methods
% in Engineering, 86(12), 1435-1456.
%
% Key Variables:
% x, y        - Nodal coordinates in x1 and x2 directions
% nod         - Nodal connectivity matrix
% xp          - Current positions of material points
% partmass    - Mass of material particles
% volpart     - Volume of particles
% defgrad     - Deformation gradient at each material point (2x2xMP)
% stress      - Cauchy stress at material points
% vpart       - Velocity of material points
% velgrad     - Velocity gradient at material points
% nodmass     - Nodal mass
% vnode       - Nodal velocity
% forcenodee  - External nodal forces
% forcenodei  - Internal nodal forces
% nodemo      - Nodal momentum
% sf          - Shape functions
% dsf         - Shape function derivatives
% -----------
% Version: 1.0
% Date: August 2025
%==========================================================================
%==========================================================================

clear all; close all; clc;
format long

% SIMULATION PARAMETERS
% Material properties
E = 2e8;               % Young's modulus (Pa)
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

% Mesh parameters
LLx = 10;              % Domain length in x-direction
LLy = 10;               % Domain length in y-direction

% Particle initialization parameters
width = 5;             % Width of particle block
height = 1;            % Height of particle block
particles_per_cell = 10; % Particles per cell in each direction

% Time parameters
t=0;
dt = 0.001;           % Time step size
T = 1;                 % Total simulation time
Tstep = round(T/dt);   % Number of time steps

% EULERIAN MESH GENERATION
% Generate structured grid of nodes
dxn = 1;
dyn = 1;
le = [dxn, dyn];

X = 0:dxn:LLx;
Y = 0:dyn:LLy;
[x1, y1] = meshgrid(Y, X);

% Reshape into column vectors
y = reshape(x1, length(Y)*length(X), 1);
x = reshape(y1, length(X)*length(Y), 1);
LOC = [x, y];          % Nodal coordinates

nnodesx = length(X);   % Number of nodes in x-direction
nnodesy = length(Y);   % Number of nodes in y-direction
nnodes = length(X)*length(Y); % Total number of nodes
nem = (LLx/dxn)*(LLy/dyn); % Number of Eulerian elements

% MATERIAL POINTS INITIALIZATION
% Create a grid of material points within the specified region
T_num_particles_x = round(width*particles_per_cell/dxn);
T_num_particles_y = round(height*particles_per_cell/dxn);
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
npart = length(xp(:, :)); % Number of particles


% Initialize particle vectors (r1, r2) defining particle domains
r1 = zeros(npart, 2);
r2 = zeros(npart, 2);
for i = 1:npart
    r1(i, :) = [(small_width/particles_per_cell)/2, 0];
    r2(i, :) = [0, (small_width/particles_per_cell)/2];
end

% Store initial particle vectors
r10 = r1;
r20 = r2;

% Calculate initial particle volumes
volpart = zeros(1, npart);
for i = 1:npart
    volpart(i) = 4 * abs(r1(i, 1)*r2(i, 2) - r1(i, 2)*r2(i, 1));
end

% Calculate particle masses
partmass = density * volpart;
volpart0 = volpart; % Store initial volumes

% Initialize particle corners
cr1 = zeros(npart, 2);
cr2 = zeros(npart, 2);
cr3 = zeros(npart, 2);
cr4 = zeros(npart, 2);
for i = 1:npart
    cr1(i, :) = xp(i, :) - r1(i, :) - r2(i, :);
    cr2(i, :) = xp(i, :) + r1(i, :) - r2(i, :);
    cr3(i, :) = xp(i, :) + r1(i, :) + r2(i, :);
    cr4(i, :) = xp(i, :) - r1(i, :) + r2(i, :);
end

% INITIALIZE FIELD VARIABLES
% Stress and strain tensors
stress = zeros(2, 2, npart); % Cauchy stress tensor
strain = zeros(2, 2, npart); % Strain tensor
nodemo=zeros(nnodes,2);
nodmass=zeros(nnodes,1);
% vnode=zeros(nnodes,2);
velgrad=zeros(2,2,npart);
forcenodee=zeros(nnodes,2);
forcenodei=zeros(nnodes,2);
forcenodetotal=zeros(nnodes,2);

% Deformation and velocity gradients
defgrad = repmat(eye(2), 1, 1, npart); % Initial deformation gradient
velgrad = zeros(2, 2, npart); % Velocity gradient

% Particle velocities
vpart = zeros(npart, 2);

% Shape functions and derivatives
shape = zeros(npart, nnodes);
dsf = zeros(2, nnodes, npart);

% Boundary conditions
ebc = 1; % Essential boundary condition flag
if ebc == 1
    ebcn = 1:nnodesx:nnodes; % Fixed nodes
else
    ebcn = [];
end

% Energy tracking
Ken = zeros(Tstep, 1); % Kinetic energy
Sen = zeros(Tstep, 1); % Strain energy
Pen = zeros(Tstep, 1); % Potential energy

% Data storage for analysis
strs = zeros(2, 2, Tstep, npart); % Stress history
strn = zeros(2, 2, Tstep, npart); % Strain history
xxp = zeros(Tstep, npart, 2); % Position history
vp = zeros(Tstep, npart, 2); % Velocity history

% Nodal quantities
nomo = zeros(Tstep, nnodes, 2); % Nodal momentum history
noms = zeros(Tstep, nnodes, 1); % Nodal mass history
vn = zeros(Tstep, nnodes, 2); % Nodal velocity history

% Force tracking
fint = zeros(Tstep, nnodes, 2); % Internal forces
fext = zeros(Tstep, nnodes, 2); % External forces
ftot = zeros(Tstep, nnodes, 2); % Total forces
trac = zeros(Tstep, nnodes, 2); % Traction forces

% Particle traction
ptraction = zeros(npart, 2);
ptr = zeros(npart, 2, Tstep); % Traction history

% Corner tracking
crr1 = zeros(npart, 2, Tstep);
crr2 = zeros(npart, 2, Tstep);
crr3 = zeros(npart, 2, Tstep);
crr4 = zeros(npart, 2, Tstep);

% PRE-COMPUTE SHAPE FUNCTIONS
% Initialize shape function arrays
sf1 = zeros(nnodes, npart);
sf2 = zeros(nnodes, npart);
sf3 = zeros(nnodes, npart);
sf4 = zeros(nnodes, npart);
sf = zeros(nnodes, npart);

% Pre-compute weight vectors for shape function derivatives
w1 = zeros(npart, 2);
w2 = zeros(npart, 2);
for i = 1:npart
    w1(i, :) = [r1(i, 2)-r2(i, 2), r2(i, 1)-r1(i, 1)];
    w2(i, :) = [r1(i, 2)+r2(i, 2), -r2(i, 1)-r1(i, 1)];
end

% MAIN SIMULATION LOOP
fprintf('Starting simulation with %d time steps...\n', Tstep);
tic

for nstep = 1:Tstep
    fprintf('Time step: %d/%d\n', nstep, Tstep);
    
    % Reset Eulerian mesh quantities
    nodemo(:) = 0;   % Nodal momentum
    nodmass(:) = 0;  % Nodal mass
    vnode(:) = 0;    % Nodal velocity
    velgrad(:) = 0; % Velocity gradient
    forcenodee(:) = 0; % External nodal forces
    forcenodei(:) = 0; % Internal nodal forces
    forcenodetotal(:) = 0; % Total nodal forces
    
    % Update particle corners based on current positions and vectors
    for i = 1:npart
        cr1(i, :) = xp(i, :) - r1(i, :) - r2(i, :);
        cr2(i, :) = xp(i, :) + r1(i, :) - r2(i, :);
        cr3(i, :) = xp(i, :) + r1(i, :) + r2(i, :);
        cr4(i, :) = xp(i, :) - r1(i, :) + r2(i, :);
    end
    
    % Store corner positions for analysis
    crr1(:, :, nstep) = cr1;
    crr2(:, :, nstep) = cr2;
    crr3(:, :, nstep) = cr3;
    crr4(:, :, nstep) = cr4;
    
    % Compute shape functions for each particle corner
    for j=1:npart
        for i=1:nnodes
            sf1(i,j)=max(0, (dxn - abs(cr1(j,1)-LOC(i,1)))/dxn).*max(0, (dyn - abs(cr1(j,2)-LOC(i,2)))/dyn);
            sf2(i,j)=max(0, (dxn - abs(cr2(j,1)-LOC(i,1)))/dxn).*max(0, (dyn - abs(cr2(j,2)-LOC(i,2)))/dyn);
            sf3(i,j)=max(0, (dxn - abs(cr3(j,1)-LOC(i,1)))/dxn).*max(0, (dyn - abs(cr3(j,2)-LOC(i,2)))/dyn);
            sf4(i,j)=max(0, (dxn - abs(cr4(j,1)-LOC(i,1)))/dxn).*max(0, (dyn - abs(cr4(j,2)-LOC(i,2)))/dyn);
        end
    end
    
    % Average shape functions from four corners
    for j = 1:npart
        for i = 1:nnodes
            sf(i, j) = 0.25 * (sf1(i, j) + sf2(i, j) + sf3(i, j) + sf4(i, j));
        end
    end
    
    % Compute shape function derivatives
    for j=1:npart
        for i=1:nnodes
            dsf(1,i,j)=1/volpart(j)*((sf1(i,j)-sf3(i,j))*w1(j,1)+(sf2(i,j)-sf4(i,j))*w2(j,1));
            dsf(2,i,j)=1/volpart(j)*((sf1(i,j)-sf3(i,j))*w1(j,2)+(sf2(i,j)-sf4(i,j))*w2(j,2));
        end
    end
    % Extract components for efficient computation
    strs11(:,:)=stress(1,1,:);
    strs22(:,:)=stress(2,2,:);
    strs12(:,:)=stress(1,2,:);
    dsfxx(:,:)=dsf(1,:,:);
    dsfyy(:,:)=dsf(2,:,:);
        
    % External force rate (zero in this simulation)
    frate = [0, -10];
    
    % MAPPING: PARTICLE TO GRID
    % Map mass to nodes
    nodmass = partmass * sf';
    
    % Map momentum to nodes
    nodemo = (partmass' .* vpart)' * sf';
    
    % External forces
    forcenodee = (frate' * nodmass)';
    
    % Internal forces (stress divergence)
    forcenodeix = ((-volpart' .* strs11)' * dsfxx' + (-volpart' .* strs12)' * dsfyy');
    forcenodeiy = ((-volpart' .* strs12)' * dsfxx' + (-volpart' .* strs22)' * dsfyy');
    forcenodei = [forcenodeix; forcenodeiy]';
    
    % Traction forces
    traction = ((volpart' .* ptraction)' * sf')';
    
    % Total nodal forces
    forcenodetotal = forcenodee + forcenodei + traction;
    
    % Update nodal momentum (impulse-momentum theorem)
    nodemo = nodemo + dt * forcenodetotal';
    
    % GRID UPDATE
    % Compute nodal acceleration and velocity
    anode=(forcenodetotal'./nodmass)';
    vnode=(nodemo./nodmass)';
    
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
    
    % MAPPING: GRID TO PARTICLE
    % Update particle velocities and positions
    vpart = vpart + sf' * dt * anode;
    xp = xp + sf' * dt * vnode;
    
    % Recompute nodal mass and momentum after update
    nodmass = partmass * sf';
    nodemo = (partmass' .* vpart)' * sf';
    
    % Update nodal velocities
    vnode=(nodemo./nodmass)';
    
    % Apply essential boundary conditions again
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

    
    % Compute deformation rate and update deformation gradient
    defrate = (velgrad + permute(velgrad, [2, 1, 3])) / 2 * dt;
    defgrad = (repmat(eye(2), 1, 1, npart) + velgrad * dt) .* defgrad;
    
    % Update particle volumes using Jacobian determinant
    J = arrayfun(@(k) det(defgrad(:, :, k)), 1:npart);
    volpart = volpart0 .* J;
    
    % Convert to Voigt notation for strain
    strain = [defrate(1, 1, :); defrate(2, 2, :); 2*defrate(1, 2, :)];
    
    % Compute stress increment (using elasticity tensor)
    stress_inc = zeros(3, 1, npart);
    for i = 1:npart
        stress_inc(:, 1, i) = C * strain(:, 1, i);
    end
    
    % Update stress (in matrix notation)
    for i = 1:npart
        stress(1, 1, i) = stress(1, 1, i) + stress_inc(1, 1, i);
        stress(2, 2, i) = stress(2, 2, i) + stress_inc(2, 1, i);
        stress(1, 2, i) = stress(1, 2, i) + stress_inc(3, 1, i);
        stress(2, 1, i) = stress(2, 1, i) + stress_inc(3, 1, i);
    end
    
    % Update particle domain vectors
    for i = 1:npart
        r1(i, 1:2) = (defgrad(1:2, 1:2, i) * r10(i, 1:2)')';
        r2(i, 1:2) = (defgrad(1:2, 1:2, i) * r20(i, 1:2)')';
    end
    
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
    
    % Update time
    t = t + dt;
end
toc


