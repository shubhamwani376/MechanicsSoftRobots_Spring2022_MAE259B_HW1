clc;clear all; close all;
%% Physical parameters
% Number of vertices
N = 50;

% Time step size
dt = 0.01; % second

% Rod length
RodLength = 1; % meter

% Discrete length
deltaL = RodLength / (N-1);

% Radius of spheres
R=0.013;
r=0.011;

% Density
rho = 2700;

% Young's modulus
Y = 70e9; % Using Y instead of E to avoid ambiguity

% Gravity
g = 9.8; % m/s^2

% Total time
totalTime = 1; % seconds

% Additional Load of 2000N downwards, at 0.75m
d = 0.75; %distance of load metres
P = 2000; %load in newtons

% Utility quantities
ne = N - 1; % Number of edges
EI = Y * pi * (R^4-r^4) / 4;
EA = Y * pi * (R^2-r^2);

% Geometry
nodes = zeros(N, 2);
for c = 1:N
    nodes(c,1) = (c-1) * deltaL;
end

% Mass matrix
M = zeros(2*N,2*N);
for c= 1:N
    M(2*c,2*c)=(pi*(R^2-r^2)*RodLength*rho)/(N-1); %from book
    M(2*c-1,2*c-1)= (pi*(R^2-r^2)*RodLength*rho)/(N-1); %from book
end
%for plotting the DER vs Euler theory.
% pcounter = 1:20;
% plotdata = zeros(3,length(pcounter));
% 
% for i = pcounter
%     P=1000*i;

% Gravity
    W = zeros(2*N,1);
    for c= 1:N
        W(2*c)=-M(2*c,2*c)*g;
    end
    W(2*round(d/deltaL,0))=W(2*round(d/deltaL,0))-P;
    
    % Initial DOF vector
    q0 = zeros(2*N,1);
    for c=1:N
        q0(2*c-1) = nodes(c,1); % x coordinate
        q0(2*c) = nodes(c,2); % y coordinate
    end
    
    % New position and velocity
    q = q0; % DOF vector
    u = (q - q0) / dt; % Velocity vector
    
    % Number of time steps
    Nsteps = round( totalTime / dt );
    ymax = zeros(Nsteps, 1); % ymax-position
    % all_mid_v = zeros( Nsteps, 1); % y-velocity of R2
    % 
    ymax(1) = 0;
    % all_mid_v(1) = u(4);
    
    % Tolerance
    tol = (EI / RodLength^2) * 1e-3;
    
    % Free and fixed index
    free_index = [3:1:(2*N-1)];
    fixed_index = [1,2,2*N];
    
    % Time marching scheme
    for c=2:Nsteps
        fprintf('Time = %f\n', (c-1) * dt );
        q_fixed = [nodes(1,1);nodes(1,2);nodes(N,2)];
        q = q0; % Guess
        
        % Updating the fixed index
        q(fixed_index) = q_fixed;
        q_free = q(free_index);
        
        % Newton Raphson
        err = 10 * tol;
        while err > tol        
            % Inertia
            f = M / dt * ( (q-q0) / dt - u );
            J = M / dt^2;
    
            % Elastic forces
    
            % Linear spring 1 between nodes 1 and 2
            for k=1:N-1
                xk = q(2*k-1);
                yk = q(2*k);
                xkp1 = q(2*k+1);
                ykp1 = q(2*k+2);
                dF = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
                dJ = hessEs(xk, yk, xkp1, ykp1, deltaL, EA);
                f(2*k-1:2*k+2) = f(2*k-1:2*k+2) + dF;
                J(2*k-1:2*k+2,2*k-1:2*k+2) = J(2*k-1:2*k+2,2*k-1:2*k+2) + dJ;
            end
            
            % Bending spring between nodes 1, 2, and 3
            for k=2:N-1
                xkm1 = q(2*k-3);
                ykm1 = q(2*k-2);
                xk = q(2*k-1);
                yk = q(2*k);
                xkp1 = q(2*k+1);
                ykp1 = q(2*k+2);
                curvature0 = 0;
                dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, deltaL, EI);
                dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, deltaL, EI);
                f(2*k-3:2*k+2) = f(2*k-3:2*k+2) + dF;
                J(2*k-3:2*k+2,2*k-3:2*k+2) = J(2*k-3:2*k+2,2*k-3:2*k+2) + dJ;
            end
    
            % Weight
            f = f - W;
            
            % Update
            f_free = f(free_index);
            J_free = J(free_index, free_index);
            q_free = q_free - J_free \ f_free;
            err = sum( abs(f_free) );
            
            % Plug free DOFs back into the full DOF vector
            q(free_index) = q_free;
    
        end
        
        % Update
        u = (q - q0) / dt; % Velocity
        q0 = q; % Old position
        ymax(c)=min(q);
        figure(1);
        plot( q(1:2:end), q(2:2:end), 'ro-');
        axis equal
        drawnow
    end
    
    figure(2);
    timeArray = (1:Nsteps) * dt;
    plot(timeArray, ymax, 'k-');
    xlabel('Time, t [sec]');
    ylabel('max displacement in y, y [meter]');
    
    cebm = min(d,RodLength-d);
    ymaxEBM=(P*cebm*(RodLength^2-cebm^2)^1.5)/(9*sqrt(3)*EI*RodLength);
    fprintf('Predicted from Euler Theory = %f\n', ymaxEBM );

%     plotdata(:,i) = [P;-1*ymax;ymaxEBM];
% end
% figure(2);
% timeArray = (1:Nsteps) * dt;
% plot(plotdata(1,:),plotdata(2:3,:));
% legend('ymax','ymaxEBM');
% xlabel('Load, P [N]');
% ylabel('max displacement in y, y [meter]');








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions below this point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI)

%
% This function returns the derivative of bending energy E_k^b with respect
% to x_{k-1}, y_{k-1}, x_k, y_k, x_{k+1}, and y_{k+1}.
%
% curvature0 is the "discrete" natural curvature [dimensionless] at node (xk, yk).
% l_k is the voronoi length of node (xk, yk).
% EI is the bending stiffness.
%

node0 = [xkm1, ykm1, 0];
node1 = [xk, yk, 0];
node2 = [xkp1, ykp1, 0];
%     m1e, 
m2e = [0 0 1];
%     m1f,
m2f = [0 0 1];
kappaBar = curvature0;

%% Computation of gradient of the two curvatures
gradKappa = zeros(6,1);

ee = node1 - node0;
ef = node2 - node1;

norm_e = norm(ee);
norm_f = norm(ef);

te = ee / norm_e;
tf = ef / norm_f;

% Curvature binormal
kb = 2.0 * cross(te, tf) / (1.0 + dot(te, tf));

chi = 1.0 + dot(te, tf);
tilde_t = (te + tf) / chi;
tilde_d2 = (m2e + m2f) / chi;

% Curvatures
kappa1 = kb(3); % 0.5 * dot( kb, m2e + m2f); % CHECKED

Dkappa1De = 1.0 / norm_e * (-kappa1 * tilde_t + cross(tf,tilde_d2));
Dkappa1Df = 1.0 / norm_f * (-kappa1 * tilde_t - cross(te,tilde_d2));

gradKappa(1:2, 1) = -Dkappa1De(1:2);
gradKappa(3:4, 1) = Dkappa1De(1:2) - Dkappa1Df(1:2);
gradKappa(5:6, 1) = Dkappa1Df(1:2);

%% Gradient of Eb
dkappa = kappa1 - kappaBar;
dF = gradKappa * EI * dkappa / l_k;
end

function F = gradEs(xk, yk, xkp1, ykp1, l_k, EA)
%
% This function returns the derivative of stretching energy E_k^s with 
% respect to x_{k-1}, y_{k-1}, x_k, and y_k.
%
F = zeros(4,1);

F(1) = -(0.1e1 - sqrt((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k) * ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1) / l_k * (-0.2e1 * xkp1 + 0.2e1 * xk);
F(2) = -(0.1e1 - sqrt((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k) * ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1) / l_k * (-0.2e1 * ykp1 + 0.2e1 * yk);
F(3) = -(0.1e1 - sqrt((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k) * ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1) / l_k * (0.2e1 * xkp1 - 0.2e1 * xk);
F(4) = -(0.1e1 - sqrt((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k) * ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1) / l_k * (0.2e1 * ykp1 - 0.2e1 * yk);

F = 0.5 * EA * l_k * F;

end
function dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI)

%
% This function returns the derivative of bending energy E_k^b with respect
% to x_{k-1}, y_{k-1}, x_k, y_k, x_{k+1}, and y_{k+1}.
%
% curvature0 is the "discrete" natural curvature [dimensionless] at node (xk, yk).
% l_k is the voronoi length of node (xk, yk).
% EI is the bending stiffness.
%

node0 = [xkm1, ykm1, 0];
node1 = [xk, yk, 0];
node2 = [xkp1, ykp1, 0];
%     m1e, 
m2e = [0 0 1];
%     m1f,
m2f = [0 0 1];
kappaBar = curvature0;

%% Computation of gradient of the two curvatures
gradKappa = zeros(6,1);

ee = node1 - node0;
ef = node2 - node1;

norm_e = norm(ee);
norm_f = norm(ef);

te = ee / norm_e;
tf = ef / norm_f;

% Curvature binormal
kb = 2.0 * cross(te, tf) / (1.0 + dot(te, tf));

chi = 1.0 + dot(te, tf);
tilde_t = (te + tf) / chi;
tilde_d2 = (m2e + m2f) / chi;

% Curvatures
kappa1 = kb(3); % 0.5 * dot( kb, m2e + m2f); % CHECKED

Dkappa1De = 1.0 / norm_e * (-kappa1 * tilde_t + cross(tf,tilde_d2));
Dkappa1Df = 1.0 / norm_f * (-kappa1 * tilde_t - cross(te,tilde_d2));

gradKappa(1:2, 1) = -Dkappa1De(1:2);
gradKappa(3:4, 1) = Dkappa1De(1:2) - Dkappa1Df(1:2);
gradKappa(5:6, 1) = Dkappa1Df(1:2);

%% Computation of hessian of the two curvatures
DDkappa1 = zeros(6, 6);
% DDkappa2 = zeros(11, 11);

norm2_e = norm_e^2;
norm2_f = norm_f^2;

tt_o_tt = tilde_t' * tilde_t; % must be 3x3. tilde_t is 1x3
tmp = cross(tf, tilde_d2);
tf_c_d2t_o_tt = tmp' * tilde_t; % must be 3x3
tt_o_tf_c_d2t = tf_c_d2t_o_tt'; % must be 3x3
kb_o_d2e = kb' * m2e; % must be 3x3
d2e_o_kb = kb_o_d2e'; % must be 3x3

Id3 = eye(3);
D2kappa1De2 ...
    = 1.0 / norm2_e * (2 * kappa1 * tt_o_tt - tf_c_d2t_o_tt - tt_o_tf_c_d2t) ...
    - kappa1 / (chi * norm2_e) * (Id3 - te'*te) ...
    + 1.0 / (4.0 * norm2_e) * (kb_o_d2e + d2e_o_kb);

tmp = cross(te, tilde_d2);
te_c_d2t_o_tt = tmp' * tilde_t;
tt_o_te_c_d2t = te_c_d2t_o_tt';
kb_o_d2f = kb' * m2f;
d2f_o_kb = kb_o_d2f';

D2kappa1Df2 ...
    = 1.0 / norm2_f * (2 * kappa1 * tt_o_tt + te_c_d2t_o_tt + tt_o_te_c_d2t) ...
    - kappa1 / (chi * norm2_f) * (Id3 - tf'*tf) ...
    + 1.0 / (4.0 * norm2_f) * (kb_o_d2f + d2f_o_kb);

D2kappa1DeDf ...
    = -kappa1/(chi * norm_e * norm_f) * (Id3 + te'*tf) ...
    + 1.0 / (norm_e*norm_f) * (2 * kappa1 * tt_o_tt - tf_c_d2t_o_tt + ...
    tt_o_te_c_d2t - crossMat(tilde_d2));
D2kappa1DfDe = D2kappa1DeDf';

% Curvature terms
DDkappa1(1:2, 1:2)  =   D2kappa1De2(1:2, 1:2);
DDkappa1(1:2, 3:4)  = - D2kappa1De2(1:2, 1:2) + D2kappa1DeDf(1:2, 1:2);
DDkappa1(1:2, 5:6) =               - D2kappa1DeDf(1:2, 1:2);
DDkappa1(3:4, 1:2)  = - D2kappa1De2(1:2, 1:2)                + D2kappa1DfDe(1:2, 1:2);
DDkappa1(3:4, 3:4)  =   D2kappa1De2(1:2, 1:2) - D2kappa1DeDf(1:2, 1:2) - ...
    D2kappa1DfDe(1:2, 1:2) + D2kappa1Df2(1:2, 1:2);
DDkappa1(3:4, 5:6) =                 D2kappa1DeDf(1:2, 1:2)                - D2kappa1Df2(1:2, 1:2);
DDkappa1(5:6, 1:2)  =                              - D2kappa1DfDe(1:2, 1:2);
DDkappa1(5:6, 3:4)  =                                D2kappa1DfDe(1:2, 1:2) - D2kappa1Df2(1:2, 1:2);
DDkappa1(5:6, 5:6) =                                               D2kappa1Df2(1:2, 1:2);

%% Hessian of Eb
dkappa = kappa1 - kappaBar;
dJ = 1.0 / l_k * EI * gradKappa * transpose(gradKappa);
temp = 1.0 / l_k * dkappa * EI;
dJ = dJ + temp * DDkappa1;

end

function A = crossMat(a)
A = [0, -a(3), a(2); ...
    a(3), 0, -a(1); ...
    -a(2), a(1), 0];
end
function J = hessEs(xk, yk, xkp1, ykp1, l_k, EA)
%
% This function returns the 4x4 hessian of the stretching energy E_k^s with
% respect to x_k, y_k, x_{k+1}, and y_{k+1}.
%
J11 = (1 / ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k ^ 2 * (-2 * xkp1 + 2 * xk) ^ 2) / 0.2e1 + (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.3e1 / 0.2e1)) / l_k * ((-2 * xkp1 + 2 * xk) ^ 2) / 0.2e1 - 0.2e1 * (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1)) / l_k;
J12 = (1 / ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k ^ 2 * (-2 * ykp1 + 2 * yk) * (-2 * xkp1 + 2 * xk)) / 0.2e1 + (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.3e1 / 0.2e1)) / l_k * (-2 * xkp1 + 2 * xk) * (-2 * ykp1 + 2 * yk) / 0.2e1;
J13 = (1 / ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k ^ 2 * (2 * xkp1 - 2 * xk) * (-2 * xkp1 + 2 * xk)) / 0.2e1 + (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.3e1 / 0.2e1)) / l_k * (-2 * xkp1 + 2 * xk) * (2 * xkp1 - 2 * xk) / 0.2e1 + 0.2e1 * (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1)) / l_k;
J14 = (1 / ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k ^ 2 * (2 * ykp1 - 2 * yk) * (-2 * xkp1 + 2 * xk)) / 0.2e1 + (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.3e1 / 0.2e1)) / l_k * (-2 * xkp1 + 2 * xk) * (2 * ykp1 - 2 * yk) / 0.2e1;
J22 = (1 / ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k ^ 2 * (-2 * ykp1 + 2 * yk) ^ 2) / 0.2e1 + (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.3e1 / 0.2e1)) / l_k * ((-2 * ykp1 + 2 * yk) ^ 2) / 0.2e1 - 0.2e1 * (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1)) / l_k;
J23 = (1 / ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k ^ 2 * (2 * xkp1 - 2 * xk) * (-2 * ykp1 + 2 * yk)) / 0.2e1 + (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.3e1 / 0.2e1)) / l_k * (-2 * ykp1 + 2 * yk) * (2 * xkp1 - 2 * xk) / 0.2e1;
J24 = (1 / ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k ^ 2 * (2 * ykp1 - 2 * yk) * (-2 * ykp1 + 2 * yk)) / 0.2e1 + (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.3e1 / 0.2e1)) / l_k * (-2 * ykp1 + 2 * yk) * (2 * ykp1 - 2 * yk) / 0.2e1 + 0.2e1 * (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1)) / l_k;
J33 = (1 / ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k ^ 2 * (2 * xkp1 - 2 * xk) ^ 2) / 0.2e1 + (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.3e1 / 0.2e1)) / l_k * ((2 * xkp1 - 2 * xk) ^ 2) / 0.2e1 - 0.2e1 * (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1)) / l_k;
J34 = (1 / ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k ^ 2 * (2 * ykp1 - 2 * yk) * (2 * xkp1 - 2 * xk)) / 0.2e1 + (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.3e1 / 0.2e1)) / l_k * (2 * xkp1 - 2 * xk) * (2 * ykp1 - 2 * yk) / 0.2e1;
J44 = (1 / ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k ^ 2 * (2 * ykp1 - 2 * yk) ^ 2) / 0.2e1 + (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.3e1 / 0.2e1)) / l_k * ((2 * ykp1 - 2 * yk) ^ 2) / 0.2e1 - 0.2e1 * (0.1e1 - sqrt(((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2)) / l_k) * (((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1)) / l_k;

J = [J11 J12 J13 J14;
     J12 J22 J23 J24;
     J13 J23 J33 J34;
     J14 J24 J34 J44];

J = 0.5 * EA * l_k * J;

end
