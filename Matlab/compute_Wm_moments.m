% =================================================================================
% Momentos of the random variable Wm
% Based on:
% Gamboa, M., López-García, M., & Lopez-Herrero, M. J. (2021)
% "A stochastic SVIR model with imperfect vaccine and external source
% of infection"European Workshop on Performance Engineering, pp. 197–209.
% https://doi.org/10.1007/978-3-030-91825-5_12
%

% ============================================================
% Example of use of momentoW_M
% ============================================================

clear; clc;

%% Model parameters (example values)
N     = 100;     % Total population
m     = 10;     % Initial vaccinated threshold
v0    = 50;      % Initial vaccinated individuals
s0    = 49;     % Initial susceptible individuals
i0    = 1;      % Initial infected individuals

beta  = 12;    % Infection rate
gamma = 1;    % Recovery rate
xi    = 0.01 ;   % External source of infection
h     = 0.2;    % Vaccine failure probability


%% Compute first and second moments of W
[K1, K2] = momentoW_M(v0, s0, i0, m, N, gamma, xi, beta, h);


%% Display results
fprintf('First moment E[W]  = %.6f\n', K1);
fprintf('Second moment E[W^2] = %.6f\n', K2);

%% variance of W
VarW = K2 - K1^2;
fprintf('Variance Var(W) = %.6f\n', VarW);



function [K1, K2] = momentoW_M(v0, s0, i0, m, N, gamma, xi, beta, h)
% INPUT:
%   v0, s0, i0 : initial conditions
%   N          : population parameters
%    m         : number of individuals get the infection
%   gamma     : recovery rate
%   xi        : external infection rate
%   beta      : transmission rate
%   h         : vaccine failure probability
%
% OUTPUT:
%   K1 : first moment of W
%   K2 : second moment of W
% =================================================================================
%% Initialization
M_Gmenos1_k1 = zeros(s0+1, m+1);   % First moment (k = 1)
M_Gmenos1_k2 = zeros(s0+1, m+1);   % Second moment (k = 2)

%% Main recursion over g
for g = (N-m+1):(v0+s0)

    % k = 0 (auxiliary matrix)
    MGK0 = ones(s0+1, N-g+1); %#ok<NASGU>

    %% =======================
    % k = 1 (First moment)
    %% =======================
    MGK1 = zeros(s0+1, N-g+1);

    for s = max(0, g-v0):s0
        for i = max(0, N-m-g+1):(N-g)

            if (g == (N-m+1)) && (i == 0)
                MGK1(s+1,i+1) = 1 / q(beta, xi, gamma, h, N, s, i, g-s);
            end

            if (g == (N-m+1)) && (i ~= 0)
                MGK1(s+1,i+1) = ...
                    (1 + gamma*i*MGK1(s+1,i)) / ...
                     q(beta, xi, gamma, h, N, s, i, g-s);
            end

            if (g ~= (N-m+1)) && (i == 0)
                if s == 0
                    MGK1(s+1,i+1) = ...
                        (1 + eta(beta, xi, h, N, g-s, i) * ...
                         M_Gmenos1_k1(s+1,i+2)) / ...
                         q(beta, xi, gamma, h, N, s, i, g-s);
                else
                    MGK1(s+1,i+1) = ...
                        (1 + lambda(beta, xi, N, s, i) * ...
                         M_Gmenos1_k1(s,i+2) + ...
                         eta(beta, xi, h, N, g-s, i) * ...
                         M_Gmenos1_k1(s+1,i+2)) / ...
                         q(beta, xi, gamma, h, N, s, i, g-s);
                end
            end

            if (g ~= (N-m+1)) && (i ~= 0)
                if s == 0
                    MGK1(s+1,i+1) = ...
                        (1 + eta(beta, xi, h, N, g-s, i) * ...
                         M_Gmenos1_k1(s+1,i+2) + ...
                         gamma*i*MGK1(s+1,i)) / ...
                         q(beta, xi, gamma, h, N, s, i, g-s);
                else
                    MGK1(s+1,i+1) = ...
                        (1 + lambda(beta, xi, N, s, i) * ...
                         M_Gmenos1_k1(s,i+2) + ...
                         eta(beta, xi, h, N, g-s, i) * ...
                         M_Gmenos1_k1(s+1,i+2) + ...
                         gamma*i*MGK1(s+1,i)) / ...
                         q(beta, xi, gamma, h, N, s, i, g-s);
                end
            end
        end
    end

    %% =======================
    % k = 2 (Second moment)
    %% =======================
    MGK2 = zeros(s0+1, N-g+1);

    for s = max(0, g-v0):s0
        for i = max(0, N-m-g+1):(N-g)

            if (g == (N-m+1)) && (i == 0)
                MGK2(s+1,i+1) = ...
                    2 * MGK1(s+1,i+1) / ...
                    q(beta, xi, gamma, h, N, s, i, g-s);
            end

            if (g == (N-m+1)) && (i ~= 0)
                MGK2(s+1,i+1) = ...
                    (2 * MGK1(s+1,i+1) + ...
                     gamma*i*MGK2(s+1,i)) / ...
                     q(beta, xi, gamma, h, N, s, i, g-s);
            end

            if (g ~= (N-m+1)) && (i == 0)
                if s == 0
                    MGK2(s+1,i+1) = ...
                        (2 * MGK1(s+1,i+1) + ...
                         eta(beta, xi, h, N, g-s, i) * ...
                         M_Gmenos1_k2(s+1,i+2)) / ...
                         q(beta, xi, gamma, h, N, s, i, g-s);
                else
                    MGK2(s+1,i+1) = ...
                        (2 * MGK1(s+1,i+1) + ...
                         lambda(beta, xi, N, s, i) * ...
                         M_Gmenos1_k2(s,i+2) + ...
                         eta(beta, xi, h, N, g-s, i) * ...
                         M_Gmenos1_k2(s+1,i+2)) / ...
                         q(beta, xi, gamma, h, N, s, i, g-s);
                end
            end

            if (g ~= (N-m+1)) && (i ~= 0)
                if s == 0
                    MGK2(s+1,i+1) = ...
                        (2 * MGK1(s+1,i+1) + ...
                         eta(beta, xi, h, N, g-s, i) * ...
                         M_Gmenos1_k2(s+1,i+2) + ...
                         gamma*i*MGK2(s+1,i)) / ...
                         q(beta, xi, gamma, h, N, s, i, g-s);
                else
                    MGK2(s+1,i+1) = ...
                        (2 * MGK1(s+1,i+1) + ...
                         lambda(beta, xi, N, s, i) * ...
                         M_Gmenos1_k2(s,i+2) + ...
                         eta(beta, xi, h, N, g-s, i) * ...
                         M_Gmenos1_k2(s+1,i+2) + ...
                         gamma*i*MGK2(s+1,i)) / ...
                         q(beta, xi, gamma, h, N, s, i, g-s);
                end
            end
        end
    end

    % Update previous-step matrices
    M_Gmenos1_k1 = MGK1;
    M_Gmenos1_k2 = MGK2;
end

%% Output moments
K1 = M_Gmenos1_k1(s0+1, i0+1);
K2 = M_Gmenos1_k2(s0+1, i0+1);

end

%% ============================================================
% Auxiliary rate functions
%% ============================================================

function x = lambda(beta, xi, N, s, i)
x = beta*s*i/N + xi*s;
end

function x = eta(beta, xi, h, N, v, i)
x = beta*v*i*h/N + xi*h*v;
end

function x = q(beta, xi, gamma, h, N, s, i, v)
x = beta*s*i/N + xi*s + beta*v*i*h/N + xi*h*v + gamma*i;
end
