%% solves local injection assuming mito network is static. Asks for a 
% cluster size s when there are ns total sites available and N total
% mito units, what is the mean concentration for each possible number of
% attached sites? So for s = 2 you can have 0, 1, or 2 sites.

function [M_bar, M2_bar] = solve_cluster_M(s, gamma, beta, N, ns)
    % SOLVE_CLUSTER_M Computes steady-state average concentrations \bar{M_j}
    % and second moments \bar{M_j^2} for a cluster of size s.
    %
    % Inputs:
    %   s     - Cluster size
    %   gamma - Ratio of decay rate to production rate (k_d / k_p)
    %   beta  - Ratio of decay rate to hop rate (k_d / k_hop)
    %   N     - Total number of particles
    %   ns    - Total number of sources (n_s)
    %
    % Outputs:
    %   M_bar  - Row vector of first moments \bar{M_j} for j in [j_min, j_max]
    %   M2_bar - Row vector of second moments \bar{M_j^2} for j in [j_min, j_max]

    % Target ratio p*
    p_star = ns / (N - ns);

    % Define state bounds for j
    j_min = max(0, s - (N - ns));
    j_max = min(s, ns);
    j_vals = j_min:j_max;
    n_states = length(j_vals);

    % Initialize linear system A * M_bar' = D
    A = zeros(n_states, n_states);
    D = zeros(n_states, 1);
    
    % Store intermediate rates and parameters for step 2
    F = zeros(1, n_states);
    G = zeros(1, n_states);
    M_inf_vec = zeros(1, n_states);

    % --- STEP 1: Construct and solve linear system for M_bar ---
    for k = 1:n_states
        j = j_vals(k);

        % Normalized jump rates (k_hop units)
        R_up   = (j < j_max) * (s - j) * p_star;
        R_down = (j > j_min) * j;
        k_esc  = R_up + R_down;
        
        % Normalized relaxation rate (k_hop units)
        k_rel = (j * beta) / (s * gamma) + beta;

        % Asymptotic steady-state concentration
        if j == 0
            M_inf = 0;
        else
            M_inf = j / (j + s * gamma);
        end
        M_inf_vec(k) = M_inf;

        % Retention factors
        if k_esc == 0
            F_j = 0;
            G_j = 0;
            alpha_j = 0;
        else
            F_j = k_esc / (k_esc + k_rel);
            G_j = k_esc / (k_esc + 2 * k_rel);
            alpha_j = R_down / k_esc;
        end
        
        F(k) = F_j;
        G(k) = G_j;

        % Tridiagonal matrix elements
        A(k, k) = 1;
        D(k) = M_inf * (1 - F_j);

        if k > 1
            A(k, k-1) = -F_j * alpha_j;
        end
        if k < n_states
            A(k, k+1) = -F_j * (1 - alpha_j);
        end
    end

    % Solve for first moment \bar{M_j}
    M_bar = (A \ D)';

    % --- STEP 2: Compute second moment \bar{M_j^2} for each state ---
    M2_bar = zeros(1, n_states);
    for k = 1:n_states
        M_inf = M_inf_vec(k);
        F_j = F(k);
        G_j = G(k);

        if F_j > 0
            delta_M = (M_bar(k) - M_inf) / F_j;
            M2_bar(k) = 2 * M_inf * M_bar(k) - M_inf^2 + (delta_M^2) * G_j;
        else
            % Static state with no escape transitions
            M2_bar(k) = M_inf^2;
        end
    end
end


% function M_bar = solve_cluster_M(s, gamma, beta, N, ns)
%     % SOLVE_CLUSTER_M Computes steady-state concentrations \bar{M_j} for a cluster of size s.
%     %
%     % Inputs:
%     %   s     - Cluster size
%     %   gamma - Ratio of decay rate to production rate (k_d / k_p)
%     %   beta  - Ratio of decay rate to hop rate (k_d / k_hop)
%     %   N     - Total number of particles
%     %   ns    - Total number of sources (n_s)
%     %
%     % Output:
%     %   M_bar - Row vector of steady-state concentrations for j in [j_min, j_max]
% 
%     % Calculate target ratio p* from total particle count and sources
%     p_star = ns / (N - ns);
% 
%     % Define state bounds for j
%     j_min = max(0, s - (N - ns));
%     j_max = min(s, ns);
%     j_vals = j_min:j_max;
%     n_states = length(j_vals);
% 
%     % Initialize linear system A * M_bar' = D
%     A = zeros(n_states, n_states);
%     D = zeros(n_states, 1);
% 
%     for k = 1:n_states
%         j = j_vals(k);
% 
%         % Dimensionless jump rates (normalized by k_hop)
%         if j < j_max
%             R_up = (s - j) * p_star;
%         else
%             R_up = 0;
%         end
% 
%         if j > j_min
%             R_down = j;
%         else
%             R_down = 0;
%         end
% 
%         k_esc_norm = R_up + R_down;
% 
%         % Asymptotic steady-state concentration M_inf for state j
%         if j == 0
%             M_inf = 0;
%         else
%             M_inf = j / (j + s * gamma);
%         end
% 
%         % Retention factor F_j and entry transition fraction alpha_j
%         if k_esc_norm == 0
%             F_j = 0;
%             alpha_j = 0;
%         else
%             F_j = k_esc_norm / (k_esc_norm + (j * beta) / (s * gamma) + beta);
%             alpha_j = R_down / k_esc_norm;
%         end
% 
%         % Construct tridiagonal system coefficients
%         A(k, k) = 1;
%         D(k) = M_inf * (1 - F_j);
% 
%         if k > 1
%             A(k, k-1) = -F_j * alpha_j;
%         end
%         if k < n_states
%             A(k, k+1) = -F_j * (1 - alpha_j);
%         end
%     end
% 
%     % Solve linear system and transpose to return a row vector
%     M_bar = (A \ D)';
% end