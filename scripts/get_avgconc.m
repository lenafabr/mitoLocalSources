function [Mtot, Msqrtot] = get_avgconc(param, s_dist)
% param = structure of kinetic parameters
% s_dist = distribution of cluster sizes
% Returns: steady-state predictions
% Mtot = average concentration
% Msqrtot = average concentration squared

% dimensionless parameters
beta = param.kdecay/param.khop;
gamma = param.kdecay/param.kpump;


% Define the Hypergeometric Model Function
% This function takes a guess for gamma, and outputs the predicted curve
hyper_model_precalc = @(gamma, beta, C0, s_array, ns, N, precalc) arrayfun(@(s) ...
    C0 * sum(hygepdf(max(0, s-(N-ns)):min(s, ns), N, ns, s) .* ...
    precalc), s_array);


Mtot = 0;
Msqrtot = 0;
for s = find(s_dist > 0)
    % solve for concentration (mean and mean squared) in clusters
    % of size s, with j sources in that cluster
    [Mj, Mjsqr] = solve_cluster_M(s, gamma, beta, param.N, param.nsource);
    % average over the number of sources j in a cluster of size s
    % and also average over cluster size distribution
    Mtot = Mtot + s/param.N*s_dist(s)*hyper_model_precalc(gamma, beta, 1, ...
        s, param.nsource, param.N, Mj);
    Msqrtot = Msqrtot + s/param.N*s_dist(s)*hyper_model_precalc(gamma, beta, 1, ...
        s, param.nsource, param.N, Mjsqr);
end