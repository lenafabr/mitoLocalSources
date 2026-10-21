function average_var_concentration
% FIX THIS
Mtot = 0;
        Msqrtot = 0;
        for s = find(s_dist > 0)
            % solve for concentration (mean and mean squared) in clusters
            % of size s, with j sources in that cluster            
            [Mj, Mjsqr] = solve_cluster_M(s, gamma, beta, N, nsources(j));
            % average over the number of sources j in a cluster of size s
            % and also average over cluster size distribution
            Mtot = Mtot + s/N*s_dist(s)*hyper_model_precalc(gamma, beta, 1, ...
                s, nsources(j), N, Mj);
            Msqrtot = Msqrtot + s/N*s_dist(s)*hyper_model_precalc(gamma, beta, 1, ...
                s, nsources(j), N, Mjsqr);
        end

end