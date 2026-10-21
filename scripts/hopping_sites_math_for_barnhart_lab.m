%% -----------------
% This code will calculate steady-state average concentration of material
% produced by local injection into a mitochondrial network
%
% The injection of material at each site follows saturating kinetics
% dC/dt = kpump*(C0-C)
% Where the saturation concentration is assumed C0 = 1 in arbitrary units
% all material concentrations are thus relative to C0
% 
% This example sweeps over different number of local injection site, for
% several different values of injection rate.
%% ------------------

% You will need to access the networktools package
% It can be pulled from:  https://github.com/lenafabr/networktools

% Edit this top point to where you networktools package is located
addpath('~/proj/networktools')

%% parameters saved in structure
param = struct();

param.khop = 1.25; % rate of hopping (units: per time) - how fast do sources move to a new mito?
param.kdecay = 1; % rate of decay (units: per time) - how fast does injected material turn over?

% injection rate at each local site (units: per time)
% here we try four different values
kmax = 4;
kpumplist = logspace(-1,2,kmax);

% number of sources, must be an array of integers
% Here we scan over several different source numbers
smax = 21;
nsourcelist = linspace(1,81,smax);

% total size of the network (number of mitochondrial units, eg: representing 0.5um-long chunks of mitochondria)
% Must be an integer. Don't change it for now unless you want to supply your own
% cluster size distributions.
param.N = 250; 

%% To set the connectivity of the network, we did a simulation sweep over two parameters:
% The two parameters C1, C2 describe the balance of fusion:fission for
% different types of events
% C1 ~ x_2/x_1^2 ~ ratio of tip-tip fusion over edge fission
% C2 ~ x_3/(x_1*x_2) ~ ratio of tip-side fusion over junction fission
% (where x_i is the number of nodes of degree i)
% In principle, these can be extracted from experimental data (would need
% to account for length density per volume as well). For now just treat
% them as knobs to control connectivity.
% The connectivity is used to define the distribution of cluster sizes,
% based on the steady state reached in simulations

% load the simulation data
load('barnhart_demo.mat')
% get the distribution of cluster sizes for each C1, C2 used in the
% simulations
meanclustersizearr = sum(mean(meanclustdistrib(:,:,:,:),4) .* reshape((1:param.N).^2, [1 1 param.N]), 3)/param.N^2;


% Decide which value of c1, c2 you want.
% Larger c2 = more percolated, branched networks.
% Larger c1 but small c2 = more snake-like networks with long tubes
% c1idx, c2idx  are indices into the arrays that store the simulated values of C1,
% C2
c1idx = 6;
c2idx = 5;

disp("C1, C2 selected values: "+c1fromnodes(c1idx,c2idx)+", "+c2fromnodes(c1idx,c2idx))

%% Visualize the parameter choices for network connectivities
% use this to play around to select the type of network you want
figure
tiledlayout(1,3)
nexttile;

% Plot (a) shows which connectivity parameters (C1,C2) you chose, 
% and the mean cluster size fraction they correspond to 
% (mean cluster size / total network size)
s = surf(c1fromnodes,c2fromnodes,meanclustersizearr);
s.EdgeAlpha = 0.2; 
hold on
scatter3(c1fromnodes(c1idx,c2idx), c2fromnodes(c1idx,c2idx), param.N, 200, "black", "filled", "pentagram")
view(2)
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'defaultTextInterpreter','latex','TickLabelInterpreter','latex','FontSize',16)
xlabel("$C_1$ (tip-tip fusion:fission balance)")
ylabel("$C_2$ (tip-side fusion:fission balance)")
cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 16;
clim([0 1]);
cb.Ticks = [0 0.5 1];
cb.Label.String = 'mean cluster size fraction';
    title('(a) Choice of network connectivity')

% must satisfy sum(s_dist.*(1:N)) = N. This is the cluster size
% distribution for your choice of (c1idx, c2idx)
s_dist = reshape(mean(meanclustdistrib(c1idx,c2idx,:,:),4), [1 param.N]);

% Plot (b) calculates the distribution of cluster sizes (if you pick a
% random mitochondrial unit, what is the probability it is in a cluster of
% a given size)
nexttile;
plot(1:param.N,s_dist.*(1:param.N)/param.N,'LineWidth',2)
set(gca,'defaultTextInterpreter','latex','TickLabelInterpreter','latex','FontSize',16)
xlabel("Cluster size, $s$")
ylabel(["for a random mito unit,"  "probability it is in cluster of this size"])
title('(b) Cluster size distribution')

xlim([0 param.N+1])

set(gcf,'color','w');
grid on

% Plot (c) shows an examp snapshot of the network
nexttile;
savenetworks(1,c1idx,c2idx).plotNetwork(); % plot example network
title("(c) Example network snapshot")

set(gcf, 'Position', [100, 400, 1800, 500]);


%% Calculate mean and variability of the steady-state concentration of material in the network

pred_conc = zeros(kmax,smax);
pred_var = zeros(kmax,smax);

for kc = 1:numel(kpumplist) % scan over different values of kpump 
    param.kpump = kpumplist(kc);
    for sc = 1:numel(nsourcelist) % scan over different total number of sources
        param.nsource = nsourcelist(sc);
        [Mtot, Msqrtot] = get_avgconc(param,s_dist);
        pred_conc(kc,sc) = Mtot;
        pred_var(kc,sc) = Msqrtot - Mtot^2;

        fprintf('kpump %g, nsource %g, avg conc %g, var conc %g\n', param.kpump, param.nsource,Mtot, pred_var(kc,sc))
    end
end

%% plot mean and variability
cm = parula(kmax);

figure

tiledlayout(1,2)
nexttile;
hold on
for kc = 1:kmax
    plot(nsourcelist,pred_conc(kc,:),'color',cm(kc,:),'linewidth',2,'DisplayName',"$k_p="+kpumplist(kc)+"$")
end
set(gca,'defaultTextInterpreter','latex','TickLabelInterpreter','latex','FontSize',16)
xlabel("$n_s$")
ylabel("$\langle M \rangle$")
legend("Interpreter","latex","FontSize",15)
title('(a) Average concentration')
xlim([0 max(nsourcelist)+1])
ylim([0 1])


nexttile;
hold on
for kc = 1:kmax
    plot(nsourcelist,sqrt(pred_var(kc,:))./pred_conc(kc,:),'color',cm(kc,:),'linewidth',2,'DisplayName',"$k_p="+kpumplist(kc)+"$")
end
set(gca,'defaultTextInterpreter','latex','TickLabelInterpreter','latex','FontSize',16)
xlabel("$n_s$")
ylabel("COV = $\sigma_M / \langle M \rangle$")
legend("Interpreter","latex","FontSize",15)
xlim([0 max(nsourcelist)+1])
ylim([0 4])
title('(b) Coefficient of variation (stdev/mean)')

set(gcf,'color','w');
grid on