%% run.m
% Hydrodynamic simulation of fluid in a wavetank, original code by Matt
% Reeves. 
clear all

%% Input Variables
savename = 'Test2D';   % name for data save file 
amp     = 0.1;          % amplitude 
N       = 124;           % # of grid points 
d      = 0.01;          % fluid depth, 0.001 = 1000:1
periods = 20;           % fundumental periods to simulate 
kappa   = 0.0;          % surface tension strength
Gamma   = 0.0;          % damping coefficicent
visc    = 0.00;         % viscocity coefficient
waveType = 'VdW';       % must be "VdW" or "gravity"

%% Simulation Initialization
d = d*3.14; % turns d into nice aspect ratio 
omega0 = sqrt(tanh(d));
amplitude = amp*d;
Ac = amplitude*omega0^2;
simparams = simparams_initialize(N,d,waveType,Ac,kappa,Gamma,visc,periods);
simparams.dt;
%simparams.dt = 1/2^5                                       %fine timestep
simparams.Nt        = (simparams.tfin/simparams.DT);       %coarse stepper
simparams.nt        = (simparams.DT/simparams.dt);         %fine stepper

simgrids = simgrids_initialize(simparams);

tempgrid = simgrids;
solution = struct; 
solution.simparams = simparams;
solution.Q0 = reshape(simgrids.Q,[N,3]);
%% Evolution 
disp(savename)
disp('Simulation parameters')
disp(struct2table(simparams))
disp('Simulation start')
tic
for jj = 1:simparams.Nt
    
    for kk = 1:simparams.nt
    tempgrid.Q = rk4(tempgrid,simparams);
    end
  
     %calculate conserved quantities
    quantities = conservedQuantities(tempgrid,simparams);
    energy.E(jj)   = quantities.E;
    energy.T(jj)   = quantities.T;
    energy.V(jj)   = quantities.V;
    volumeFlux(jj) = quantities.volumeFlux;
    eta_bar(jj)    = quantities.Ybar;
    
    %Save the solution for later processing
    solution.Q(1:N,1:3,jj) = reshape(tempgrid.Q,[N,3]);
    
    %Plot as we go
    figure(999)
    clf
    plot(tempgrid.Q(1:simparams.N),tempgrid.Q((1:simparams.N)+simparams.N),'.-')
    hold on
    ylim([-d 2*simparams.A]), yticks([-d 0]), yticklabels({'0','d'})
    xlim([-pi 0]), xticks([-pi -pi/2 0]), xticklabels({'0','\pi/2','\pi'})
    hold off
    drawnow
end
toc
disp('Simulation finished')
%% Simulation Plots
TimePlots(solution,simparams,tempgrid,energy)
%WaveAnimation(simparams,solution,simgrids,savename)
SurfacePlot(simparams,solution)
%% Simulation Output Save
save(savename)
disp(['Data saved as ',savename,'.mat'])
  