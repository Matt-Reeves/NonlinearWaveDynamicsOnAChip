%% Simulation.m
% Hydrodynamic simulation of fluid in a wavetank, original code by Matt
% Reeves. 

function Simulation(savename,savedir,amp,N,d,periods,tstep,Gamma)
%% Input Variables
% savename = 'Testing_1';   % name for data save file 
% amp     = 0.1;          % amplitude 
% N       = 256;           % # of grid points 
% d      = 0.1;          % fluid depth, 0.001 = 1000:1
% periods = 40;           % fundumental periods to simulate 
kappa   = 0.0;          % surface tension strength
% Gamma   = 0.002;          % damping coefficicent
visc    = 0.00;         % viscocity coefficient
waveType = 'VdW';       % must be "VdW" or "gravity"

%% Simulation Initialization
omega0 = sqrt(tanh(d));
amplitude = amp*d;
Ac = amplitude*omega0^2;
simparams = simparams_initialize(N,d,waveType,Ac,kappa,Gamma,visc,periods,tstep);
simparams.dt;
% simparams.dt = 1/2^5                                       %fine timestep
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
    xlim([-pi pi])%, xticks([-pi -pi/2 0]), xticklabels({'0','\pi/2','\pi'})
    hold off
    drawnow
end
toc
disp('Simulation finished')
%% Simulation Output Save
savefile = strcat(savedir,savename);
save(savefile)

disp(['Data saved as ',savename,'.mat in ',savefile])
  