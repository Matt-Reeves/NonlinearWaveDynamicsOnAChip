function simparams = simparams_initialize(N,d,waveType,Ac,kappa,Gamma,visc,periods,tstep)
simparams = struct; 

%number of points
simparams.N = N;

%fluid parameters
simparams.waveType = waveType;                             % 'gravity' or 'VdW'
simparams.d        = d;                                    % fluid depth
simparams.AspectRatio = pi/d;                              % aspect ratio of wave tank
simparams.depth    = (2*10^-4)*(pi/d);                     % physical depth 
simparams.kappa    = kappa;                                % surface tension (ignoring this for now) 
simparams.Gamma    = Gamma;                                % Damping coefficient
simparams.visc     = visc;                                 % viscocity coefficient
simparams.periods     = periods;                           % periods to simulate for
simparams.tstep    = tstep;
% residual-related parameters.
if     strcmp(simparams.waveType,'VdW'),     simparams.Act   = +Ac;                               % target crest acceleteration
elseif strcmp(simparams.waveType,'gravity'), simparams.Act   = -Ac;                
end
simparams.alpha = -1e-4;                                % potential offset 

% time-related parameters
simparams.omega0    = sqrt((1 + kappa)*tanh(simparams.d)); % linear frequency
simparams.omega_max = sqrt((N/2 + kappa*(N/2)^3)*tanh(simparams.d*N/2)); % linear frequency
simparams.T         = 2*pi/simparams.omega0;               % linear period
simparams.tfin      = simparams.T*simparams.periods;                      % evolve 10T
simparams.DT        = simparams.T/simparams.tstep;                      % coarse timestep, 50 by default
simparams.dt        = 1/2/simparams.omega_max;           % fine timestep
simparams.Nt        = (simparams.tfin/simparams.DT);       % coarse stepper
simparams.nt        = (simparams.DT/simparams.dt);         % fine stepper
simparams.t         = 0;                                   % current time 

% Guess wave height for starting point, (linearized approximation)
simparams.A     =  simparams.Act/simparams.omega0.^2;

return




