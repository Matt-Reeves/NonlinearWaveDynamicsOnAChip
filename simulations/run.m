%% SimSweep
% Rund the hydrodynamic simulation for given parameters and saves each
% output. 
clear 
%% Input Variables
%savename = 'QuickTest';   % name for data save file 
%savedir = 'C:\Users\s4606711\OneDrive - The University of Queensland\PhD\Hydrodynamics\Simulations\Data\Test\';
% amp     = 0.1;          % amplitude 
% N       = 2048;           % # of grid points 
% d      = 0.5;          % fluid depth, 0.001 = 1000:1
% periods = 20;           % fundumental periods to simulate 

% savename = 'RecurranceTest_005';   % name for data save file 
% savedir = 'C:\Users\uqrhar19\OneDrive - The University of Queensland\Documents\PhD\Hydrodynamics\Simulations\Data\Test\';
% amp     = 0.1;          % amplitude 
% N       = 256;           % # of grid points 
% d      = 0.02;          % fluid depth, 0.001 = 1000:1
% periods = 40;           % fundumental periods to simulate 
% tstep=50;
% Simulation(savename,savedir,amp,N,d,periods,tstep)

%% Damping Tests
%savename = '002_DampingTest';   % name for data save file 
savedir = 'C:\Users\uqrhar19\OneDrive - The University of Queensland\Documents\PhD\Hydrodynamics\Simulations\Data\Test\';
amp     = 0.1;          % amplitude 
N       = 2048;           % # of grid points 
periods = 40;           % fundumental periods to simulate 
tstep = 100;

% Gamma = 0.002;
% aspect_ratio1 = 50;
% d1 = pi/aspect_ratio1; 
% savename = '009_DampingTest';
% Simulation(savename,savedir,amp,N,d1,periods,tstep,Gamma)
% Gamma = 0.002;
% aspect_ratio2 = 100;
% d2 = pi/aspect_ratio2; 
% savename = '010_DampingTest';
% Simulation(savename,savedir,amp,N,d2,periods,tstep,Gamma)
% Gamma = 0.002;
% aspect_ratio3 = 500;
% d3 = pi/aspect_ratio3; 
% savename = '011_DampingTest';
% Simulation(savename,savedir,amp,N,d3,periods,tstep,Gamma)
% Gamma = 0.002;
% aspect_ratio4 = 1000;
% d4 = pi/aspect_ratio4; 
% savename = '012_DampingTest';
% Simulation(savename,savedir,amp,N,d4,periods,tstep,Gamma)

Gamma = 0.0005;
aspect_ratio4 = 1000;
d4 = pi/aspect_ratio4; 
savename = '013_DampingTest';
Simulation(savename,savedir,amp,N,d4,periods,tstep,Gamma)
Gamma = 0.0;
aspect_ratio4 = 1000;
d4 = pi/aspect_ratio4; 
savename = '014_DampingTest';
Simulation(savename,savedir,amp,N,d4,periods,tstep,Gamma)



%% 2D surface Plot Testing

%Simulation(savename,savedir,amp,N,d,periods)
% savename = '00_Testing_timestep2';   % name for data save file 
% d    = 0.001;          % fluid depth, 0.001 = 1000:1
% N    = 512;           % # of grid points 
% amp = 0.20;
% periods = 5;
% tstep = 200;
% Simulation(savename,savedir,amp,N,d,periods,tstep)



%% Soliton Scaling
%savename = 'QuickTest';   % name for data save file 
% savedir = 'C:\Users\s4606711\OneDrive - The University of Queensland\PhD\Hydrodynamics\Simulations\Data\SolitonSweep\';
% amp     = 0.1;          % amplitude 
% %N       = 64;           % # of grid points 
% periods = 20;           % fundumental periods to simulate 
% Simulation(savename,savedir,amp,N,d,periods)
% 
% savename = '01_depth_0p5';   % name for data save file 
% d      = 0.5;          % fluid depth, 0.001 = 1000:1
% N       = 512;           % # of grid points 
% Simulation(savename,savedir,amp,N,d,periods)
% 
% savename = '02_depth_0p1';   % name for data save file 
% d      = 0.1;          % fluid depth, 0.001 = 1000:1
% N       = 512;           % # of grid points 
% Simulation(savename,savedir,amp,N,d,periods)
% 
% savename = '03_depth_0p05';   % name for data save file 
% d      = 0.05;          % fluid depth, 0.001 = 1000:1
% N       = 1024;           % # of grid points 
% Simulation(savename,savedir,amp,N,d,periods)
% 
% savename = '04_depth_0p01';   % name for data save file 
% d      = 0.01;          % fluid depth, 0.001 = 1000:1
% N       = 1024;           % # of grid points 
% Simulation(savename,savedir,amp,N,d,periods)
% 
% savename = '05_depth_0p005';   % name for data save file 
% d      = 0.005;          % fluid depth, 0.001 = 1000:1
% N       = 2048;           % # of grid points 
% Simulation(savename,savedir,amp,N,d,periods)
% 
% savename = '07_depth_0p001';   % name for data save file 
% d    = 0.001;          % fluid depth, 0.001 = 1000:1
% N    = 2048;           % # of grid points 
% periods = 20;
% tstep = 50;
% amp = 0.1;
% Simulation(savename,savedir,amp,N,d,periods,tstep)

%% Recurrance Time

% savedir = 'C:\Users\uqrhar19\OneDrive - The University of Queensland\Documents\PhD\Hydrodynamics\Simulations\Data\RecurranceSweep\LongTraces\';
% periods = 200;
% tstep = 50;
% amp = 0.1;
% Simulation('recurrence_001_p300',savedir,amp,256,0.30,periods,tstep)
% Simulation('recurrence_002_p200',savedir,amp,256,0.20,periods,tstep)
% Simulation('recurrence_003_p100',savedir,amp,256,0.10,periods,tstep)
% Simulation('recurrence_004_p090',savedir,amp,512,0.09,periods,tstep)
% Simulation('recurrence_005_p080',savedir,amp,512,0.08,periods,tstep)
% Simulation('recurrence_006_p070',savedir,amp,512,0.07,periods,tstep)
% Simulation('recurrence_007_p060',savedir,amp,512,0.06,periods,tstep)
% Simulation('recurrence_008_p050',savedir,amp,512,0.05,periods,tstep)
% Simulation('recurrence_009_p040',savedir,amp,512,0.04,periods,tstep)
% Simulation('recurrence_010_p030',savedir,amp,512,0.03,periods,tstep)
% Simulation('recurrence_011_p020',savedir,amp,1024,0.02,periods,tstep)
% Simulation('recurrence_012_p010',savedir,amp,1024,0.01,periods,tstep)
% Simulation('recurrence_013_p008',savedir,amp,2048,0.008,periods,tstep)
% Simulation('recurrence_014_p006',savedir,amp,2048,0.006,periods,tstep)
% Simulation('recurrence_015_p004',savedir,amp,2048,0.004,periods,tstep)
% Simulation('recurrence_015_p002',savedir,amp,4096,0.002,periods,tstep)























