%% Recurrance Analysis
% for loading and plotting data from Matt Hysdrodynbamics simulations. 
%% Index of subroutines
% [01] simdata = DataImport(datafile, datadir, lscale)
% [02] SurfacePlot(simdata)
% [03] TimeTracePlot(simdata,n)
% [04] ModeAvgPlot(simdata)
% [05] surface = SurfaceEvolution(solution,simparams)
% [06] AssymetryPlot(simdata,window)
% [07] aspect_ratio = AspectRatio(simparams)
% [08] depth  = ScaleDepth(simparams,lscale)
% [09] time_axis = ScaleTime(solution,simparams,lscale)
% [10] spatial_average = SpatialAverage(simparams,simdata.surface)
% [W] Walter Subroutines (for assymetry) 

%% Variables 
% clear simparams tempgrid DataName solution volumeFlux Ac amp amplitude ans d energy eta_bar Gamma jj kappa kk N omega0 periods quantities simgrids visc volumeFlux waveType
clear all

% lambda is the fundumental resonator wavelength
lambda = 200*10^(-6); % 200 um for 1D crystal
% lscale is half the fundumental wavelength

%% Load data
% Select Directory
directory = 'C:\Users\uqrhar19\OneDrive - The University of Queensland\Documents\PhD\Hydrodynamics\Simulations\Data\RecurranceSweep';
folder = 'LongTraces';
% "C:\Users\uqrhar19\OneDrive - The University of Queensland\Documents\PhD\Hydrodynamics\Simulations\Data\RecurranceSweep\LongTraces\recurrence_011_p020.mat"
data = 'recurrence_011_p020.mat';
simdata = DataImport(data,folder,directory,lambda/2);

% data14 = '014_DampingTest';
% simdata14 = DataImport(data14,folder,directory,lambda/2);

%% Plots

dir = strcat(directory,'\',folder,'\');
WaveAnimation(simdata)

%% Plot for Paper
window = 300;
figure
sgtitle(['Aspect ratio = 1:',num2str(simdata.ar),', Depth = ',num2str(round(simdata.depth*10^9)),' [nm]'])
subplot(4,1,1)
SurfacePlot(simdata)
subplot(4,1,2)
ModeAvgPlot(simdata)
subplot(4,1,3)
TimeTracePlot(simdata)
subplot(4,1,4)
AssymetryPlot(simdata,window)
hold on
title('Assymetry of average')
hold off

figure
sgtitle(['\Gamma = 0, Aspect ratio = 1:',num2str(simdata14.ar),', Depth = ',num2str(round(simdata14.depth*10^9)),' [nm]'])
subplot(4,1,1)
SurfacePlot(simdata14)
subplot(4,1,2)
ModeAvgPlot(simdata14)
subplot(4,1,3)
TimeTracePlot(simdata14)
subplot(4,1,4)
AssymetryPlot(simdata14,window)
hold on
title('Assymetry of average')

%% [00] Ray Subroutines
% Below are functions for importing data from simulation and
% analysing/plotting 
%% [01]  DataImport
function simdata = DataImport(datafile,folder, datadir, lscale)
% Import the data file, assign label and create structure for required
% paramters for analysis. 
    data = strcat(datadir,'\',folder,'\',datafile);
    load(data)
    simdata.surface = SurfaceEvolution(solution,simparams);
    simdata.time = ScaleTime(solution,simparams,lscale);
    simdata.average = SpatialAverage(simparams,simdata.surface);
    simdata.ar = round(AspectRatio(simparams));
    simdata.depth = ScaleDepth(simparams,lscale);
    simdata.simparams = simparams; % sve simparams in case useful
    simdata.n = simparams.N; % points tracking the surface
    simdata.solution = solution.Q; % x, y, u, v data
end

%% [02] SurfacePlot
function SurfacePlot(simdata)
% plot the full spatial information agaisnt time (simulated time steps)
    realsurface = simdata.surface(1:simdata.simparams.N/2,:);
    Nopt = round(0.06*simdata.simparams.N);
    imagesc(realsurface)
    line([1,2000], [Nopt,Nopt], 'Color', 'b');
    % hold on
    colormap winter
    xlabel('simulation time [steps]')
    yticks([simdata.n/simdata.n simdata.n/8 simdata.n/4 3*simdata.n/8 simdata.n/2])
    yticklabels({'0','25','50','75','100'}), ylabel('Position x (\mu m)')
    % title(['Surface evolution, Aspect ratio = 1:',num2str(simdata.ar),', Depth = ',num2str(round(simdata.depth*10^9)),' [nm]'])
    % hold off
end

%% [03] TimeTracePlot
function TimeTracePlot(simdata)
% Plots a time trace from edge of resonator, n is the point on the surface
% to take the time trace from. 
Nopt = round(0.06*simdata.simparams.N);
    plot(simdata.time.steps,simdata.surface(Nopt,:))
    % hold on, box on
        ylabel('Amplitude')
        xlabel('simulation time [steps]')
        % title(['Time Trace, Aspect ratio = 1:',num2str(simdata.ar),', Depth = ',num2str(round(simdata.depth*10^9)),' [nm]'])
    % hold off
end

%% [04] ModeAvgPlot
function ModeAvgPlot(simdata)
    plot(simdata.time.steps,simdata.average)
    %hold on, box on
        ylabel('Mean Amplitude'), xlabel('simulation time [steps]'),
        % title(['Optical Mode Average, Aspect ratio = 1:',num2str(simdata.ar),', Depth = ',num2str(round(simdata.depth*10^9)),' [nm]'])
    %hold off
end 

%% [05] SurfaceEvolution
function surface = SurfaceEvolution(solution,simparams)
% formats array for sptial-temproal information of simulation
    clear surface
    surface = zeros(simparams.N,round(simparams.Nt-1));
    for i = 1:simparams.Nt
        % surface(:,i-1+1) = solution.Q(1:simparams.N,2,i);
        surface(:,i) = solution.Q(1:simparams.N,2,i);
    end
end

%% [06] AssymetryPlot
function AssymetryPlot(simdata,window)
plot(simdata.time.steps,convAsy(simdata.average,window))
% hold on, box on
ylabel('Assymetry'), xlabel('simulation time [steps]')
% title(['Aspect ratio = 1:',num2str(simdata.ar),', Depth = ',num2str(round(simdata.depth*10^9)),' [nm]'])
end
%% [07] AspectRatio
function aspect_ratio = AspectRatio(simparams)
% computes the aspect ratio of the simulation for scaling
    aspect_ratio = pi/simparams.d; % 1/(d/pi) = (\lambda/2)/depth
end

%% [08] ScaleDepth
function depth  = ScaleDepth(simparams,lscale)
% calculates scaled depth based on resontaotr length (fixed resonator length = 200 um)
    depth = (lscale)/AspectRatio(simparams);
end

%% [09] ScaleTime
function time = ScaleTime(solution,simparams,lscale)
% Scaled axis for time trace plots (fixed resonator length = 200 um)
    alpha = 3.5*10^(-24); % vdw coefficient for Si
    depth = ScaleDepth(simparams,lscale); % for depth dependance
    f_zero =@(h) sqrt((3*alpha/(h^4))*(2*pi/lscale)*tanh(2*pi*h/lscale))/(2*pi); % fundumental frequency (angular)
    tpoints = length(solution.Q(1,1,:));
    time.time = (simparams.periods/f_zero(depth))/tpoints; % simulated time per time step
    time.time_axis = linspace(0,(simparams.periods/f_zero(depth)),length(solution.Q(1,1,:)))*1000; % time vector for time trace plots
    time.steps = linspace(1,length(solution.Q(1,1,:)),length(solution.Q(1,1,:))); 
end

%% [10] SpatialAverage
function spatial_average = SpatialAverage(simparams,surface)
% Spatial average, 3% of space, 6% from edge (corresnds to location of
% optical mode vvolume). 
    Nopt = 0.06; %optical mode location %age from end
    OptModeW = 0.03; % %optical mode width % 0.003 for flat profile on entire optical mode. 
    Nopt1 = (Nopt-OptModeW/2)*simparams.N; % start of optical mode
    Nopt2 = (Nopt+OptModeW/2)*simparams.N;  % start of optical mode
    Nopt1 = round(Nopt1);% Q index for start of mode
    Nopt2 = round(Nopt2); % Q index for end of mode
    temp = size(surface);
    spatial_average = zeros(1,temp(2));
    for i = 1:temp(2)
        spatial_average(i) = mean(surface(Nopt1:Nopt2,i));
    end 
end



%% [WW] walter subroutines
%The main one of interest here is the convSkew(x,k), which forces a
%blackman window of either k or k+1 elements long (whichever is odd). The
%other ones are mostly helper functions for this.  movingSkewness is my old
%way of doing it: it forces a rectangular window though.  nonparametric
%skew is generally a touch faster/more obvious and but requires a square
%window.

%Reminder that the asymmetry is the skewness of the negative imaginary
%component of the analytic signal: asymmetry = skewness(-imag(hilbert(x))).

function Asy = convAsy(x,k)
    temp_analytical = hilbert(x-mean(x)); %The subtraction of the mean "shouldn't" matter
    Asy =  -convSkew(imag(temp_analytical),k);
end

function S = convsum(x,k)
% window = ones(k,1);
window = mywindow(k);
S = conv(x,window,'same');
end

function V = convvar(x,k)
%  This ends up being significantly bigger than movvar: why?
% window = ones(k,1);
window = mywindow(k);
% k = length(window);
window_mean = window/sum(window);
xbar = convsum(x,k); %if averaging already, won't change anything, otherwise will average
% xbar = movmean(x,k);
V = convsum((x-xbar).^2,k);

% V = convsum(x.^2,window_mean) - convsum(xbar.^2,window_mean);

end

function S = convSkew(x,k)
% k = length(window);
% M = sum(window);
window = mywindow(k);
window_mean = window./sum(window);
xbar = convsum(x,k);
xvar = convvar(x,k);
% S = convsum((x-xbar).^3,window_mean)./(xvar.^3);

% V = convvar(x,window_mean);
term1 = convsum(x.^3,k);
term2 = -3*xbar.*convsum(x.^2,k);
term3 = 2*xbar.^3;
S = (term1+term2+term3)./(xvar.^(3/2));
end

function W=mywindow(k)

if mod(floor(k),2)==0
    W=blackman(floor(k+1));
else
    W=blackman(floor(k));
end
W = W/sum(W(:));
end

function Snp = nonparametricSkew(x,k)
moving_med = movmedian(x,k);
moving_var = movvar(x,k);
moving_mean = movmean(x,k);
Snp = (moving_mean-moving_med)./(sqrt(moving_var));
end

function S = movingSkewness(x,k)
meanx = movmean(x,k);
stdx = movstd(x,k);
S = movsum((x-meanx).^3,k)./(k.*stdx.^3);
end