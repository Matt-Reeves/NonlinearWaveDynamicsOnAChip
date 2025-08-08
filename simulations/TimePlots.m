%% Simulation Plots
function TimePlots(solution,simparams,tempgrid,energy,savename,savedir)
% Time plots, plots values of interest over time, such as time traces and
% energy conservation during the simulation. 
%% Plot Testing


%% Plotting scale Math
% for finding point of optical cavity along wavetank.  
depthscale = (100/pi)*simparams.d; % scale in microns of the depth of the fluid
depthscale = round(depthscale*1000);

Nopt = 0.06*simparams.N/2; %optical mode location % from end
Nopt = round(Nopt);
OptModeW = 0.03;
Nopt1 = (Nopt-OptModeW/2)*simparams.N; % start of optical mode
Nopt2 = (Nopt+OptModeW/2)*simparams.N;  % start of optical mode
Nopt1 = round(Nopt1);% Q index for start of mode
Nopt2 = round(Nopt2); % Q index for end of mode

Ntrac = (2*pi*Nopt)/simparams.N; % test of which point is being plotted compared to where the opticalmode should be. 
%% Spatial Averaging
 xspace = squeeze(solution.Q(:,1,:));
 myindex = (xspace>(-3.002))&(xspace<(-2.9060));
 yspace = squeeze(solution.Q(:,2,:));
 temporary = yspace;
 temporary(~myindex)=0;
 mymean = sum(temporary,1)./sum(myindex,1);
 
 %figure()
 %plot(mymean,'-b')
%% Time Trace Plots

figure
   % Wave Profile Plotting
    subplot(3,1,1)
        plot(tempgrid.Q(1:simparams.N),tempgrid.Q((1:simparams.N)+simparams.N),'.-')
        hold on
        ylim([-simparams.d 2*simparams.A])%ylim([-2 2]*simparams.A)
        yticks([-simparams.d  0]), yticklabels({'0','d'})
        xlim([-pi 0])
        xline(-(pi-0.06*pi),'--k') %Cavity Marker
        xline(-(pi-Ntrac),'--r') % tracking the point on the surface I am equating to the mode
        xline(-(pi-0.045*pi),'--b'), xline(-(pi-0.075*pi),'--b') % optical mode width
        xticks([-pi -(pi-0.06*pi) -pi/2 0])
        xticklabels({'0','cav','50','100'})
        ylabel('Height'), xlabel('x position (\mu m)'), title('Fluid profile at end')
        legend('Wave','Mode Center','Plot Point','Mode')
    hold off

    % Energy Conservation Plotting
    subplot(3,1,2)
        plot(energy.E)
        hold on
        plot(energy.T)
        plot(energy.V)
        legend('E','T','V'), ylabel('Energy'), xlabel('Time')%, title('Energy')
        title(['d = ', num2str(depthscale),' nm' ,',   N = ', num2str(simparams.N),',  \Gamma = ',num2str(simparams.Gamma),',  \mu = ',num2str(simparams.visc)])
    hold off
    
    % Time trace of amplitude at cavity plotting
    subplot(3,1,3)
        plot(squeeze(solution.Q(Nopt,2,:)),'k') % history of height of in centre of optical mode
        hold on
        plot(mymean,'-b')
        %ylim([-1.2 1.2]*simparams.A) % acrage over optical mode
        %plot(squeeze(solution.Q(1,2,:)),'b') % history of height of film at point
        ylabel('Amp @ Cav'), xlabel('Time')%, title('Time trace from wave profile')
    hold off
    
%     % Plot of acerage film heigher over width of opticalmode
%     subplot(4,1,4)
%         plot(mymean,'-b')
%     hold off
