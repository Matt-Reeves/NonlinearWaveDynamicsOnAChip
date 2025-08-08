%% SurfacePlot.m
% plots the wave profile in a 2D plot with time, to see the solotions
% fission and evolution 

function SurfacePlot(simparams,solution,savename,savedir)

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
% Spatial Averaging
 xspace = squeeze(solution.Q(:,1,:));
 myindex = (xspace>(-3.002))&(xspace<(-2.9060));
 yspace = squeeze(solution.Q(:,2,:));
 temporary = yspace;
 temporary(~myindex)=0;
 mymean = sum(temporary,1)./sum(myindex,1);
 
 % Surface Profile 
clear surface
start = 1;
point = 10;
finish = simparams.Nt;
surface = zeros(simparams.N/2,round(simparams.Nt-start));
for i = start:simparams.Nt
surface(:,i-start+1) = solution.Q(1:simparams.N/2,2,i);
end
vertsurface = surface';
%surface = rescale(surface).*255;

 figure
 subplot(2,1,1)
        imagesc(surface(:,start:finish))
        hold on
        yline(Nopt,'--k','LineWidth',2) %Cavity Marker
        colormap winter;
        title(['d = ', num2str(depthscale),' nm' ,',   N = ', num2str(simparams.N),',  \Gamma = ',num2str(simparams.Gamma),',  \mu = ',num2str(simparams.visc)])
        yticks([simparams.N/simparams.N Nopt simparams.N/8 simparams.N/4 3*simparams.N/8 simparams.N/2])
        yticklabels({'0','Cav','25','50','75','100'}), ylabel('Position x (\mu m)')
        xlabel('time ')
        %colorbar('Ticks',[0],'TickLabels',{'\eta -d'})
        hold off
 subplot(2,1,2)
        plot(mymean(start:finish),'-k')
       hold on
       xlim([start finish])
        ylabel('Amp @ Cav'), xlabel('Time')%, title('Time trace from wave profile')
        hold off
        
%         figure
%  subplot(1,2,1)
%         imagesc(vertsurface(start:finish,:))
%         hold on
%         yline(start+point,'--k','LineWidth',2) %Cavity Marker
%         colormap winter;
%         title(['d = ', num2str(depthscale),' nm' ,',   N = ', num2str(simparams.N),',  \Gamma = ',num2str(simparams.Gamma),',  \mu = ',num2str(simparams.visc)])
%         xticks([simparams.N/simparams.N Nopt simparams.N/8 simparams.N/4 3*simparams.N/8 simparams.N/2])
%         xticklabels({'0','Cav','25','50','75','100'}), ylabel('Position x (\mu m)')
%         ylabel('time ')
%         colorbar('Ticks',[0],'TickLabels',{'\eta -d'})
%         hold off
%  subplot(1,2,2)
%         plot(solution.Q(1:simparams.N/2,1,start+point),solution.Q(1:simparams.N/2,2,start+point),'.-k')
%        
%         hold on
%         ylabel('Amplitude'), xlabel('Position')%, title('Time trace from wave profile')
%         hold off

%% Saveing
        
        
        