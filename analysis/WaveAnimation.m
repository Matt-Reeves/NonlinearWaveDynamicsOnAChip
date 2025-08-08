%% WaveAnimation.m
% Plots movie of the wave profile over time
function WaveAnimation(data)


%%
N = data.simparams.N;
Nt = data.simparams.Nt;
time = (1:Nt)*data.simparams.DT;

x = squeeze(data.solution(:,1,:));
y = squeeze(data.solution(:,2,:));


ymax = max(y(:));
ymin = min(y(:));

depthscale = (100/pi)*data.simparams.d; % scale in microns of the depth of the fluid
depthscale = round(depthscale*1000);
%% Calculate conserved quantities from the solution
volumeFlux = zeros(1,round(Nt));
eta_bar = zeros(1,round(Nt));

for jj = 1:Nt
    temp = data.solution(:,:,jj);
    temp = temp(:);
    %tempgrids.Q = temp;
    %tempgrids.k = simgrids.k;
% 
%     % quantities = conservedQuantities(tempgrids,simparams);
%     % energy.E(jj)   = quantities.E;
%     % energy.T(jj)   = quantities.T;
%     % energy.V(jj)   = quantities.V;
%     % volumeFlux(jj) = quantities.volumeFlux;
%     % eta_bar(jj)    = quantities.Ybar;
 end
% %% Plot dynamics
% savefile = strcat(savedir,savename);
% save(savefile)
% v = VideoWriter(savefile,'MPEG-4');
% v.FrameRate = 10;
% open(v);
% disp('Compiling video of wave simulation')
% tic
% for jj = 1:Nt
% 
%     figure(999)
%     clf
% 
%     % Waveprofile Plots
%     %subplot(121)
%     plot(x(:,jj),y(:,jj),'b')
%        ylim([-simparams.d 2*simparams.A])%ylim([-2 2]*simparams.A)
%        yticks([-simparams.d  0]), yticklabels({'0','d'})
%        xlim([-pi 0])
%        %xline(-(pi-0.06*pi),'--k') %Cavity Marker
%        %xline(-(pi-Ntrac),'--r') % tracking the point on the surface I am equating to the mode
%        xline(-(pi-0.045*pi),'--k'), xline(-(pi-0.075*pi),'--k')
%        xticks([-pi -(pi-0.06*pi) -3*pi/4 -pi/2 -pi/4 0]), xticklabels({'0','cav','25','50','75','100'})
%        ylabel('Height'), xlabel('x position (\mu m)')
%        title(['d = ', num2str(depthscale),' nm' ,',   N = ', num2str(simparams.N),',  \Gamma = ',num2str(simparams.Gamma),',  \mu = ',num2str(simparams.visc)])
% 
%     % Energy Plots
% %     subplot(122)
% %     plot(time(1:jj),energy.E(1:jj))
% %     hold on
% %     plot(time(1:jj),energy.V(1:jj))
% %     plot(time(1:jj),energy.T(1:jj))
% %     xlim([0 time(end)])
% %     ylim([0 max(energy.E)*1.2])
% 
%     frame = getframe(gcf);
%     writeVideo(v,frame);
% 
% end
% close(v)
% toc
% disp(['Video compiled, saved as ',savename,'.MP4'])
% 
% %%
% % figure()
% % subplot(211)
% % plot(time,y(3,:))
% % ylabel('Height')
% % title('Driven Response (\omega=0.145)')
% % subplot(212)
% % plot(time,energy.E)
% % ylabel('Energy')
% 
% %% Save Video Output
% % save(['d' num2str(d) '_amp' num2str(ampvals) '.mat'])
% %saveas(gcf,'DrivenResponsew145.png')