%% simgrids_initialize.m
% initializes the the aray fro tracking the surface of the wave during the
% simulation. 

function simgrids = simgrids_initialize(simparams)

N = simparams.N; 

J = (0:N-1);
x0 = J*(2*pi/N) - pi;
y0 = -simparams.A*cos(x0);
phi0 = zeros(size(x0));

simgrids = struct;
simgrids.Q = [x0(:);y0(:);phi0(:)];
simgrids.v = [x0(2:N/2).'; y0(1:N/2+1).'; simparams.T; simparams.alpha ];
clear x0 y0 phi0

%The highest fourier mode must be supressed [see Roberts, Trefethen]. 
k = ((-N/2:N/2-1)*2*pi/N);
k(1) =0;
k = fftshift(k);
simgrids.k = k; 
clear k; 

return
