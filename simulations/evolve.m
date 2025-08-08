%% evolve.m
% Defines the differential equation and derivatives at each time step for
% the evolution of the surface of the fluid. 

function dQdt = evolve(t,Q,k,N,d,kappa,Gamma,visc) %#ok
%% Derivatives
 
x = Q(1:N).';  
x(1) = -pi; 
x(N/2+1) = 0; 
x(N/2+2:N) = - x(N/2:-1:2);

x = project(x - (0:N-1)/N*2*pi ) + (0:N-1)/N*2*pi;

y = Q(N+1:2*N).'; 
y(N/2+2:N) = fliplr(y(2:N/2));
y = project(y);

phi = Q(2*N+1:3*N).'; 
phi(N/2+2:N) = fliplr(phi(2:N/2));
phi = project(phi);

if (any(diff(x) <0) )
    error('Non-monatonic x. Aborting.')
end

%x(1) = 0; x(N/2+1) = pi; x(N/2+2:end) = 2*pi-fliplr(x(2:N/2));

Z = x+1i*y;
ZF = fft(Z - (0:N-1)*(2*pi/N));

%k = fftshift((-N/2:N/2-1)*2*pi/N);
ZD = 1i*k.*ZF; 
ZD(1) = 2*pi;
ZD = ifft(ZD);

ZDD = -k.^2.*ZF; 
ZDD = ifft(ZDD);

%derivative of phi
dphi_dj = real(ifft(1i*k.*fft(phi)));

%Matrices for vortex strengths
K = matrixD(Z,ZD,ZDD,N,d);

%solve vortex strengths and derivatives
a = (K\dphi_dj.').';
ad = real(ifft(1i*k.*fft(a)));

%calculate velocities
w = velocity(a,ad,Z,ZD,ZDD,N,d);
u = real(w); v = -imag(w);

%u = project(u); v = project(v);

Xd = real(ZD);
Yd = imag(ZD);
Xdd = real(ZDD);
Ydd = imag(ZDD);
Rs_inv = (Xd.*Ydd - Yd.*Xdd)./(Xd.^2 + Yd.^2).^(3/2);

%% Equation of motion for phi. 
%switch simparams.waveType
%    case 'VdW'
        dphi_dt = d/3*(1./(1+y/d).^3 - 1) ... % van der Waals restoring force
                + (u.^2 + v.^2)/2 ... % velocities, x and y
                + kappa*Rs_inv ... % surface tension
                - Gamma*phi; % damping rate
                - visc.*(ifft(1i*k.*fft(u)))*(Xd'); % attempt at viscocity parameter

%    case 'gravity'
%        dphi_dt = -y ...
%                  + (u.^2 + v.^2)/2 ...
%                  + kappa*Rs_inv;
%end

%u(1) = 0; u(N/2+1) = 0;
dQdt = [u.';v.';dphi_dt.'];

%% (OLD)Equation of motion for phi. 
%switch simparams.waveType
%    case 'VdW'
%         dphi_dt = d/3*(1./(1+y/d).^3 - 1) ... % van der Waals restoring force
%                 + (u.^2 + v.^2)/2 ... % velocities, x and y
%                 + kappa*Rs_inv); % surface tension
%
%    case 'gravity'
%        dphi_dt = -y ...
%                  + (u.^2 + v.^2)/2 ...
%                  + kappa*Rs_inv;
%end

%u(1) = 0; u(N/2+1) = 0;
dQdt = [u.';v.';dphi_dt.'];
end
%% stability
function x = project(x)
 %p = ones(size(x));
 N = length(x);
 
 %The analysis in Roberts shows that the highest mode needs to be zeroed
 %for the method to be stable. For large amplitude dynamics, instability 
 %develops, but exponential filter [Wilkening (2021)] restores stability ... 
 k = -N/2:N/2-1;
 p = exp(-36*(k/(N/2)).^36);
 p = fftshift(p);
 x = ifft(p.*fft(x));
 %keyboard
end
 


%% Edited Code
% %% evolve.m
% % Defines the differential equation and derivatives at each time step for
% % the evolution of the surface of the fluid. 
% 
% function dQdt = evolve(t,Q,k,N,d,kappa,Gamma,visc) %#ok
% %% Derivatives
% 
% x = Q(1:N).';  
% %x(1)
% 
% x(N/2+1) = 0; 
% x(N/2+2:N) = - x(N/2:-1:2);
% 
% x = project(x - (0:N-1)/N*2*pi ) + (0:N-1)/N*2*pi;
% 
% y = Q(N+1:2*N).'; 
% y(N/2+2:N) = fliplr(y(2:N/2));
% y = project(y);
% 
% phi = Q(2*N+1:3*N).'; 
% phi(N/2+2:N) = fliplr(phi(2:N/2));
% phi = project(phi);
% 
% if (any(diff(x) <0) )
%     error('Non-monatonic x. Aborting.')
% end
% 
% %x(1) = 0; x(N/2+1) = pi; x(N/2+2:end) = 2*pi-fliplr(x(2:N/2));
% 
% Z = x+1i*y;
% ZF = fft(Z - (0:N-1)*(2*pi/N));
% 
% %k = fftshift((-N/2:N/2-1)*2*pi/N);
% ZD = 1i*k.*ZF; 
% ZD(1) = 2*pi;
% ZD = ifft(ZD);
% 
% ZDD = -k.^2.*ZF; 
% ZDD = ifft(ZDD);
% 
% %derivative of phi
% dphi_dj = real(ifft(1i*k.*fft(phi)));
% 
% %Matrices for vortex strengths
% K = matrixD(Z,ZD,ZDD,N,d);
% 
% %solve vortex strengths and derivatives
% a = (K\dphi_dj.').';
% ad = real(ifft(1i*k.*fft(a)));
% 
% %calculate velocities
% w = velocity(a,ad,Z,ZD,ZDD,N,d);
% u = real(w); v = -imag(w);
% 
% %u = project(u); v = project(v);
% 
% Xd = real(ZD);
% Yd = imag(ZD);
% Xdd = real(ZDD);
% Ydd = imag(ZDD);
% Rs_inv = (Xd.*Ydd - Yd.*Xdd)./(Xd.^2 + Yd.^2).^(3/2);
% 
% %% Equation of motion for phi. 
% %switch simparams.waveType
% %    case 'VdW'
%         dphi_dt = d/3*(1./(1+y/d).^3 - 1) ... % van der Waals restoring force
%                 + (u.^2 + v.^2)/2 ... % velocities, x and y
%                 + kappa*Rs_inv;% ... % surface tension
%                 %- Gamma*phi; % damping rate
%                 %- visc.*(ifft(1i*k.*fft(u)))*(Xd'); % attempt at viscocity parameter
% 
% %    case 'gravity'
% %        dphi_dt = -y ...
% %                  + (u.^2 + v.^2)/2 ...
% %                  + kappa*Rs_inv;
% %end
% 
% u(1) = 0; u(N/2+1) = 0;
% %u(1) = cos(t); u(N/2+1) = 0;
% dQdt = [u.';v.';dphi_dt.'];
% 
% %% (OLD)Equation of motion for phi. 
% %switch simparams.waveType
% %    case 'VdW'
% %         dphi_dt = d/3*(1./(1+y/d).^3 - 1) ... % van der Waals restoring force
% %                 + (u.^2 + v.^2)/2 ... % velocities, x and y
% %                 + kappa*Rs_inv); % surface tension
% %
% %    case 'gravity'
% %        dphi_dt = -y ...
% %                  + (u.^2 + v.^2)/2 ...
% %                  + kappa*Rs_inv;
% %end
% 
% %u(1) = 0; u(N/2+1) = 0;
% dQdt = [u.';v.';dphi_dt.'];
% end
% %% stability
% function x = project(x)
%  %p = ones(size(x));
%  N = length(x);
%  
%  %The analysis in Roberts shows that the highest mode needs to be zeroed
%  %for the method to be stable. For large amplitude dynamics, instability 
%  %develops, but exponential filter [Wilkening (2021)] restores stability ... 
% %  k = -N/2:N/2-1;
% %  p = exp(-36*(k/(N/2)).^36);
% %  p = fftshift(p);
% %  x = ifft(p.*fft(x));
%  %keyboard
% end
%  
% 
