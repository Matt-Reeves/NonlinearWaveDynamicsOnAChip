function [residuals,solution] = NewtonRaphson(simgrids,simparams,tol)

N = simparams.N;
residuals = calculateResiduals(simgrids,simparams);

figure(1)
hold on
plot(simgrids.Q(1:N/2+1),simgrids.Q((1:N/2+1)+N),'.-b')

while residuals.SSE > tol


Jacobi = calculateJacobian(simgrids,simparams,residuals);
delta = -4/7*(Jacobi.'\residuals.all);
[x,y,phi] = unpack(simgrids.Q);
deltaX = delta(1:simparams.N/2-1);
deltaY = delta(simparams.N/2:simparams.N);
deltaT = delta(end-1);
deltaAlpha = delta(end);

xtest = x(1:simparams.N/2+1) + [0 deltaX.' 0];
ytest = y(1:simparams.N/2+1) + deltaY.';
[xtest,ytest,phi] = imposeSymmetry(xtest,ytest,phi,simparams.N);
Q = pack(xtest,ytest,phi);
simgrids.Q = Q;

simparams.T = simparams.T + deltaT;
simparams.tfin   = simparams.T*0.25; 
simparams.DT     = simparams.tfin;
simparams.dt     = simparams.DT/300;
simparams.Nt     = (simparams.tfin/simparams.DT);
simparams.nt     = (simparams.DT/simparams.dt);

simparams.alpha = simparams.alpha + deltaAlpha;

residuals = calculateResiduals(simgrids,simparams)

figure(1)
plot(x(1:N/2+1)+[0 deltaX.', 0],y(1:N/2+1)+deltaY.','.-')
hold on
drawnow
figure(2)
semilogy(abs(residuals.all))
hold on
drawnow

end


residuals = calculateResiduals(simgrids,simparams);
solution = struct; 
solution.Q = simgrids.Q;
solution.x = xtest;
solution.y = ytest;
solution.T = simparams.T;
solution.alpha = simparams.alpha;
end

function residuals = calculateResiduals(simgrids,simparams)

eta_bar = meanHeight(simgrids,simparams);
Ac = getCrestAcceleration(simgrids,simparams);
%x-residuals

tempgrid = simgrids;
for jj = 1:simparams.Nt
    for kk = 1:simparams.nt
    tempgrid.Q = rk4(tempgrid,simparams);
    end
%     plot(simgrids.Q(1:simparams.N),simgrids.Q((1:simparams.N)+simparams.N),'.-')
%     ylim([-2 2]*simparams.A)
%     drawnow
end

%% Residuals at T/4 guess
N = simparams.N; 
Q = tempgrid.Q;
x = Q(1:N).'; y = Q(N+1:2*N).'; phi = Q(2*N+1:3*N).';
resx = x(2:N/2) - 2*pi/N*(1:N/2-1);

%y-residuals
j = 1:N/4; 
resy = y(j) - y(N/2+2-j);

%phi-residuals
j = 1:N/4+1;
resphi = phi(j) + phi(N/2+2-j) - simparams.alpha;

%constructs the residuals vector of the N+2 unknowns
residuals.all = [eta_bar; simparams.Act+Ac; resx.'; resy.'; resphi.' ];
residuals.eta_bar = eta_bar;
residuals.Ac = Ac;
residuals.resx = resx.'; 
residuals.resy = resy.';
residuals.phi = resphi.';
residuals.SSE = sum(residuals.all.^2);

end
function eta_bar = meanHeight(simgrids,simparams)
N = simparams.N;
k = simgrids.k;
x = simgrids.Q(1:N).'; y = simgrids.Q(N+1:2*N).';
xF = 1i*k.*fft(x - (0:N-1)*(2*pi/N)); xF(1) = 2*pi; dxdj = ifft(xF);
eta_bar = sum(y.*dxdj);
end
function [x,y,phi] = unpack(Q)
 N = length(Q)/3;
 x = Q(1:N).'; 
 y = Q(N+1:2*N).';
 phi = Q(2*N+1:3*N).'; 
end
function Q = pack(x,y,phi)
Q = [x(:);y(:);phi(:)]; 
end
%{ 
function Jacobi = calculateJacobian(simgrids,simparams,residuals)
N = simparams.N;

[x0,y0,phi0] = unpack(simgrids.Q); 

%% x-part of the Jacobian
disp('Calculating X-Jacobian')
x = x0; y= y0; phi = phi0;
temp = simgrids;
delta = 0.01*2*pi/N;
xvals = x(2:N/2);
Jacobi_x = zeros(N/2-1,N+2);

for kk = 1:length(xvals)
   
    xtemp = xvals;
    xtemp(kk) = xtemp(kk)+delta;
    x = [0, xtemp, pi, 2*pi-fliplr(xtemp)];
    [x,y,phi] = imposeSymmetry(x,y,phi,N);
    temp.Q = pack(x,y,phi);
   
    residuals_new = calculateResiduals(temp,simparams);
    Jacobi_x(kk,:) = (residuals_new.all - residuals.all)/delta;
end
%% y-part of the Jacobian
disp('Calculating Y-Jacobian')
x= x0; phi = phi0;
temp = simgrids; 

Jacobi_y = zeros(N/2+1,N+2);
for kk = 1:N/2+1
    delta = 0.01*simparams.A;
    y = y0;
    y(kk) = y(kk) + delta; 
    [x,y,phi] = imposeSymmetry(x,y,phi,N);
    temp = simgrids;
    temp.Q = pack(x,y,phi);
    residuals_new = calculateResiduals(temp,simparams);
    Jacobi_y(kk,:) = (residuals_new.all - residuals.all)/delta;
end
%% T-part
x= x0; y =y0; phi = phi0;
simgrids.Q = pack(x,y,phi); 
temp = simparams; 
delta = 0.005*temp.tfin; 
temp.tfin = temp.tfin + delta;
temp.DT = temp.tfin/20;
temp.dt = 0.05*temp.DT;
temp.nt = (temp.DT/temp.dt);
temp.Nt = (temp.tfin/temp.DT);
residuals_new = calculateResiduals(simgrids,temp);
Jacobi_T(1,1:N+2) = (residuals_new.all - residuals.all)/delta;
%% alpha part
x= x0; y =y0; phi = phi0;
Q = pack(x,y,phi); 
delta = 0.001;
temp = simparams;
temp.alpha = temp.alpha + delta; 

residuals_new = calculateResiduals(simgrids,temp);
Jacobi_alpha(1,1:N+2) = (residuals_new.all - residuals.all)/delta;

Jacobi = [Jacobi_x; Jacobi_y; Jacobi_T; Jacobi_alpha];

 

end
%}

function Jacobi = calculateJacobian(simgrids,simparams,residuals)
N = simparams.N;

[x0,y0,phi0] = unpack(simgrids.Q); 

%% x-part of the Jacobian
disp('Calculating X-Jacobian')

x = x0;
change = 0.01;
delta = change*2*pi/N;
xvals = x(2:N/2);
Jacobi_x = zeros(N/2-1,N+2);
old_residuals = residuals.all;

parfor kk = 1:length(xvals)
   
    xtemp = xvals;
    xtemp(kk) = xtemp(kk)+delta;
    x = [0, xtemp, pi, 2*pi-fliplr(xtemp)]; y = y0; phi = phi0;
    [x,y,phi] = imposeSymmetry(x,y,phi,N);
    
    temp = simgrids;
    temp.Q = pack(x,y,phi);
   
    residuals_new = calculateResiduals(temp,simparams);
    Jacobi_x(kk,:) = (residuals_new.all - old_residuals)/delta;
end
%% y-part of the Jacobian
disp('Calculating Y-Jacobian')

Jacobi_y = zeros(N/2+1,N+2);
A = max(y0);
parfor (kk = 1:N/2+1,8)
    delta = change*A;
    x= x0; y = y0; phi = phi0;
    y(kk) = y(kk) + delta; 
    [x,y,phi] = imposeSymmetry(x,y,phi,N);
    temp = simgrids;
    temp.Q = pack(x,y,phi);
    residuals_new = calculateResiduals(temp,simparams);
    Jacobi_y(kk,:) = (residuals_new.all - old_residuals)/delta;
end
%% T-part
x= x0; y =y0; phi = phi0;
simgrids.Q = pack(x,y,phi); 
temp = simparams; 
delta = change*simparams.T; 
temp.tfin = temp.tfin + delta;
temp.DT = temp.tfin;
temp.dt = temp.DT/400;
temp.nt = (temp.DT/temp.dt);
temp.Nt = (temp.tfin/temp.DT);
residuals_new = calculateResiduals(simgrids,temp);
Jacobi_T(1,1:N+2) = (residuals_new.all - old_residuals)/delta;
Jacobi_T(1,1:2) = 0;
%% alpha part

delta = max(sqrt(eps),change*simparams.alpha);
temp = simparams;
temp.alpha = temp.alpha + delta; 

residuals_new = calculateResiduals(simgrids,temp);
Jacobi_alpha(1,1:N+2) = (residuals_new.all - old_residuals)/delta;
Jacobi_alpha(1,1:2) =0;

Jacobi = [Jacobi_x; Jacobi_y; Jacobi_T; Jacobi_alpha];


end