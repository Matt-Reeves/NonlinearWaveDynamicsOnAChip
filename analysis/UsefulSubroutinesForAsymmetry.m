
%% Subroutines
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
xbar = mean(x);
window_weight = sum(window);
S = xbar*window_weight + conv(x-xbar,window,'same');
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