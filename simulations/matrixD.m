%% Original Code
% function K = matrixD(Z,ZD,ZDD,N,d)
% K = zeros(N);
% %todo: vectorize this 
% f = 0.25/pi;
% for aa = 1:N
%      bb = (1:N ~= aa);
%      K(aa,aa) =     0.5 + f*imag(ZDD(aa)/ZD(aa));
%      K(aa,bb) =         + f*imag(ZD(aa)*cot( 0.5*(Z(aa)-Z(bb))) );   
%      K(aa,1:N)  = K(aa,1:N) - f*imag(ZD(aa)*cot( 0.5*(Z(aa) - conj(Z) + 2i*d) ));
% end
% end

%% Vectorized Code

function K = matrixD(Z,ZD,ZDD,N,d)
% K = zeros(N);
% %todo: vectorize this 
% f = 0.25/pi;
% for aa = 1:N
%      bb = (1:N ~= aa);
%      K(aa,aa) =     0.5 + f*imag(ZDD(aa)/ZD(aa));
%      K(aa,bb) =         + f*imag(ZD(aa)*cot( 0.5*(Z(aa)-Z(bb))) );   
%      K(aa,1:N)  = K(aa,1:N) - f*imag(ZD(aa)*cot( 0.5*(Z(aa) - conj(Z) + 2i*d) ));
% end

f = 0.25/pi;
%Get the off diagonal terms (K(aa,bb) terms from above)
%      Note: diagonals are infinities: removing them now is actually slower I
%             think? Can set them to zero here and then fill them in below if desired
K = cot(0.5*(Z.'-Z)); %Difference matrix by broadcasting, nonconjugate transpose is used. Infinities at diagonals though
% K = 1./(tan(0.5*(Z.'-Z)));
K = f*imag(ZD.'.*K); %row-wise (columnwise?) multiplication of the nondiagonal elements

%Add the diagonal terms back in (K(aa,aa) line above)
K(1:N+1:end) = 0.5+f*imag(ZDD./ZD); %To do: can do this at the same time as removing diagonals?

%Now there is that full matrix row-wise (column-wise?) term: (K(aa,1:N) line above)
K = K -f*imag(ZD.'.*cot(0.5*(Z.'-conj(Z) + 2i*d)));
% K = K -f*imag(ZD.'.*(1./(tan(0.5*(Z.'-conj(Z)+2i*d)))));

end