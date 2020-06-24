function [ phi ] = NOON_Wave_XY( N,theta )
%NOON_WAVE generates a NOON state for N photons in 2 modes with the photons
%in an equal superposition of N photons 
if nargin<1
    error('Not enough input arguments')
elseif nargin==1
    theta=0;
end

index=2.^(0:N)+1;
weights=1/sqrt(2)*[1 exp(1i*theta)];
phi=zeros(2^(N+1),1);
phi(index(1))=weights(1);
phi(index(end))=weights(2);
end