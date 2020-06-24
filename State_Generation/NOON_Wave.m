function [ phi ] = NOON_Wave( n,N,m,theta )
%NOON_WAVE generates a NOON state for n photons in N modes with the photons
%in an equal superposition of n photons in mode m or N-m with a relative
%phase of N*theta
if nargin<2
    error('Not enough input arguments')
elseif nargin==2
    m=1;
    theta=0;
elseif nargin==3
    theta=0;
end
if m>N
    m=mod(m,N)+1;
end

weights=1/sqrt(2)*[1 exp(1i*theta*n)];
fock=zeros(2,N);
fock(1,m)=n;
fock(2,N-m+1)=n;
phi=fock2vec(n,N,fock,weights);
end

