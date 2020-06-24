function [ phi ] = NOON_Wave_Ladder( n,N,m,theta )
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
fock=zeros(2,(n+1)*N);
fock(1,m*(n+1))=1;
for i=1:N
    if i~=m
        fock(1,(i-1)*(n+1)+1)=1;
    end
end
fock(2,(N-m+1)*(n+1))=1;
for i=1:N
    if i~=(N-m+1)
        fock(2,(i-1)*(n+1)+1)=1;
    end
end
fock(:,1:end)=fock(:,end:-1:1);
phi=fock2vec(1,(n+1)*N,fock,weights);
end