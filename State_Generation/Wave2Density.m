function [ rho ] = Wave2Density( wave,p )
%WAVE2DENSITY transform an ensemble of wavefunctions wave=[w_1,w_2,...]
%with probabilities p=[p_1,p_2,...] into a density matrix rho
if min(size(wave))==1
    if size(wave,1)==1
        wave=conj(wave');
    end
end

if nargin==1
    p=ones(size(wave,2),1)/size(wave,2);
end

len=size(wave,1);
rho=zeros(len,len);

for i=1:size(wave,2)
    rho2=wave(:,i)*wave(:,i)';
    rho=rho+rho2/norm(wave(:,i))*p;
end
if sum(p)~=1
    rho=rho/sum(p);
end
end

