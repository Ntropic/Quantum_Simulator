function [ phi ] = RandWave( n,how_many )
%Generates random wavefunction
%n=dim of wave vector
%how_many is the number of wave vectors
if nargin==1
    how_many=1;
end

for i=1:how_many
    p=rand(n,1);
    a=sum(p);
    p=p/a;
    c=sqrt(p);
    phi(:,i)=c.*exp(1i*2*pi*rand(n,1));
end
end

