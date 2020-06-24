function [ rho ] = CZ_Twirl_Step( rho, K, p, index )
%EFFICIENT_CZ_TWIRL applies CZ twirl after U_CZ *rho*U'_CZ

rho2=rho;
n=log2(length(rho));
rho=rho*p(1);
for i=1:length(K)
    K2=Embed_Gate(K{i},index,n);
    rho=rho+p(i+1)*K2*rho*K2';
end
end

