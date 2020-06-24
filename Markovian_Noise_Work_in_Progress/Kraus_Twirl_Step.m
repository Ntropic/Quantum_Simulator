function [ rho ] = Kraus_Twirl_Step( rho,n,K,U )
%KRAUS_TWIRL_STEP creates a noise step with single qubit Kraus matrices K

for j=1:n
    rho2=rho;
    rho=rho*0;
    for i=1:size(K,2)
        rho=rho+K{j,i}*U*K{j,i}*rho2*K{j,i}';
    end
end

