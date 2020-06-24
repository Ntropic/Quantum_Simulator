function [H]= Exchange_Hamiltonian_Particles(n,s)
%PARTICLE_MODE_CHANGE Generates the unitary operator for the exchange of a 
%particle between two modes. 
%   - n particles per mode for storing particle number (via gray code)
%   - 2n qubits in total for interaction
if nargin==1
    s=ceil(log2(n+1));
end

H=zeros(2^(2*s),2^(2*s));

for num1=1:n+1
    for num2=1:n+2-num1
        bin_num=[num1,num2];
        bin=dec2bin(bin_num-1,s)';
        bin=bin(:)';
        index=bin2dec(bin)+1;
        if num1>1 && num2<2^n
            bin_num2=[num1-1,num2+1];
            bin2=dec2bin(bin_num2-1,s)';
            bin2=bin2(:)';
            index2=bin2dec(bin2)+1;
            H(index,index2)=sqrt(num2)*sqrt(num1-1);
        end
        if num2>1 && num1<2^n
            bin_num3=[num1+1,num2-1];
            bin3=dec2bin(bin_num3-1,s)';
            bin3=bin3(:)';
            index3=bin2dec(bin3)+1;
            H(index,index3)=sqrt(num2-1)*sqrt(num1);
        end
    end
end
end