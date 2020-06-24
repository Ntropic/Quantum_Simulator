function [ H ] = Exchange_Hamiltonian(n)
%Exchange_Hamiltonian Create n qubit exchange matrix for N modes with
%exchange between modes ij and ij+1


H=zeros(2^(2*n),2^(2*n));

for num1=1:2^n
    for num2=1:2^n-num1+1
        bin_num=[num1,num2];
        bin=dec2bin(bin_num-1,n)';
        bin=bin(:)';
        index=bin2dec(bin)+1;
        if num1>1 && num2<2^n
            bin_num2=[num1-1,num2+1];
            bin2=dec2bin(bin_num2-1,n)';
            bin2=bin2(:)';
            index2=bin2dec(bin2)+1;
            H(index,index2)=sqrt(num2)*sqrt(num1-1);
        end
        if num2>1 && num1<2^n
            bin_num3=[num1+1,num2-1];
            bin3=dec2bin(bin_num3-1,n)';
            bin3=bin3(:)';
            index3=bin2dec(bin3)+1;
            H(index,index3)=sqrt(num2-1)*sqrt(num1);
        end
    end
end
end
