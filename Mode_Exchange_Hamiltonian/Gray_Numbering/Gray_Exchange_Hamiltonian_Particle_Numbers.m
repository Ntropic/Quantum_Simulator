function [H,s]= Gray_Exchange_Hamiltonian_Particle_Numbers(n,m)
%PARTICLE_MODE_CHANGE Generates the unitary operator for the exchange of a 
%particle between two modes. 
%   - n particles per mode for storing particle number (via gray code)
%   - 2n qubits in total for interaction
%   - only m particle states
s=ceil(log2(n+1));

H=sparse(2^(2*s),2^(2*s));

%Determine gray code of nums
tree=bintree(s,1);
gray=all_num_algorithm(tree);

if m>n
    error('m must be smaller or equal to n');
end

for num1=1:m+1
    for num2=m+2-num1
        gray_num=gray([num1,num2]);
        graybin=dec2bin(gray_num-1,s)';
        graybin=graybin(:)';
        index=bin2dec(graybin)+1;
        if num1>1 
            gray_num2=gray([num1-1,num2+1]);
            graybin2=dec2bin(gray_num2-1,s)';
            graybin2=graybin2(:)';
            index2=bin2dec(graybin2)+1;
            H(index,index2)=sqrt(num2)*sqrt(num1-1);
        end
        if num2>1 
            gray_num3=gray([num1+1,num2-1]);
            graybin3=dec2bin(gray_num3-1,s)';
            graybin3=graybin3(:)';
            index3=bin2dec(graybin3)+1;
            H(index,index3)=sqrt(num2-1)*sqrt(num1);
        end
    end
end
end