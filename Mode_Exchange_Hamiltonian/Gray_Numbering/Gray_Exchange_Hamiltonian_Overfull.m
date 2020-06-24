function [H]= Gray_Exchange_Hamiltonian_Overfull(n)
%PARTICLE_MODE_CHANGE Generates the unitary operator for the exchange of a 
%particle between two modes. 
%   - n quibts per mode for storing particle number (via gray code)
%   - 2n qubits in total for interaction

H=sparse(2^(2*n),2^(2*n));

%Determine gray code of nums
tree=bintree(n,1);
gray=all_num_algorithm(tree);


for num1=1:2^n
    for num2=1:2^n
        gray_num=gray([num1,num2]);
        graybin=dec2bin(gray_num-1,n)';
        graybin=graybin(:)';
        index=bin2dec(graybin)+1;
        if num1>1 && num2<2^n
            gray_num2=gray([num1-1,num2+1]);
            graybin2=dec2bin(gray_num2-1,n)';
            graybin2=graybin2(:)';
            index2=bin2dec(graybin2)+1;
            H(index,index2)=sqrt(num2)*sqrt(num1-1);
        end
        if num2>1 && num1<2^n
            gray_num3=gray([num1+1,num2-1]);
            graybin3=dec2bin(gray_num3-1,n)';
            graybin3=graybin3(:)';
            index3=bin2dec(graybin3)+1;
            H(index,index3)=sqrt(num2-1)*sqrt(num1);
        end
    end
end
end