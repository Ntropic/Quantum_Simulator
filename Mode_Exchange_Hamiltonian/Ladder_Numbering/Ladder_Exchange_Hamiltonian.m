function [H]= Ladder_Exchange_Hamiltonian(n)
%LADDER_MODE_CHANGE Generates the unitary operator for the exchange of a 
%particle between two modes in ladder coding. for n photons

s=2*(n+1);
H=sparse(2^s,2^s);

for num1=0:n
    for num2=0:n-num1
        index=2^num1+2^(n+1+num2)+1;
        if num1>0 
            index2=2^(num1-1)+2^(n+1+num2+1)+1;
            H(index,index2)=sqrt(num2+1)*sqrt(num1);
        end
        if num2>0
            index2=2^(num1+1)+2^(n+1+num2-1)+1;
            H(index,index2)=sqrt(num2)*sqrt(num1+1);
        end
    end
end
end