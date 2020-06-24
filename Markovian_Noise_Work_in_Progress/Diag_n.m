function [ D ] = Diag_n( n )
%DIAG_N creates diagonal one matrices for up to n qubits
D={};
for i=1:n
    D={D{:},diag(ones(2^i,1))};
end
end

