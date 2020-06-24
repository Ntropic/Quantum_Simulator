function [ U ] = H2U( H,t )
%H2U creates a unitary operator from a Hamilton operator
U=expm(-1i*H*t);
end

