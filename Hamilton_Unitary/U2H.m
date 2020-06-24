function [ H ] = U2H( U )
%H2U creates a unitary operator from a Hamilton operator
if isa(U,'sym')
    diff(U,sym('t'))*conj(U)*i;
else
    error('U needs to be symbolic.')
end
end

