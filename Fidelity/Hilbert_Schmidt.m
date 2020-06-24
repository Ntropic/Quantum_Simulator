function [ HS HS_ind ] = Hilbert_Schmidt( A_ex,A_ap,indexes )
%HILBERT_SCHMIDT calculates the hilbert-schmidt distance also known as
%trace norm -> it's a generalization of the Kolmogorov distance to quantum
%probability distributions
if size(A_ex,1)~=size(A_ex,2)
    error('Matrices must be square (A_ex)')
end
if size(A_ap,1)~=size(A_ap,2)
    error('Matrices must be square (A_ap)')
end
if size(A_ex,1)~=size(A_ap,1)
    error('Matrices must be of equal size')
end

A=A_ex-A_ap;
HS=1/2*trace(sqrt(A'*A));

if nargin==3
    B=A_ap(indexes,indexes);
    nB=trace(sqrt(B'*B));
    A_ap=A_ap/nB;
    
    B=A_ex(indexes,indexes);
    nB=trace(sqrt(B'*B));
    A_ex=A_ex/nB;
end
A=A_ex-A_ap;
if nargin==3
	A=A(indexes,indexes);
end
HS_ind=1/2*trace(sqrt(A'*A));

end

