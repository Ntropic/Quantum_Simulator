function [ indexes ] = Gray_Indexes_Const_N( n,N )
%GRAY_INDEXES_CONST_N creates a list of indexes in the subspace of the gray
%coding for two modes (or N) and constant total occupation number N=n1+n2
if nargin==1
    indexes=Multi_Gray_Indexes_Const_N( 2,n );
else
    indexes=Multi_Gray_Indexes_Const_N( N,n );
end
end

