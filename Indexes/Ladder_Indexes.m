function [ indexes ] = Ladder_Indexes( n,N )
%LADDER_INDEXES creates a list of indexes in the subspace of the ladder
%coding for two modes
indexes=[];

if nargin==1
    for j=0:n
        for i=0:n-j
            ind=2^i+2^(j+n+1)+1;
            indexes=[indexes,ind];
        end
    end
else
    indexes=Multi_Ladder_Indexes(N,n);
end
end

