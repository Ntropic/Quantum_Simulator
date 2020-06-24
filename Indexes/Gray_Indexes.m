function [ indexes ] = Gray_Indexes( n,N )
%GRAY_INDEXES creates a list of indexes in the subspace of the gray
%coding for two modes
if nargin==1
    s=ceil(log2(n+1));
    tree=bintree(s,1);
    gray=all_num_algorithm(tree);

    gray_red=gray(1:n+1)-1;

    s_m=2^s;
    indexes=[];
    for i=1:n+1
        for j=1:n+2-i
            indexes=[indexes,s_m*gray_red(i)+gray_red(j)];
        end
    end
    indexes=indexes+1;
else
    indexes=Multi_Gray_Indexes( N,n );
    
end
end

