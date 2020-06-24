function [ indexes ] = Multi_Gray_Indexes( m,n )
%MULTI_GRAY_INDEXES creates a list of indexes in the subspace of the gray
%coding for m modes
s=ceil(log2(n+1));
tree=bintree(s,1);
gray=all_num_algorithm(tree);

gray_red=gray(1:n+1)-1;

s_m=2^(s*(m-1));
indexes=[];
if m>2
    for i=1:n+1
        indexes=[indexes,s_m*gray_red(i)+Rec_Multi_Gray_Indexes( m-1,n-i+1,s)];
    end
else
    for i=1:n+1
        for j=1:n+2-i
            indexes=[indexes,s_m*gray_red(i)+gray_red(j)];
        end
    end
end
indexes=indexes+1;
end

