function [ indexes,numbers ] = Rec_Multi_Gray_Indexes_Const_N( m,n,s )
%REC_MULTI_GRAY_INDEXES creates a list of indexes in the subspace of the 
%gray coding for m modes it is used within Multi_Gray_Indexes

tree=bintree(s,1);
gray=all_num_algorithm(tree);

gray_red=gray(1:n+1)-1;

s_m=2^(s*(m-1));
indexes=[];
numbers=[];
if m>2
    for i=1:n+1
        [recs,nums]=Rec_Multi_Gray_Indexes_Const_N( m-1,n-i+1,s);
        indexes=[indexes,s_m*gray_red(i)+recs];
        numbers=[numbers,i-1+nums];
    end
else
    for i=1:n+1
        j=n+2-i;
        indexes=[indexes,s_m*gray_red(i)+gray_red(j)];
        numbers=[numbers,i+j-2];
    end
end
end

