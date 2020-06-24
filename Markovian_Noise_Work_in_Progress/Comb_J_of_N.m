function [ list ] = Comb_J_of_N( j,n )
%COMB_J_OF_N creates a list with all possibilities of j numbers among the n
%first integers
%Base n numbers with length j

list=[];
%Use recursion
if j==1
    list=1:n;
    list=list';
else
    for i=1:n
        short_list=Comb_J_of_N(j-1,n);
        list=[list;[repmat(i,size(short_list,1),1),short_list]];
    end
end

end

