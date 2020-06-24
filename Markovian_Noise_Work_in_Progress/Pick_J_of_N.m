function [ list ] = Pick_J_of_N( j,n )
%PICK_J_OF_N creates a list with all possibilities for picking j integers 
% among a list of the first n integers

list=[];
%Use recursion
if j==1
    list=1:n;
    list=list';
else
    for i=1:n-j+1
        short_list=Pick_J_of_N(j-1,n-i);
        list=[list;[repmat(i,size(short_list,1),1),short_list+i]];
    end
end

end

