function [ paths ] = Gray_NumberPath( m )
%Gray_NumberPaths creates the gray permutation for all m bit long numbers
%(starting at 1)

paths(:)=Rec_BinTree(m,1); %10 times as fast -> doesn't require huge matrix
%tree=bintree(m,1);
%paths(:)=all_num_algorithm(tree);
end

