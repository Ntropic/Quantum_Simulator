function [ paths ] = NumberPaths( m )
%NumberPaths creates a permutations of the numbers n to 2^m (for all n=1:2^m-1), where every
%subsequent number differs in only one binary digit from it's neighbors

paths=zeros(2^m,2^m);
for i=1:2^m-1
    tree=bintree(m,i);
    paths(i:2^m,i)=all_num_algorithm(tree);
end
end

