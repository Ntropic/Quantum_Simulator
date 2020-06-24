function [ perm ] = Rec_BinTree( N,m )
%Like Recursive_BinTree (which it utilizes) but returns decimals not
%binaries
m=m-1;
perm=Recursive_BinTree_2(N,m,N,m);
perm=perm+1;
end

