function [ perm ] = Rec_BinTree( N,m )
%Like Recursive_BinTree (which it utilizes) but returns decimals not
%binaries

p=Recursive_BinTree(N,N,m);
perm=bin2dec(num2str(p))+1;
end

