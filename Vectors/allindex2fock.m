function [ list,s ] = allindex2fock( n,N,print )
%ALLINDEX2FOCK creates a list of indices and corresponding fock states

%The hilbert space therefore has a dimension of 2^(n*N)
%Every mode m is represented by s quibits as |m>=|s>|s-1>...|2>|1>_m
%And every vector in the Hilbert space consists of N modes as
%|vec>=|N>|N-1>...|m>...|2>|1>
if nargin==2
    print='n';
end
list='';
s=ceil(log(n+1)/log(2));
dim=2^(s*N);
index_length=['%0' num2str(length(num2str(dim))) 'd'];

for index=1:dim
    fock=[];
    i=index-1;
    bin=dec2bin(i,s*N);
    for k=1:N
        bin_comp=bin2dec(bin(s*(k-1)+1:s*k));
        fock=[fock,bin_comp];
    end
    list=[list sprintf(index_length,index) '\t'];
    for k=1:N
        list=[list '|' num2str(fock(k)) '>'];
    end
    list=[list '\n'];
end
if print=='y'
    fprintf(list)
end
end

