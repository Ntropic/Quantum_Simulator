function [ vec,index ] = fock2vec( n,N,fock,weights )
%multiFOCK2VEC transforms multiple Fock states into vector states
%fock is a vector with the fock=[n_N,...,n_2,n_1] or a matrix of multiple focks 
%n=2^s-1 is the maximum amount of particles that can be represented by the 
%   Hilbert space  
%N is the number of modes represented by the Hilbert space

%The hilbert space therefore has a dimension of 2^(s*N)
%Every mode m is represented by s quibits as |m>=|s>|s-1>...|2>|1>_m
%And every vector in the Hilbert space consists of N modes as
%|vec>=|N>|N-1>...|m>...|2>|1>

s=ceil(log(n+1)/log(2));
vec=sparse(2^(s*N),1);

if size(fock,2)>N
    warning('Fock vector is longer than the number of modes')
end
for f=1:size(fock,1)
    vecstring='';
    for i=1:N
        vecstring=[vecstring dec2bin(fock(f,i),s)];
    end

    index=bin2dec(vecstring)+1;
    vec(index)=weights(f);
end
end

