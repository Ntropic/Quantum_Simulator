function [ fock coeff ] = vec2fock( n,N,vec )
%VEC2FOCK transforms vector states into Fock states
%fock is a vector with the fock=[n_N,...,n_2,n_1] 
%n=2^s is the maximum amount of particles that can be represented by the Hilbert space  
%N is the number of modes represented by the Hilbert space

%The hilbert space therefore has a dimension of 2^(n*N)
%Every mode m is represented by s quibits as |m>=|s>|s-1>...|2>|1>_m
%And every vector in the Hilbert space consists of N modes as
%|vec>=|N>|N-1>...|m>...|2>|1>

s=ceil(log(n+1)/log(2));
fock=[];
coeff=[];

for index=1:length(vec)
    if vec(index)>0
        %Add fock component
        i=index-1;
        bin=dec2bin(i,s*N);
        fock_new=[];
        for k=1:N
            bin_comp=bin2dec(bin(s*k-(s-1):s*k));
            fock_new=[fock_new,bin_comp];
        end
        fock=[fock;fock_new];
        coeff=[coeff,vec(index)];
    end
end
end

