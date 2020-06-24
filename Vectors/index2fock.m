function [ fock ] = index2fock( n,N,index )
%INDEX2FOCK transforms index's into Fock states
%fock is a vector with the fock=[n_N,...,n_2,n_1] 
%n (number of bits per mode) -> s (with 2^s-1=n) is the maximum amount of particles that can be represented by the Hilbert space  
%N is the number of modes represented by the Hilbert space

%The hilbert space therefore has a dimension of 2^(n*N)
%Every mode m is represented by s quibits as |m>=|s>|s-1>...|2>|1>_m
%And every vector in the Hilbert space consists of N modes as
%|vec>=|N>|N-1>...|m>...|2>|1>
index2=index;
fock=[];
for j=1:length(index);
    index=index2(j);
    s=ceil(log(n+1)/log(2));
    focke=[];
    if length(N)==2
        N_anc=N(2);
    else
        N_anc=0;
    end
    N=N(1);
    q_num=s*N+N_anc;

    index=index-1;
    fock2=dec2bin(index,q_num)-'0';
    focke=zeros(N+N_anc,1);
    focke(1:N_anc)=fock2(1:N_anc);
    for i=1:N
        focke(N_anc+i)=bin2dec(num2str(fock2(N_anc+1+(i-1)*s:N_anc+i*s)));
    end
    fock=[fock,focke];
end
end

