function [ fock ] = index2gray_fock_multi( n,N,index )
%INDEX2GRAY_FOCK transforms index's into Gray_Fock states
%fock is a vector with the fock=[n_N,...,n_2,n_1] 
%n (number of bits per mode) -> s (with 2^s-1=n) is the maximum amount of particles that can be represented by the Hilbert space  
%N is the number of modes represented by the Hilbert space

%The hilbert space therefore has a dimension of 2^(n*N)
%Every mode m is represented by s quibits as |m>=|s>|s-1>...|2>|1>_m
%And every vector in the Hilbert space consists of N modes as
%|vec>=|N>|N-1>...|m>...|2>|1>

perm=Gray_NumberPath(n)-1;

s=ceil(log(n+1)/log(2));
fock=cell(length(index),1);
if length(N)==2
    N_anc=N(2);
else
    N_anc=0;
end
N=N(1);
q_num=s*N+N_anc;
    
index=index-1;
for i=1:length(index)
    fock2=dec2bin(index(i),q_num)-'0';
    fock{i}=zeros(N+N_anc,1);
    fock{i}(1:N_anc)=fock2(1:N_anc);
    for j=1:N
        fock{i}(N_anc+j)=find(perm==bin2dec(num2str(fock2(N_anc+1+(j-1)*s:N_anc+j*s))))-1;
    end
end
end

