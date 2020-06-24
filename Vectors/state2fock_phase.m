function [ list,s ] = state2fock_phase( n,N,phi,eps,print,name )
%STATE2FOCK creates a list of fockstates in a state

%The hilbert space therefore has a dimension of 2^(n*N)
%Every mode m is represented by s quibits as |m>=|s>|s-1>...|2>|1>_m
%And every vector in the Hilbert space consists of N modes as
%|vec>=|N>|N-1>...|m>...|2>|1>

if nargin==3
    eps=10^(-6);
    print='n';
end
if nargin==4
    print='n';
end
list='';
if length(name)>0
    list=['|' name '>=\n'];
end
s=ceil(log(n+1)/log(2));
dim=2^(s*N);


for index=1:dim
    fock=[];
    i=index-1;
    bin=dec2bin(i,s*N);
    if abs(phi(index))^2>eps
        real_phi=sprintf('%0.4f',full(real(phi(index))));
        imag_phi=sprintf('%0.4f',full(imag(phi(index))));
        if real_phi(1)~='-'
            real_phi=[' ' real_phi];
        end
        if imag_phi(1)~='-'
            imag_phi=[' + ' imag_phi];
        else
            imag_phi=[' - ' imag_phi(2:end)];
        end
        for k=1:N
            bin_comp=bin2dec(bin(s*(k-1)+1:s*k));
            fock=[fock,bin_comp];
        end
        list=[list '\t' real_phi imag_phi '\t|'];
        for k=1:N
            list=[list num2str(fock(k)) ','];
        end
        list=[list(1:end-1) '>\n'];   
    end
end
list=[list '\n'];
if print=='y'
    fprintf(list)
end
end

