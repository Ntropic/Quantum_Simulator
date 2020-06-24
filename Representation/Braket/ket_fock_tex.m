function [ ket_string ket anc_ket ] = ket_fock_str( index,mode )
%KET_FOCK_STR creates a ket fock string 

if nargin==2
    fock_comp_len=0;
end

ket_str='';

%Create Fock-states
n=mode{2};
N=mode{3};
if length(mode)==4
    N_anc=mode{4};
else
    N_anc=0;
end

if any(strfind(mode{1},'fock_bin'))==0
    fock_comp_len=0;
    for i=1:length(index)
        if length(mode)==4
            fock=index2fock(n,[N N_anc],index(i));
        else
            fock=index2fock(n,N,index(i));
        end
        for j=1:length(fock)
            st=num2str(fock(j));
            if length(st)>fock_comp_len
                fock_comp_len=length(st);
            end
        end
    end
else
    fock_comp_len=ceil(log(n+1)/log(2));
    s=fock_comp_len;
end

fock_comp_str=['%0' num2str(fock_comp_len) 'd'];
    
fock_bin=length(strfind(mode{1},'fock_bin'));
fock_gray=length(strfind(mode{1},'fock_gray'));

ket=cell(length(index),1);
anc_ket=cell(length(index),1);
ket_string=cell(length(index),1);

if fock_gray>0
    if length(mode)==4
        focks=index2gray_fock_multi(n,[N N_anc],index);
    else
        focks=index2gray_fock_multi(n,N,index);
    end
end
for i=1:length(index)
    if length(mode)==4
        if fock_gray>0
            fock=focks{i};
        else
            fock=index2fock(n,[N N_anc],index(i));
        end
    else
        if fock_gray>0
            fock=focks{i};
        else
            fock=index2fock(n,N,index(i));
        end
    end 
    
    anc_ket_str='\\ket{';
    if length(strfind(mode{1},'fock_bin'))==0
        for j=1:N_anc
            anc_ket_str=[anc_ket_str num2str(fock(j)) ','];
        end
    elseif fock_bin>0
        for j=1:N_anc
            anc_ket_str=[anc_ket_str dec2bin(fock(j)) ','];
        end
    end
    anc_ket_str=[anc_ket_str(1:end-1) '}'];
    
    ket_str='\\ket{';
    if length(strfind(mode{1},'fock_bin'))==0
        for j=N_anc+1:length(fock)
            ket_str=[ket_str num2str(fock(j),fock_comp_str) ','];
        end
    elseif fock_bin>0
        for j=N_anc+1:length(fock)
            ket_str=[ket_str dec2bin(fock(j),s) ','];
        end
    end
    ket_str=[ket_str(1:end-1) '}'];

    ket{i}=ket_str;
    anc_ket{i}=anc_ket_str;
    if N_anc>0
        ket_string{i}=[anc_ket_str ket_str];
    else
        ket_string{i}=ket_str;
    end
end
end

