function [ ket_string] = vertfock_str( index,mode,mode2 )
%VERTFOCK_STR creates a vertical fock string (Use
%strrep(vert_string,'\n','') in conjunction, to get rid of the '\n') 

ket_vert_str='';

if nargin==2
    mode2='';
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

fock_bin=length(strfind(mode{1},'fock_bin'));
fock_gray=length(strfind(mode{1},'fock_gray'));

fock_comp_str=['%0' num2str(fock_comp_len) 'd'];
    
ket=cell(length(index),1);
anc_ket=cell(length(index),1);
ket_string=cell(length(index),1);
for i=1:length(index)
    if length(mode)==4
        if fock_gray>0
            fock=index2gray_fock(n,[N N_anc],index(i));
        else
            fock=index2fock(n,[N N_anc],index(i));
        end
    else
        if fock_gray>0
            fock=index2gray_fock(n,N,index(i));
        else
            fock=index2fock(n,N,index(i));
        end
    end 

    if length(strfind(mode2,'mathbf'))==0
        anc_ket_str='\\text{\\begin{sideways}$\\ket{';
    else
        anc_ket_str='\\text{\\begin{sideways}$\\ket{\\mathbf{';
    end
    if length(strfind(mode{1},'fock_bin'))==0
        for j=1:N_anc
            anc_ket_str=[anc_ket_str num2str(fock(j)) ','];
        end
    elseif length(strfind(mode{1},'fock_bin'))==1
        for j=1:N_anc
            anc_ket_str=[anc_ket_str dec2bin(fock(j)) ','];
        end
    end
    if N_anc>0
        ket_str='';
    else
        if length(strfind(mode2,'mathbf'))==0
            anc_ket_str=[anc_ket_str(1:end-1) '}$\\end{sideways}}'];
            ket_str='\\text{\\begin{sideways}$\\ket{';
        else
            anc_ket_str=[anc_ket_str(1:end-1) '}}$\\end{sideways}}'];
            ket_str='\\text{\\begin{sideways}$\\ket{\\mathbf{';
        end
    end
    if length(strfind(mode{1},'fock_bin'))==0
        for j=N_anc+1:length(fock)
            ket_str=[ket_str num2str(fock(j),fock_comp_str) ','];
        end
    elseif length(strfind(mode{1},'fock_bin'))==1
        for j=N_anc+1:length(fock)
            ket_str=[ket_str dec2bin(fock(j),s) ','];
        end
    end
    if length(strfind(mode2,'mathbf'))==0
        ket_str=[ket_str(1:end-1) '}$\\end{sideways}}'];
    else
        ket_str=[ket_str(1:end-1) '}}$\\end{sideways}}'];
    end
    ket{i}=ket_str;
    anc_ket{i}=anc_ket_str;
    if N_anc>0
        ket_string{i}=[anc_ket_str ket_str];
    else
        ket_string{i}=ket_str;
    end
end
end

