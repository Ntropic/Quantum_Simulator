function [ vert_string,anc_ket_vert,ket_vert,fock_comp_len ] = vertfock_str( index,mode )
%VERTFOCK_STR creates a vertical fock string (Use
%strrep(vert_string,'\n','') in conjunction, to get rid of the '\n') 

ket_vert_str='';

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
    
ket_vert=cell(length(index),1);
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
    
    anc_ket_vert_str='';
    anc_ket_vert_str(1:fock_comp_len)='_';
    anc_ket_vert_str=[anc_ket_vert_str '\n'];
    
    if length(strfind(mode{1},'fock_bin'))==0
        for j=1:N_anc
            anc_ket_vert_str=[anc_ket_vert_str num2str(fock(j),fock_comp_str) '\n'];
        end
    elseif length(strfind(mode{1},'fock_bin'))==1
        for j=1:N_anc
            anc_ket_vert_str=[anc_ket_vert_str dec2bin(fock(j),s) '\n'];
        end
    end
    
    if fock_comp_len==1
        anc_ket_vert_str=[anc_ket_vert_str char(709)];
    elseif mod(fock_comp_len,2)==0
        for l=1:fock_comp_len/2
            for k=1:l-1
                anc_ket_vert_str=[anc_ket_vert_str ' '];
            end
            anc_ket_vert_str=[anc_ket_vert_str '\\'];
            for k=1:(fock_comp_len)-(l-1)*2-2
                anc_ket_vert_str=[anc_ket_vert_str ' '];
            end
            anc_ket_vert_str=[anc_ket_vert_str char(47)];
            for k=1:l-1
                anc_ket_vert_str=[anc_ket_vert_str ' '];
            end
            if l~=fock_comp_len/2
                anc_ket_vert_str=[anc_ket_vert_str '\n'];
            end
        end
    else %Uneven
        for l=1:(fock_comp_len-1)/2
            for k=1:l-1
                anc_ket_vert_str=[anc_ket_vert_str ' '];
            end
            anc_ket_vert_str=[anc_ket_vert_str '\\'];
            for k=1:(fock_comp_len)-(l-1)*2-2
                anc_ket_vert_str=[anc_ket_vert_str ' '];
            end
            anc_ket_vert_str=[anc_ket_vert_str char(47)];
            for k=1:l-1
                anc_ket_vert_str=[anc_ket_vert_str ' '];
            end
            anc_ket_vert_str=[anc_ket_vert_str '\n'];
        end
        for k=1:((fock_comp_len)-1)/2
            anc_ket_vert_str=[anc_ket_vert_str ' '];
        end
        ket_vert_str=[ket_vert_str char(709)];
        for k=1:((fock_comp_len)-1)/2
            anc_ket_vert_str=[anc_ket_vert_str ' '];
        end
    end
    
    ket_vert_str='';
    ket_vert_str(1:fock_comp_len)='_';
    ket_vert_str=[ket_vert_str '\n'];
    if length(strfind(mode{1},'fock_bin'))==0
        for j=N_anc+1:length(fock)
            ket_vert_str=[ket_vert_str num2str(fock(j),fock_comp_str) '\n'];
        end
    elseif length(strfind(mode{1},'fock_bin'))==1
        for j=N_anc+1:length(fock)
            ket_vert_str=[ket_vert_str dec2bin(fock(j),s) '\n'];
        end
    end
        
    if fock_comp_len==1
        ket_vert_str=[ket_vert_str char(709)];
    elseif mod(fock_comp_len,2)==0
        for l=1:fock_comp_len/2
            for k=1:l-1
                ket_vert_str=[ket_vert_str ' '];
            end
            ket_vert_str=[ket_vert_str '\\'];
            for k=1:(fock_comp_len)-(l-1)*2-2
                ket_vert_str=[ket_vert_str ' '];
            end
            ket_vert_str=[ket_vert_str char(47)];
            for k=1:l-1
                ket_vert_str=[ket_vert_str ' '];
            end
            if l~=fock_comp_len/2
                ket_vert_str=[ket_vert_str '\n'];
            end
        end
    else %Uneven
        for l=1:(fock_comp_len-1)/2
            for k=1:l-1
                ket_vert_str=[ket_vert_str ' '];
            end
            ket_vert_str=[ket_vert_str '\\'];
            for k=1:(fock_comp_len)-(l-1)*2-2
                ket_vert_str=[ket_vert_str ' '];
            end
            ket_vert_str=[ket_vert_str char(47)];
            for k=1:l-1
                ket_vert_str=[ket_vert_str ' '];
            end
            ket_vert_str=[ket_vert_str '\n'];
        end
        for k=1:((fock_comp_len)-1)/2
            ket_vert_str=[ket_vert_str ' '];
        end
        ket_vert_str=[ket_vert_str char(709)];
        for k=1:((fock_comp_len)-1)/2
            ket_vert_str=[ket_vert_str ' '];
        end
    end
    anc_ket_vert{i}=anc_ket_vert_str;
    ket_vert{i}=ket_vert_str;
    if N_anc>0
        vert_string{i}=[anc_ket_vert_str '\n' ket_vert_str];
    else
        vert_string{i}=ket_vert_str;
    end
end
end

