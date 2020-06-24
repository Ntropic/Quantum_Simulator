function [ PU_gates, Cn ] = CreateCompleteSet_Cn_U2( n, elem_gates, comp_gates, save_load_anc, printing)
%CREATECOMPLETESET_CN_U2 generates a complete set of n(2^n-1) conditionally
%controlled gates CN_U2 that are triggered by a bit sequence of n-1 bits
%save_load determines if gates are to be saved (1) or in case they exist to be
%loaded (2) and if anc_matrix is to be determined(3)
ver=1.2;    %Version 1.1 adds support for circuit diagrams
if nargin==3
    save_load_anc=[0 0 0];
    printing=0;
elseif nargin==4
    printing=0;
end

loaded=0;
if save_load_anc(2)==1 %Try to load data
    path=which('quantum_gates');
    path_elem=strsplit(path,'\');
    path=strjoin(path_elem((1:length(path_elem)-1)),'\');
    path=[path,'\C',num2str(n),'_U2.mat'];
    if exist(path)==2
        S=load(path);
        if S.ver==ver
            loaded=1;
        end
        PU_gates=S.PU_gates;
        Cn=S.Cn;
    end
end

if loaded==0
    seq_nums=0:2^n-1;
    seq=dec2bin(seq_nums,n)-'0';
    if save_load_anc(3)==1
        seq_anc=dec2bin(1:2^(n-1)-1)-'0';
    end
    
    cn_u2_names=cell(size(seq,1));
    Ph_names=cell(1,2^n-1);
    Cn=struct('names',Ph_names);
    
    %First create the basis matrices via plugin
    U=Gate_by_Name('U',elem_gates,comp_gates);
    diag_mat=sym(diag(repmat(1,2^(n+1),1)));
    
    if save_load_anc(3)==1
        diag_anc_mat=sym(diag(repmat(1,2^(n*2),1)));
    end
    for l=1:size(seq,1)
        cnu2=Cn_U2(elem_gates,comp_gates,seq(l,n:-1:1),[]);
        for fn=fieldnames(cnu2)'
            Cn(l).(fn{1})=cnu2.(fn{1});
        end
        matrix=diag_mat;

        mat_list=2*l-1:2*l;
        matrix(mat_list,mat_list)=U.matrix;
        if save_load_anc(3)==1
            anc_matrix=diag_anc_mat;
            anc_matrix(mat_list,mat_list)=U.matrix;
        end
        %Find other possible combinations
        if save_load_anc(3)==1
            if n>=2
                for i=1:length(seq_anc)
                    for j=1:length(seq)
                        curr_seq_anc=seq_anc(i,(n-1):-1:1);
                        curr_seq=seq(j,n:-1:1);
                        comp_seq=seq(l,n:-1:1);
                        if curr_seq(1:2)==comp_seq(1:2)
                            seq_anc(1)=mod(seq_anc(1)+1,2);
                        end
                        for k=3:n
                            if curr_seq(k)==comp_seq(k) && seq_anc(k-2)==1
                                seq_anc(k-1)=mod(seq_anc(k-1)+1,2);
                            end
                        end
                        if seq_anc(n-1)==1
                            curr_index=bin2dec(num2str([seq_anc(i,(n-1):-1:1) seq(j,n:-1:1) 0]));
                            anc_matrix(curr_index:curr_index+1,curr_index:curr_index+1)=U.matrix;
                        end
                    end
                    
                end
            end
            Cn(l).anc_matrix=anc_matrix;
        end
        Cn(l).matrix=matrix;
    end
    k=1:n+1;

    cn_u2_names=cell(length(seq),n+1);
    PU_gates=struct('names',cn_u2_names);
    timer=tic();
    for i=1:size(seq,1)
        if printing==1
            if i>1
                elapsed_time=toc(timer);
                fprintf('-> Constructing: i=%1d / %1d - rem: %1.2d\n',i,size(seq,1),elapsed_time/(i-1)*(size(seq,1)-i));
            else
                fprintf('-> Constructing: i=%1d / %1d \n',i,size(seq,1));
            end
        end
        k3=dec2bin(seq_nums(i),n);
        for j=1:n+1
            k=[j 1:(j-1) (j+1):(n+1)];
            k2=[k3(1:(n+1-j)) 'x' k3((n+2-j):n)];

            matrix=Embed_Gate(Cn(i).matrix,k,n+1);
            if save_load_anc(3)==1
                anc_matrix=Embed_Gate(Cn(i).anc_matrix,[k (1:Cn(i).anc_size)+length(k)],n+1+Cn(i).anc_size);
            end
            for fn=fieldnames(Cn(i))'
                PU_gates(i,j).(fn{1})=Cn(i).(fn{1});
            end
            PU_gates(i,j).names=['PU_' num2str(n+1) '_' k2];
            PU_gates(i,j).matrix=matrix;
            if save_load_anc(3)==1
                PU_gates(i,j).anc_matrix=anc_matrix;
            end
            for l=1:length(PU_gates(i,j).steps.index)
                o_list=find(PU_gates(i,j).steps.index{l}==1);
                j_list=find(PU_gates(i,j).steps.index{l}==j);
                PU_gates(i,j).steps.index{l}(o_list)=j;
                PU_gates(i,j).steps.index{l}(j_list)=1;
            end
            
            k2_x=strfind(k2,'x');
            k2_y=strfind(k2,'y');
            k2_y_list=dec2bin(0:2^length(k2_y)-1,length(k2_y));
            PU_gates(i,j).fock_index=zeros(2,size(k2_y_list,1));
            for q=1:size(k2_y_list,1)
                k2_sub0=k2;
                k2_sub0(k2_y)=k2_y_list(q,:);
                k2_sub1=k2_sub0;
                k2_sub0(k2_x)='0';
                k2_sub1(k2_x)='1';
                PU_gates(i,j).fock_index(1,q)=bin2dec(k2_sub0)+1;
                PU_gates(i,j).fock_index(2,q)=bin2dec(k2_sub1)+1;
            end
           	PU_gates(i,j).circuit={[1:length(k);ones(1,length(k))]};
            for q=length(k2):-1:1
                if k2(q)=='1'
                    PU_gates(i,j).circuit{2-q+length(k2)}=['\ctrl{_' num2str(length(k2)-k2_x+1) '_}'];
                elseif k2(q)=='x'
                    PU_gates(i,j).circuit{2-q+length(k2)}='\gate{U}';
                elseif k2(q)=='0'
                    PU_gates(i,j).circuit{2-q+length(k2)}=['\ctrlo{_' num2str(length(k2)-k2_x+1) '_}'];
                end     
            end
        end
    end
    if save_load_anc(1)==1 %Save results
        path=which('quantum_gates');
        path_elem=strsplit(path,'\');
        path=strjoin(path_elem((1:length(path_elem)-1)),'\');
        path=[path,'\C',num2str(n),'_U2.mat'];
        
        save(path,'ver','PU_gates','Cn');
    end
end

end

