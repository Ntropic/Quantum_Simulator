function [ new_gate ] = Cn_U2( elem_gates, comp_gates, seq, U, make_mat )
%CN_U creates a C_n-U gate from toffoli gates. It uses n-1 ancilla
%qbits and 2*(n-1) Toffoli gates and one U gate. The sequence in seq
%specifies which state of control bits triggers the U gate. It can be given
%as a fock state or index. U is given by [alpha,beta,delta,theta] or a 
%matrix 
%make_mat=[bin(matrix) bin(anc_matrix)] determines if the (anc_matrix) shall be created on the fly 
if nargin==4
    make_mat=[0 0];
elseif nargin==5
    if length(make_mat)==1
        make_mat(2)=0;
    end
end

n=length(seq);
N=n+1;
anc_N=max([n-1,0]);
if length(seq)==0
    %Prepare Toffoli sequence
    new_gate.size=1;
    new_gate.anc_size=0;
    new_gate.step_num=1;
    new_gate.steps.index={};
    new_gate.steps.gates={};
    new_gate.steps.param={};
    new_gate.steps.gates={new_gate.steps.gates{:},'U'};

    if length(U)==0
        new_gate.names=['C' num2str(n) '_' num2str(seq(end:-1:1),'%1d') '_U'];
        new_gate.steps.index={new_gate.steps.index{:},[1]};
        new_gate.steps.param={new_gate.steps.param{:},[]};
    elseif length(U)==4
        new_gate.names=['C' num2str(n) '_' num2str(seq(end:-1:1),'%1d') '_U_a' num2str(U(1)) '_b_' num2str(U(2)) '_d_' num2str(U(3)) '_t_' num2str(U(1))];
        new_gate.steps.index={new_gate.steps.index{:},[1]};
        new_gate.steps.param={new_gate.steps.param{:},U(:)'};
    elseif size(U)==[2,2]
        %matrix_small=Gate2Matrix( elem_gates,comp_gates, U_gate,2);
        new_gate.steps.index={new_gate.steps.index{:},[1]};
        [alpha,beta,delta,theta]=Unitary2Angles(U);
        new_gate.steps.param={new_gate.steps.param{:},[alpha,beta,delta,theta]};
        new_gate.names=['C' num2str(n) '_' num2str(seq(end:-1:1),'%1d') '_U_a' num2str(U(1)) '_b_' num2str(U(2)) '_d_' num2str(U(3)) '_t_' num2str(U(1))];
    end
    
    new_gate.circuit={[1;1]};
    new_gate.circuit{2}='\gate{U}';
else
    %Prepare Toffoli sequence
    new_gate.size=n+1;
    new_gate.anc_size=n-1;
    new_gate.step_num=(n-1)*2+1;
    new_gate.steps.index={};
    new_gate.steps.gates={};
    new_gate.steps.param={};

    %First Toffoli
    if n>=2
        if seq(1:2)==[1 1]
            new_gate.steps.gates={new_gate.steps.gates{:},'toffoli11'};
            new_gate.steps.index={new_gate.steps.index{:},[n+2,2,3]};
            new_gate.steps.param={new_gate.steps.param{:},[]};
        elseif seq(1:2)==[0 1]
            new_gate.steps.gates={new_gate.steps.gates{:},'toffoli01'};
            new_gate.steps.index={new_gate.steps.index{:},[n+2,2,3]};
            new_gate.steps.param={new_gate.steps.param{:},[]};
        elseif seq(1:2)==[1 0]
            new_gate.steps.gates={new_gate.steps.gates{:},'toffoli01'};
            new_gate.steps.index={new_gate.steps.index{:},[n+2,3,2]};
            new_gate.steps.param={new_gate.steps.param{:},[]};
        elseif seq(1:2)==[0 0]
            new_gate.steps.gates={new_gate.steps.gates{:},'toffoli00'};
            new_gate.steps.index={new_gate.steps.index{:},[n+2,2,3]};
            new_gate.steps.param={new_gate.steps.param{:},[]};
        else
            error('Sequence has to consist of 0s and 1s.')
        end
    end

    %The remaining toffolis
    for i=4:n+1
        if seq(i-1)==0
            new_gate.steps.gates={new_gate.steps.gates{:},'toffoli01'};
            new_gate.steps.index={new_gate.steps.index{:},[n+i-1,i,n+i-2]};
            new_gate.steps.param={new_gate.steps.param{:},[]};
        elseif seq(i-1)==1
            new_gate.steps.gates={new_gate.steps.gates{:},'toffoli11'};
            new_gate.steps.index={new_gate.steps.index{:},[n+i-1,i,n+i-2]};
            new_gate.steps.param={new_gate.steps.param{:},[]};
        end
    end

    %The unitary Trafo
    %Create small matrix from information in U (U.name and U.param)
    if n==1 && seq==0
        new_gate.steps.gates={new_gate.steps.gates{:},'C0U'};
    else
        new_gate.steps.gates={new_gate.steps.gates{:},'CU'};
    end
    if length(U)==0
        new_gate.names=['C' num2str(n) '_' num2str(seq(end:-1:1),'%1d') '_U'];
        new_gate.steps.index={new_gate.steps.index{:},[1,n*2]};
        new_gate.steps.param={new_gate.steps.param{:},[]};
    elseif length(U)==4
        new_gate.names=['C' num2str(n) '_' num2str(seq(end:-1:1),'%1d') '_U_a' num2str(U(1)) '_b_' num2str(U(2)) '_d_' num2str(U(3)) '_t_' num2str(U(1))];
        new_gate.steps.index={new_gate.steps.index{:},[1,n*2]};
        new_gate.steps.param={new_gate.steps.param{:},U(:)'};
    elseif size(U)==[2,2]
        %matrix_small=Gate2Matrix( elem_gates,comp_gates, U_gate,2);
        new_gate.steps.index={new_gate.steps.index{:},[1,n*2]};
        [alpha,beta,delta,theta]=Unitary2Angles(U);
        new_gate.steps.param={new_gate.steps.param{:},[alpha,beta,delta,theta]};
        new_gate.names=['C' num2str(n) '_' num2str(seq(end:-1:1),'%1d') '_U_a' num2str(U(1)) '_b_' num2str(U(2)) '_d_' num2str(U(3)) '_t_' num2str(U(1))];
    end

    new_gate.steps.gates={new_gate.steps.gates{:},new_gate.steps.gates{(n-1):-1:1}};
    new_gate.steps.index={new_gate.steps.index{:},new_gate.steps.index{(n-1):-1:1}};
    new_gate.steps.param={new_gate.steps.param{:},new_gate.steps.param{(n-1):-1:1}};
    
    new_gate.circuit={[1:(length(seq)+1);ones(1,length(seq)+1)]};
    new_gate.circuit{2}='\gate{U}';
    for q=1:length(seq)
        if seq(q)==1
            new_gate.circuit{2+q}=['\ctrl{_1_}'];
        elseif seq(q)==0
            new_gate.circuit{2+q}=['\ctrlo{_1_}'];
        end     
    end
end

%% --------------------------------------------------------------------
if make_mat(1)==1 || make_mat(2)==1    
    U_sym=0;
    if length(U)==0
        Um=Gate_by_Name('U',elem_gates,comp_gates);
        U=Um.matrix;
        U_sym=1;
    elseif length(U)==4
        U=[U(1) U(2);U(3) U(4)];
    end
end
if make_mat(1)==1
    ind=bin2dec(num2str([seq(n:-1:1) 0]))+1;
    ind=[ind,ind+1];
    if U_sym
        matrix=sym(diag(repmat(1,2^(n+1),1)));
    else
        matrix=zeros(2^(n+1));
    end
    matrix(ind,ind)=U;
    new_gate.matrix=matrix;
              
end
if make_mat(2)==1
    ind=bin2dec(num2str([seq(n:-1:1) 0]))+1;
    ind=[ind,ind+1];
    if U_sym
        matrix=sym(diag(repmat(1,2^(n+1+anc_N),1)));
    else
        matrix=zeros(2^(n+1+anc_N));
    end
    matrix(ind,ind)=U;
    
    %Find other possible combinations
    seqs=dec2bin(0:2^n-1,n)-'0';
    seq_anc=dec2bin(1:2^(n-1)-1)-'0';
    if n>=2
        for i=1:length(seq_anc)
            for j=1:length(seqs)
                curr_seq_anc=seq_anc(i,(n-1):-1:1);
                curr_seq=seqs(j,n:-1:1);
                if curr_seq(1:2)==seq(1:2)
                    seq_anc(1)=mod(seq_anc(1)+1,2);
                end
                for k=3:n
                    if curr_seq(k)==seq(k) && seq_anc(k-2)==1
                        seq_anc(k-1)=mod(seq_anc(k-1)+1,2);
                    end
                end
                if seq_anc(n-1)==1
                    curr_index=bin2dec(num2str([seq_anc(i,(n-1):-1:1) seqs(j,n:-1:1) 0]));
                    anc_matrix(curr_index:curr_index+1,curr_index:curr_index+1)=U;
                end
            end
        end
    end
    new_gate.anc_matrix=anc_matrix;
end
end