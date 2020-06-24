function [ other_gates sizle anc_sizle names names2 names_sum ] = MWG_Gauss_Jordan_Decomposition_Prepare( H,t,l,rec_order,other_gates )
%Prepares all necessary gates for the Trotter expansion to represent the
%Hamiltionian H's unitary evolution steps precisely
%
% rec_order defines the trotter sheme
    % 1-> quintupling
    % 0-> tripling
    % [] -> symmetric trotter

if length(rec_order)>0
    [seq,d_seq]=Higher_Trotter_Time_Steps(rec_order); %Which time steps have to be created?
    [d_seq2,~,position]=uniquetol(d_seq,10^-14);       % Time steps ,~,position     
    names2={};
    names={};
    names_sum={};
    for nam=1:length(d_seq2)
        %Full steps
        name2=['U_f_t_' num2str(d_seq2(nam)*t/l)];
        circ_name2=['U_{f}(t=' num2str(d_seq2(nam)*t/l) ')'];
        U_2=expm(-1i*H*d_seq2(nam)*t/l);
        gj_dec_gate=Fast_Gauss_Jordan_Decomposition(U_2,0,name2);
        other_gates=Comp_Gate_Merger(other_gates,gj_dec_gate); 
        names2{end+1}=name2;
        
        %Half steps
        name=['U_h_t_' num2str(d_seq2(nam)*t/l/2)];
        circ_name2=['U_{h}(t=' num2str(d_seq2(nam)*t/l/2) ')'];
        U_2=expm(-1i*H*d_seq2(nam)*t/l/2);
        gj_dec_gate=Fast_Gauss_Jordan_Decomposition(U_2,0,name);
        other_gates=Comp_Gate_Merger(other_gates,gj_dec_gate); 
        names{end+1}=name;
    end
    %Sum steps
    d_seq_sum=[d_seq(2:end) d_seq(1)]+d_seq(1:end);
    [d_seq_sum2,~,position_sum]=uniquetol(d_seq_sum,10^-14);
    for nam=1:length(d_seq_sum2)
        name_s=['U_sum_t_' num2str(d_seq_sum2(nam)*t/l/2)];
        circ_name2=['U_{s}(t=' num2str(d_seq_sum2(nam)*t/l/2) ')'];
        U_2=expm(-1i*H*d_seq_sum2(nam)*t/l/2);
        gj_dec_gate=Fast_Gauss_Jordan_Decomposition(U_2,0,name_s);
        other_gates=Comp_Gate_Merger(other_gates,gj_dec_gate); 
        names_sum{end+1}=name_s;
    end
    names2={names2{position}};
    names={names{position}};
    names_sum={names_sum{position_sum}};

    
else %Symmetric Trotter
        
    name2=['U_t_' num2str(t/l)];
    circ_name2=['U(t=' num2str(t/l) ')'];
    U_2=expm(-1i*H*t/l);
    gj_dec_gate2=Fast_Gauss_Jordan_Decomposition(U_2,0,name2);

    name=['U_t_' num2str(t/l/2)];
    circ_name2=['U(t=' num2str(t/l/2) ')'];
    U=expm(-1i*H*t/l/2);              
    gj_dec_gate=Fast_Gauss_Jordan_Decomposition(U,0,name);
    
    other_gates=Comp_Gate_Merger(other_gates,gj_dec_gate2,gj_dec_gate);
    names={name,name2};
end

sizle=gj_dec_gate.size;
anc_sizle=gj_dec_gate.anc_size;
    
end

