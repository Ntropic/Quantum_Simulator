function [ init main connector outit gate_names] = Split_Trot_Pauli_Ladder_Gates( n,name,circ_string )
%SPLIT_TROT_PAULI_LADDER_GATES (withlinear Trotter) Creates the gates for 
%Pauli ladder encoding and 2 mode interactions for n photons but splits the 
%gate into a init gate the main gate and an outit,
if nargin==3
    circ_string_in= [circ_string '^{(1),(i)}'];
    circ_string_out=[circ_string '^{(1),(f)}'];
    circ_string_con=[circ_string '^{(1),(c)}'];
    circ_string=    [circ_string '^{(1),(m)}'];
    name_in=[name '_in'];
    name_out=[name '_out'];
    name_con=[name '_con'];
    name=[name '_cen'];
elseif nargin==2
    circ_string_in= ['L_{n={' num2str(n) '}}^{(1),(i)}'];
    circ_string_out=['L_{n={' num2str(n) '}}^{(1),(f)}'];
    circ_string_con=['L_{n={' num2str(n) '}}^{(1),(c)}'];
    circ_string=    ['L_{n={' num2str(n) '}}^{(1),(m)}'];
    name_in=[name '_in'];
    name_out=[name '_out'];
    name_con=[name '_con'];
    name=[name '_cen'];
elseif nargin==1
    circ_string_in= ['L_{n={' num2str(n) '}}^{(1),(i)}'];
    circ_string_out=['L_{n={' num2str(n) '}}^{(1),(f)}'];
    circ_string_con=['L_{n={' num2str(n) '}}^{(1),(c)}'];
    circ_string=    ['L_{n={' num2str(n) '}}^{(1),(m)}'];
    name_in=['Ladder_n_' num2str(n) '_step_in'];
    name_out=['Ladder_n_' num2str(n) '_step_out'];
    name_con=['Ladder_n_' num2str(n) '_step_con'];
    name=['Ladder_n_' num2str(n) '_step_cen'];
end
gate_names={name_in,name_out,name,name_con};

s=(n+1)*2; %n+1 qubits for n photons *2 modes

init=Create_Empty_Comp_Gate( name_in,s,0);
init=Generate_Gate_Circuit(init,circ_string_in,[1:s]);
outit=Create_Empty_Comp_Gate( name_out,s,0 );
outit=Generate_Gate_Circuit(outit,circ_string_out,[1:s]);
connector=Create_Empty_Comp_Gate( name_con,s,0 );
connector=Generate_Gate_Circuit(connector,circ_string_con,[1:s]);
main=Create_Empty_Comp_Gate( name,s,0 );
main=Generate_Gate_Circuit(main,circ_string,[1:s]);

phi=sym('phi');

%Difference (d) based model for the ladder
block_nums=[];
prefactor=[];
subsi=[];
% for d=0:n
%     i_max=floor((n+d+1)/2);
%     for i=d+1:i_max
%         j=i-d;
%         block_nums=[block_nums;i j];
%         prefactor=[prefactor;2*sqrt(i*j)];  %Strength of interaction 
%         subsi=[subsi;i*j];
%     end
%     if d~=0
%         for i=d+1:i_max
%             j=i-d;
%             block_nums=[block_nums;j i];
%             prefactor=[prefactor;2*sqrt(i*j)];  %Strength of interaction 
%             subsi=[subsi;i*j];
%         end
%     end
% end
%Even uneven ordering model
for n_now=1:n
    for i=1:2:n_now
        j=n_now+1-i;
        block_nums=[block_nums;i j];
        prefactor=[prefactor;2*sqrt(i*j)];  %Strength of interaction 
        subsi=[subsi;i*j];
    end
end
for n_now=1:n
    for i=2:2:n_now
        j=n_now+1-i;
        block_nums=[block_nums;i j];
        prefactor=[prefactor;2*sqrt(i*j)];  %Strength of interaction 
        subsi=[subsi;i*j];
    end
end

for i=1:size(block_nums,1)
    if prefactor(i)~=1
        subs={[{phi,phi*prefactor(i)}]};
    else
        subs={[]};
    end
    i2=block_nums(i,1);
    j2=block_nums(i,2);
    main=Add_Gate(main,{'L_{+}'},{[i2 i2+1 n+1+j2 n+2+j2]},subs);
    if subsi(i)~=1
        if abs(sqrt(subsi(i))-round(sqrt(subsi(i))))<10^-14
            main.circuit_subs{i}=[{'\phi',[num2str(round(sqrt(subsi(i)))) '\phi' ]}];
        else
            main.circuit_subs{i}=[{'\phi',['\sqrt{' num2str(subsi(i)) '}\phi']}];
        end
    else
        main.circuit_subs{i}=[];
    end
end

end