function [ gate ] = MWG_Trotter( gate_names,connections,step_num,name,circ_string,sizle,anc_sizle,t )
%MWG_TROTTER_INOUT does the Lie-Suzuki-Trotter decomposition for
%gates_names={init,outit,center,connect gates}
%Create sequence (seq) and step width d(seq) für multiple waveguides
%weights defines the interaction strengths of N modes   
%size(connections)=[2,interactions];

siz=sizle/2;
if siz~=round(siz)
    error('Size needs to be even qubits big.')
end

if length(connections)==1
    WG=connections;
    connections=[(1:WG-1)',(2:WG)'];
else
    WG=max(connections(:));
end
if isa(gate_names,'char')
    g=gate_names;
    gate_names={gate_names};
    for i=2:size(connections,1)
        gate_names={gate_names{:},g};
    end
end

if length(gate_names)==1
    g=gate_names{1};
    gate_names={gate_names};
    for i=2:size(connections,1)
        gate_names={gate_names{:},g};
    end
end

s=siz*WG;
anc_s=(1:anc_sizle)+s;

gate_steps={};
index_steps={};
param_steps={};

for i=1:step_num
    for j=1:size(connections,1)
        in=[(connections(j,1)-1)*siz+(1:siz),(connections(j,2)-1)*siz+(1:siz),anc_s];

        gate_steps={gate_steps{:},gate_names{j}};
        index_steps={index_steps{:},in};
        param_steps={param_steps{:},[{}]};
    end
end




if length(name)==0
    name=['T_{N=' num2str(WG) ',n=' num2str(sizle) '}^{(' num2str(length(step_num)) ')}(' latex(sym(t)) ')'];
end
if length(circ_string)==0
    circ_string=['T_{N=' num2str(WG) ',n=' num2str(sizle) '}^{(' num2str(length(step_num)) ')}(' latex(sym(t)) ')'];
end

gate=Create_Comp_Gate(name,s,anc_sizle,gate_steps,index_steps,param_steps);
gate=Generate_Gate_Circuit(gate,circ_string,1:sizle,[],[]);
end