function [ gate ] = MWG_Trotter_InOut( t,connections,weights,step_num,name,circ_string,sizle,anc_sizle,gates_names )
%MWG_TROTTER_INOUT does the Lie-Suzuki-Trotter decomposition for
%gates_names={init,outit,center,connect gates}
%Create sequence (seq) and step width d(seq) für multiple waveguides
%weights defines the interaction strengths of N modes   
%size(connections)=[2,interactions];
t1=sym('t');
phi=sym('phi');

siz=sizle/2;
if siz~=round(siz)
    error('Size needs to be even qubits big.')
end

if length(connections)==1
    WG=connections;
    weights=ones(WG-1,1);
    connections=[(1:WG-1)',(2:WG)'];
else
    WG=max(connections(:));
end

s=siz*WG;
anc_s=(1:anc_sizle)+s;

init=gates_names{1};
outit=gates_names{2};
center=gates_names{3};
connect=gates_names{4};


gate_steps={};
index_steps={};
if length(t)==1 && length(step_num)==1
    
    param_steps={};
    for i=1:step_num
        for j=[1:2:length(weights) 2:2:length(weights)]
            in=[(connections(j,1)-1)*siz+(1:siz),(connections(j,2)-1)*siz+(1:siz),anc_s];
            pa=[{phi,weights(j)*t/step_num}];
            gate_steps={gate_steps{:},init,center,outit};
            index_steps={index_steps{:},in,in,in};
            param_steps={param_steps{:},pa,pa,pa};
        end
    end

elseif length(t)==0 && length(step_num)==1
    param_steps={};
    for i=1:step_num
        for j=[1:2:length(weights) 2:2:length(weights)]
            in=[(connections(j,1)-1)*siz+(1:siz),(connections(j,2)-1)*siz+(1:siz),anc_s];
            pa=[{phi,weights(j)*t1/step_num}];
            gate_steps={gate_steps{:},init,center,outit};
            index_steps={index_steps{:},in,in,in};
            param_steps={param_steps{:},pa,pa,pa};
        end
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