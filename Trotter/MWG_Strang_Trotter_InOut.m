function [ gate, name ] = MWG_Strang_Trotter_InOut( t,modes,weights,step_num,name,circ_string,sizle,anc_sizle,gates_names )
%MWG_TROTTER_InOut does the Lie-Suzuki-Trotter decomposition 
%Create sequence (seq) and step width d(seq)
%Pattern:
% A/2+B/2+C/2+D+C/2+B/2+A+B/2+C/2+D+C/2+B/2+A/2
%    |_________________| |_________________|
%        =step_gate          =step_gate

t1=sym('t');
phi=sym('phi');

siz=sizle/2;
if siz~=round(siz)
    error('Size needs to be even qubits big.')
end

WG=modes;
if length(weights)==1
    weights=weights*ones(WG-1,1);
elseif length(weights)==0
    weights=ones(WG-1,1);
end

s=siz*WG;
anc_s=(1:anc_sizle)+s;

ie=gates_names{1};
iu=gates_names{2};
fe=gates_names{3};
fu=gates_names{4};
u=gates_names{5};
e=gates_names{6};

% Decompose into even and uneven versions of all modes
%gate names
U_I={ie,fu};
U_F={iu,fe};
U_C={iu,e,fu};
E_C={ie,u,fe};

gate_steps={};
index_steps={};
if length(t)==1 && length(step_num)==1
    t=t;
else
    t=t1;
end
param_steps={};
%% Initial gates (Uneven)
if modes>2
    for j=1:2:modes-1
        gate_steps={gate_steps{:},U_I{:}};
        in=[(j-1)*siz+(1:siz),(j)*siz+(1:siz),anc_s];
        pa=[{phi,weights(j)*t/step_num}];
        index_steps={index_steps{:},in,in};
        param_steps={param_steps{:},pa,pa};
    end
    for i=1:step_num
        %% Even Center Piece
        for j=2:2:modes-1
            gate_steps={gate_steps{:},E_C{:}};
            in=[(j-1)*siz+(1:siz),(j)*siz+(1:siz),anc_s];
            pa=[{phi,weights(j)*t/step_num}];
            index_steps={index_steps{:},in,in,in};
            param_steps={param_steps{:},pa,pa,pa};
        end
        if i<step_num
            for j=1:2:modes-1
                gate_steps={gate_steps{:},U_C{:}};
                in=[(j-1)*siz+(1:siz),(j)*siz+(1:siz),anc_s];
                pa=[{phi,weights(j)*t/step_num}];
                index_steps={index_steps{:},in,in,in};
                param_steps={param_steps{:},pa,pa,pa};
            end
        end
        if i==step_num
            for j=1:2:modes-1
                gate_steps={gate_steps{:},U_F{:}};
                in=[(j-1)*siz+(1:siz),(j)*siz+(1:siz),anc_s];
                pa=[{phi,weights(j)*t/step_num}];
                index_steps={index_steps{:},in,in};
                param_steps={param_steps{:},pa,pa};
            end
        end
    end
else %2 mode case -> only uneven mode interaction exists
    gate_steps={gate_steps{:},iu,e};
    in=[(1:siz),siz+(1:siz),anc_s];
    pa=[{phi,weights(1)*t/step_num}];
    index_steps={index_steps{:},in,in};
    param_steps={param_steps{:},pa,pa};
    for i=1:step_num-1
  	    gate_steps={gate_steps{:},u,e};
        in=[(1:siz),siz+(1:siz),anc_s];
        pa=[{phi,weights(1)*t/step_num}];
        index_steps={index_steps{:},in,in};
        param_steps={param_steps{:},pa,pa};
    end
   	gate_steps={gate_steps{:},fu};
    in=[(1:siz),siz+(1:siz),anc_s];
    pa=[{phi,weights(1)*t/step_num}];
    index_steps={index_steps{:},in};
    param_steps={param_steps{:},pa};
end

if length(name)==0
    name=['S_{M=' num2str(WG) ',N\leq =' num2str(sizle) '}^{(' num2str(length(step_num)) ')}(' latex(sym(t)) ')'];
end
if length(circ_string)==0
    circ_string=['S_{M=' num2str(WG) ',N\leq ' num2str(sizle) '}^{(' num2str(length(step_num)) ')}(' latex(sym(t)) ')'];
end

gate=Create_Comp_Gate(name,s,anc_sizle,gate_steps,index_steps,param_steps);
gate=Generate_Gate_Circuit(gate,circ_string,1:s,[],[]);
end