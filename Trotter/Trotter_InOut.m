function [ gate ] = Trotter_InOut( t,step_num,name,circ_string,sizle,anc_sizle,gates_names )
%TROTTER_INOUT does the Lie-Suzuki-Trotter decomposition for
%gates_names={init,outit,center,connect gates}
%Create sequence (seq) and step width d(seq)

t1=sym('t');
phi=sym('phi');


init=gates_names{1};
outit=gates_names{2};
center=gates_names{3};
connect=gates_names{4};

in=[1:sizle+anc_sizle];
gate_steps={init};
index_steps={in};

if length(t)==1 && length(step_num)==1
    pa=[{phi,t/step_num}];
    %pa2=[{phi,t*2/step_num}];
    param_steps={pa};
    for i=1:step_num-1
        gate_steps={gate_steps{:},center,connect};
        index_steps={index_steps{:},in,in};
        param_steps={param_steps{:},pa,pa};
    end
    gate_steps={gate_steps{:},center,outit};
    index_steps={index_steps{:},in,in};
    param_steps={param_steps{:},pa,pa};
elseif length(t)==0 && length(step_num)==1
    pa=[{phi,t1/step_num}];
    %pa2=[{phi,t1*2/step_num}];
    param_steps={pa};
    for i=1:step_num-1
        gate_steps={gate_steps{:},center,connect};
        index_steps={index_steps{:},in,in};
        param_steps={param_steps{:},pa,pa};
    end
    gate_steps={gate_steps{:},center,outit};
    index_steps={index_steps{:},in,in};
    param_steps={param_steps{:},pa,pa};
end


if length(name)==0
    name=['T_{n=' num2str(sizle) '}^{(' num2str(length(step_num)) ')}(' latex(sym(t)) ')'];
end
if length(circ_string)==0
    circ_string=['T_{n=' num2str(sizle) '}^{(' num2str(length(step_num)) ')}(' latex(sym(t)) ')'];
end
gate=Create_Comp_Gate(name,sizle,anc_sizle,gate_steps,index_steps,param_steps);
gate=Generate_Gate_Circuit(gate,circ_string,1:sizle,[],[]);
end