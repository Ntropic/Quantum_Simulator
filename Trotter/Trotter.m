function [ gate , step_gate , name, step_name ] = Trotter( t,step_num,name,circ_string,circ_string_step,sizle,anc_size,gates_names,gates_indexes,gates_param )
%TROTTER does the Lie-Suzuki-Trotter decomposition 
%Create sequence (seq) and step width d(seq)
step_name=[name '_step'];

if length(name)==0
    name=['Gray_s_' num2str(step_num) '_trott'];
end
if length(circ_string)==0
    circ_string=['G_{s=' num2str(step_num) '}'];
end

if length(gates_names)>1 && isa(gates_names,'str')==0
    if nargin==9
        step_gate=Trotter_Step(step_name,circ_string_step,sizle,anc_size,gates_names,gates_indexes);
    elseif nargin==10
        step_gate=Trotter_Step(step_name,circ_string_step,sizle,anc_size,gates_names,gates_indexes,gates_param);
    end
else
    %Only one gate 
    if length(gates_names)==1
        step_gate.names=gates_names{1};
    else
	    step_gate.names=gates_names;
    end
end

gate_steps={};
index_steps={};
param_steps={};
if length(t)==1 && length(step_num)==1
    for i=1:step_num
        gate_steps={gate_steps{:},step_gate.names};
        index_steps={index_steps{:},[1:sizle+anc_size]};
        param_steps={param_steps{:},[{sym('phi'),t/step_num}]};
    end
elseif length(t)==0 && length(step_num)==1
    phi=sym('phi');
    t1=sym('t');
    for i=1:step_num
        gate_steps={gate_steps{:},step_gate.names};
        index_steps={index_steps{:},[1:sizle+anc_size]};
        param_steps={param_steps{:},[{phi,t1/step_num}]};
    end
else
    error('Step number needs to be specified, and time can maximally be length 1.');
end

gate=Create_Comp_Gate(name,sizle,anc_size,gate_steps,index_steps,param_steps);
gate=Generate_Gate_Circuit(gate,circ_string,1:sizle,[],[]);
end