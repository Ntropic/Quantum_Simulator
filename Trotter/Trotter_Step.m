function [ gate ] = Trotter_Step( name,circ_string,size,anc_size,gates_names,gates_indexes,gates_param )
%Trotter_Step creates a trotter step from operators given by gates_names,
%and gates_indexes
%   -> Operators Uij need to be dependent on sym('phi'), so that this can be
%   -> substituted by sym('t') or sym('t/2')

if nargin==6
    gate_steps={gates_names{1:end}};
    index_steps={gates_indexes{1:end}};

    param_steps={};
    for i=1:length(gates_names)

        param_steps={param_steps{:},[]};
    end

    gate=Create_Comp_Gate(name,size,anc_size,gate_steps,index_steps,param_steps);
    gate=Generate_Gate_Circuit(gate,circ_string,1:size,[],[]);
elseif nargin==7
    gate_steps={gates_names{1:end}};
    index_steps={gates_indexes{1:end}};
    param_steps={gates_param{1:end}};

    gate=Create_Comp_Gate(name,size,anc_size,gate_steps,index_steps,param_steps);
    gate=Generate_Gate_Circuit(gate,circ_string,1:size,[],[]);
end
end