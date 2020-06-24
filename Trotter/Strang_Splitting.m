function [ gate ] = Strang_Splitting( name,circ_string,sizle,anc_size,gates_names,gates_indexes,gates_param)
%STRANG_SPLITTING splits a set of operators 'gates_names' symmetrically, 
%erasing the second order correction terms 
%   -> Operators Uij need to be dependent on sym('phi'), so that this can be
%   -> substituted by sym('t') or sym('t/2')

if nargin==6
    gate_steps={gates_names{1:end} gates_names{end-1:-1:1}};
    index_steps={gates_indexes{1:end} gates_indexes{end-1:-1:1}};

    param_steps={};
    phi=sym('phi');
    for i=1:length(gates_names)-1
        param_steps={param_steps{:},[{phi,phi/2}]};
    end
    param_steps={param_steps{:},[]};
    for i=1:length(gates_names)-1
        param_steps={param_steps{:},[{phi,phi/2}]};
    end

    gate=Create_Comp_Gate(name,sizle,anc_size,gate_steps,index_steps,param_steps);
    gate=Generate_Gate_Circuit(gate,circ_string,1:sizle,[],[]);
elseif nargin==7
    gate_steps={gates_names{1:end} gates_names{end-1:-1:1}};
    index_steps={gates_indexes{1:end} gates_indexes{end-1:-1:1}};

    param_steps={};
    phi=sym('phi');
    for i=1:length(gates_names)-1
        if length(gates_param{i})==0
            param_steps={param_steps{:},[{phi,phi/2}]};
        else
            param_steps={param_steps{:},[{phi,phi/2*gates_param{i}}]};
        end
    end
    param_steps={param_steps{:},[]};
    for i=1:length(gates_names)-1
      	if length(gates_param{i})==0
            param_steps={param_steps{:},[{phi,phi/2}]};
        else
            param_steps={param_steps{:},[{phi,phi/2*gates_param{i}}]};
        end
    end

    gate=Create_Comp_Gate(name,sizle,anc_size,gate_steps,index_steps,param_steps);
    gate=Generate_Gate_Circuit(gate,circ_string,1:sizle,[],[]);    
    
end
end