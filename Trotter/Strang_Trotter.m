function [ gate, step_gate, name, step_name ] = Strang_Trotter( t,step_num,name,circ_string,circ_string_step,sizle,anc_size,gates_names,gates_indexes,gates_param )
%TROTTER does the Lie-Suzuki-Trotter decomposition 
%Create sequence (seq) and step width d(seq)
%Pattern:
% A/2+B/2+C/2+D+C/2+B/2+A+B/2+C/2+D+C/2+B/2+A/2
%    |_________________| |_________________|
%        =step_gate          =step_gate

step_name=[name '_step'];

%no_step_gate if length(gates_names)==2 only when it's bigger
if length(gates_names)>2
    if length(gates_names)>1
        if nargin==9
            step_gate=Strang_Splitting(step_name,circ_string_step,sizle,anc_size,gates_names(2:end),gates_indexes(2:end));
        elseif nargin==10
            step_gate=Strang_Splitting(step_name,circ_string_step,sizle,anc_size,gates_names(2:end),gates_indexes(2:end),gates_param);
        end
    else
        %Only one gate 
        step_gate=Gate_by_Name(gates_names);
    end

    gate_steps={};
    index_steps={};
    param_steps={};
    if length(t)==1 && length(step_num)==1
        phi=sym('phi');
        gate_steps={gate_steps{:},gates_names{1}};
        index_steps={index_steps{:},gates_indexes{1}};
        param_steps={param_steps{:},[{phi,t/step_num/2}]};
        for i=1:step_num
            gate_steps={gate_steps{:},step_gate.names,gates_names{1}};
            index_steps={index_steps{:},[1:sizle+anc_size],gates_indexes{1}};
            param_steps={param_steps{:},[{phi,t/step_num}],[{phi,t/step_num}]};
        end
        gate_steps={gate_steps{:},gates_names{1}};
        index_steps={index_steps{:},gates_indexes{1}};
        param_steps={param_steps{:},[{phi,t/step_num/2}]};
    elseif length(t)==0 && length(step_num)==1
        phi=sym('phi');
        t1=sym('t');
        gate_steps={gate_steps{:},gates_names{1}};
        index_steps={index_steps{:},gates_indexes{1}};
        param_steps={param_steps{:},[{phi,t1/step_num/2}]};
        for i=1:step_num
            gate_steps={gate_steps{:},step_gate.names,gates_names{1}};
            index_steps={index_steps{:},[1:sizle+anc_size],gates_indexes{1}};
            param_steps={param_steps{:},[{phi,t1/step_num}],[{phi,t1/step_num}]};
        end
        gate_steps={gate_steps{:},gates_names{1}};
        index_steps={index_steps{:},gates_indexes{1}};
        param_steps={param_steps{:},[{phi,t1/step_num/2}]};
    else
        error('Step number needs to be specified, and time can maximally be length 1.');
    end
else %prepare step_gate development

    if step_num>=1
        gate_steps={};
        index_steps={};
        param_steps={};
        if length(t)==1 && length(step_num)==1
            phi=sym('phi');
            gate_steps={gate_steps{:},gates_names{1},gates_names{2}};
            index_steps={index_steps{:},gates_indexes{1},gates_indexes{2}};
            if nargin==10 && length(gates_param)>0
                param_steps={param_steps{:},[{phi,t/step_num/2*gates_param{1}}],[{phi,t/step_num*gates_param{2}}]};
            else
                param_steps={param_steps{:},[{phi,t/step_num/2}],[{phi,t/step_num}]};
            end
            for i=2:step_num
                gate_steps={gate_steps{:},gates_names{1},gates_names{2}};
                index_steps={index_steps{:},gates_indexes{1},gates_indexes{2}};
            	if nargin==10 && length(gates_param)>0
                    param_steps={param_steps{:},[{phi,t/step_num*gates_param{1}}],[{phi,t/step_num*gates_param{2}}]};
                else
                    param_steps={param_steps{:},[{phi,t/step_num}],[{phi,t/step_num}]};
                end
            end
            gate_steps={gate_steps{:},gates_names{1}};
            index_steps={index_steps{:},gates_indexes{1}};
            if nargin==10 && length(gates_param)>0
                param_steps={param_steps{:},[{phi,t/step_num/2*gates_param{1}}]};
            else
                param_steps={param_steps{:},[{phi,t/step_num/2}]};
            end
        elseif length(t)==0 && length(step_num)==1
            phi=sym('phi');
            t1=sym('t');
            gate_steps={gate_steps{:},gates_names{1},gates_names{2}};
            index_steps={index_steps{:},gates_indexes{1},gates_indexes{2}};
            if nargin==10 && length(gates_param)>0
                param_steps={param_steps{:},[{phi,t1/step_num/2*gates_param{1}}],[{phi,t1/step_num*gates_param{2}}]};
            else
                param_steps={param_steps{:},[{phi,t1/step_num/2}],[{phi,t1/step_num}]};
            end
            for i=2:step_num
                gate_steps={gate_steps{:},gates_names{1},gates_names{2}};
                index_steps={index_steps{:},gates_indexes{1},gates_indexes{2}};
                if nargin==10 && length(gates_param)>0
                    param_steps={param_steps{:},[{phi,t1/step_num*gates_param{1}}],[{phi,t1/step_num*gates_param{2}}]};
                else
                    param_steps={param_steps{:},[{phi,t1/step_num}],[{phi,t1/step_num}]};
                end
            end
            gate_steps={gate_steps{:},gates_names{1}};
            index_steps={index_steps{:},gates_indexes{1}};
            if nargin==10 && length(gates_param)>0
                param_steps={param_steps{:},[{phi,t1/step_num/2*gates_param{1}}]};
            else
                param_steps={param_steps{:},[{phi,t1/step_num/2}]};
            end
        else
            error('Step number needs to be specified, and time can maximally be length 1.');
        end
    else
        error('step_num must be >=1');
    end
    
    step_gate=[];
end

gate=Create_Comp_Gate(name,sizle,anc_size,gate_steps,index_steps,param_steps);
gate=Generate_Gate_Circuit(gate,circ_string,1:sizle,[],[]);
end