function [ gate ] = Create_Comp_Gate( name,sizle,anc_size,gate_steps,index_steps,param_steps,circ,circuit_subs  )
%Create_Comp_Gate creates a new gate
if nargin==5
    param_steps=[];
end

if length(gate_steps)==0
    
    if nargin>6
        gate=Create_Empty_Comp_Gate( name,sizle,anc_size,circ );
    else
        gate=Create_Empty_Comp_Gate( name,sizle,anc_size);
    end
    
else

    gate=struct('names',name);
    gate.size=sizle;
    gate.anc_size=anc_size;
    gate.step_num=length(gate_steps);

    gate.steps.index=index_steps;

    max_ind=max([index_steps{:}]);
    if max_ind>anc_size+sizle
        fprintf('Gate uses more ancilla qubits, corrected description of gate (Create_Comp_Gate)');
        gate.anc_size=max([index_steps{:}])-sizle;
    end

    gate.steps.gates=gate_steps;

    if length(gate_steps)~=length(index_steps)
        error('Gate steps and Index Steps must be of same length! No output was created.');
    end

    if length(param_steps)>0
        gate.steps.param=param_steps;

        if length(gate_steps)~=length(param_steps)
            error('Gate steps and Parameter Steps must be of same length! No output was created.');
        end
        
        gate.steps.subs_names={};
        for j=1:length(param_steps)
            if length(param_steps{j})>0
                subs_name=gate_steps{j};
                for i=1:length(param_steps{j})
                    if isa(param_steps{j},'cell')
                        param_j_i=param_steps{j}{i};
                    else
                        param_j_i=param_steps{j}(i);
                    end
                    if isa(param_j_i,'sym')
                        subs_name=[subs_name '_' char(param_j_i)];
                    elseif isnumeric(param_j_i)
                        subs_name=[subs_name '_' num2str(param_j_i)];
                    else
                        subs_name=[subs_name '_' param_j_i];
                    end
                end
                gate.steps.subs_names{j}=subs_name;
            else
                gate.steps.subs_names{j}='';
            end
        end        
    else
        gate.steps.param={[]};
        for i=2:length(index_steps)
            gate.steps.param={gate.steps.param{:},[]};
        end
    end
    

    if nargin>6
        if length(circ)>0
            gate.circuit=circ;
        end
        if nargin>7
            gate.circuit_subs=circuit_subs;
        end
    end
    
end
end

