function [ seq , d_seq , gate , step_gate ] = Higher_Trotter( t,l,steps,name,gate_circ_string,circ_string_step,sizle,anc_size,gates_names,gates_indexes,gates_param )
%Higher_Trotter generates higher order Trotter decomposition formulas
% It requires
%       steps   -   to determine the higher order sheme (3 or 5 sheme)
%       S1      -   the Strang Splitting of the operators
% It outputs
%       seq     -   time sequence of the integration sheme
%       gate    -   new composite gate for one higher order Trotter step
%       S_i_gates - higher order Strang_splitting_operators as a gate list

%Create sequence (seq) and step width d(seq)
[seq,d_seq]=Higher_Trotter_Time_Steps(steps);

if length(name)>0
    step_name=[name '_step'];
else
    name=['H' num2str(steps(:)','_%i')];
    step_name=[name '_step'];
end

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

    step_gate_names=step_gate.names;

    in=[1:sizle+anc_size];

    gate_steps={};
    index_steps={};
    param_steps={};

    phi=sym('phi');
    t1=sym('t');

    l_d_s=length(d_seq);
    
    gate_steps={gate_steps{:},gates_names{1}};
    index_steps={index_steps{:},gates_indexes{1}};
    if length(t)==1
        param_steps={param_steps{:},[{phi,t/2*d_seq(1/l)}]};
    else
        param_steps={param_steps{:},[{phi,t1/2*d_seq(1)/l}]};
    end
    for j=1:l
        for i=1:l_d_s
            if i<l_d_s
                gate_steps={gate_steps{:},step_gate_names,gates_names{1}};
                index_steps={index_steps{:},in,gates_indexes{1}};
                if length(t)==1
                    param_steps={param_steps{:},[{phi,t*d_seq(i)/l}],[{phi,t/2*(d_seq(i)+d_seq(i+1))/l}]};
                else
                    param_steps={param_steps{:},[{phi,t1*d_seq(i)/l}],[{phi,t1/2*(d_seq(i)+d_seq(i+1))/l}]};
                end
            else
                if j<l
                    gate_steps={gate_steps{:},step_gate_names,gates_names{1}};
                    index_steps={index_steps{:},in,gates_indexes{1}};
                    if length(t)==1
                        param_steps={param_steps{:},[{phi,t*d_seq(i)/l}],[{phi,t1/2*(d_seq(i)+d_seq(1))/l}]};
                    else
                        param_steps={param_steps{:},[{phi,t1*d_seq(i)/l}],[{phi,t1/2*(d_seq(i)+d_seq(1))/l}]};
                    end
                else
                    gate_steps={gate_steps{:},step_gate_names};
                    index_steps={index_steps{:},in};
                    if length(t)==1
                        param_steps={param_steps{:},[{phi,t*d_seq(i)/l}]};
                    else
                        param_steps={param_steps{:},[{phi,t1*d_seq(i)/l}]};
                    end                    
                end
            end
        end
    end
    gate_steps={gate_steps{:},gates_names{1}};
    index_steps={index_steps{:},gates_indexes{1}};
    if length(t)==1
        param_steps={param_steps{:},[{phi,t/2*d_seq(end)/l}]};
    else
        param_steps={param_steps{:},[{phi,t1/2*d_seq(end)/l}]};
    end
else %prepare step_gate development (2 cunstituent ates)
    in=[1:sizle+anc_size];

    gate_steps={};
    index_steps={};
    param_steps={};

    phi=sym('phi');
    t1=sym('t');

    l_d_s=length(d_seq);
    
    gate_steps={gate_steps{:},gates_names{1}};
    index_steps={index_steps{:},gates_indexes{1}};
    if length(t)==1
        param_steps={param_steps{:},[{phi,t/2*d_seq(1)/l}]};
    else
        param_steps={param_steps{:},[{phi,t1/2*d_seq(1)/l}]};
    end
    for j=1:l
        for i=1:l_d_s
            if i<l_d_s
                gate_steps={gate_steps{:},gates_names{2},gates_names{1}};
                index_steps={index_steps{:},gates_indexes{2},gates_indexes{1}};
                if length(t)==1
                    param_steps={param_steps{:},[{phi,t*d_seq(i)/l}],[{phi,t/2*(d_seq(i)+d_seq(i+1))/l}]};
                else
                    param_steps={param_steps{:},[{phi,t1*d_seq(i)/l}],[{phi,t1/2*(d_seq(i)+d_seq(i+1))/l}]};
                end
            else
                if j<l
                    gate_steps={gate_steps{:},gates_names{2},gates_names{1}};
                    index_steps={index_steps{:},gates_indexes{2},gates_indexes{1}};
                    if length(t)==1
                        param_steps={param_steps{:},[{phi,t*d_seq(i)/l}],[{phi,t/2*(d_seq(i)+d_seq(1))/l}]};
                    else
                        param_steps={param_steps{:},[{phi,t1*d_seq(i)/l}],[{phi,t1/2*(d_seq(i)+d_seq(1))/l}]};
                    end
                else
                   	gate_steps={gate_steps{:},gates_names{2}};
                    index_steps={index_steps{:},gates_indexes{2}};
                    if length(t)==1
                        param_steps={param_steps{:},[{phi,t*d_seq(i)/l}]};
                    else
                        param_steps={param_steps{:},[{phi,t1*d_seq(i)/l}]};
                    end
                end
            end
        end
    end
    gate_steps={gate_steps{:},gates_names{1}};
    index_steps={index_steps{:},gates_indexes{1}};
    if length(t)==1
        param_steps={param_steps{:},[{phi,t/2*(d_seq(end))/l}]};
    else
        param_steps={param_steps{:},[{phi,t1/2*(d_seq(end))/l}]};
    end
    step_gate=[];
end
gate=Create_Comp_Gate(name,sizle,anc_size,gate_steps,index_steps,param_steps);
gate=Generate_Gate_Circuit(gate,gate_circ_string,1:sizle,[],[]);
end