function [ gate ] = Higher_Trotter_InOut( t,l,pattern,name,circ_string,size,anc_size,gates_names  )
%Higher_Trotter generates higher order Trotter decomposition formulas
% It requires
%       steps   -   to determine the higher order sheme (3 or 5 sheme)
%       gates_names={init,outit,center,connect gates}
% It outputs
%       seq     -   time sequence of the integration sheme
%       gate    -   new composite gate for one higher order Trotter step

%Create sequence (seq) and step width d(seq)
[seq,d_seq]=Higher_Trotter_Time_Steps(pattern);

init=gates_names{1};
outit=gates_names{2};
center=gates_names{3};
connect=gates_names{4};

in=[1:size+anc_size];

gate_steps={init};
index_steps={in};

phi=sym('phi');
t1=sym('t');

l_d_s=length(d_seq);

for j=1:l
    for i=1:l_d_s
        %Coefficients for different types of recursive approach
        if i<l_d_s
            if length(t)==1
                pa=[{phi,t*d_seq(i)/l}];
                pa2=[{phi,t*(d_seq(i)+d_seq(i+1))/2/l}];
            else
                pa=[{phi,t1*d_seq(i)/l}];
                pa2=[{phi,t1*(d_seq(i)+d_seq(i+1))/2/l}];
            end
            if i==1 && j==1
                param_steps={pa};
            end
            gate_steps={gate_steps{:},center,connect};
            index_steps={index_steps{:},in,in};
            param_steps={param_steps{:},pa,pa2};
        else
            if length(t)==1
                pa=[{phi,t*d_seq(i)/l}];
                pa2=[{phi,t*(d_seq(i)+d_seq(1))/2/l}];
            else
                pa=[{phi,t1*d_seq(i)/l}];
                pa2=[{phi,t1*(d_seq(i)+d_seq(1))/2/l}];
            end
            if i==1 && j==1
                param_steps={pa};
                gate_steps={gate_steps{:},center,outit};
                index_steps={index_steps{:},in,in};
                param_steps={param_steps{:},pa,pa};
            else
                if j==l
                    gate_steps={gate_steps{:},center,outit};
                    index_steps={index_steps{:},in,in};
                    param_steps={param_steps{:},pa,pa};
                else
                    gate_steps={gate_steps{:},center,connect};
                    index_steps={index_steps{:},in,in};
                    param_steps={param_steps{:},pa,pa2};                  
                end
            end
        end
    end
end

gate=Create_Comp_Gate(name,size,anc_size,gate_steps,index_steps,param_steps);
gate=Generate_Gate_Circuit(gate,circ_string,1:size,[],[]);
end

