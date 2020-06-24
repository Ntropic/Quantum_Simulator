function [ gate ] = MWG_Higher_Trotter_InOut( t,l,pattern,connections,weights,name,circ_string,sizle,anc_sizle,gates_names  )
%Higher_Trotter generates higher order Trotter decomposition formulas
% It requires
%       steps   -   to determine the higher order sheme (3 or 5 sheme)
%       gates_names={init,outit,center,connect gates}
% It outputs
%       seq     -   time sequence of the integration sheme
%       gate    -   new composite gate for one higher order Trotter step
% l=:Trotter steps

%Create sequence (seq) and step width d(seq)
[seq,d_seq]=Higher_Trotter_Time_Steps(pattern);

if length(connections)==1
    WG=connections;
    if length(weights)==1
        weights=weights*ones(WG-1,1);
    elseif length(weights)==0
        weights=ones(WG-1,1);
    end
    connections=[(1:WG-1)',(2:WG)'];
else
    WG=max(connections(:));
end
if length(weights)==1
    weights=weights*ones(WG-1,1);
elseif length(weights)==0
    weights=ones(WG-1,1);
end

siz=sizle/2;
if siz~=round(siz)
    error('Size needs to be even qubits big.')
end

s=siz*WG;
anc_s=(1:anc_sizle)+s;

init=gates_names{1};
init2=gates_names{2};
outit=gates_names{3};
outit2=gates_names{4};
symit=gates_names{5};
con=gates_names{6};

phi=sym('phi');
if length(t)==1 && length(l)==1
    t=t;
else
    t1=sym('t');
    t=t1;
end
l_d_s=length(d_seq);

gate_steps={};
index_steps={};
param_steps={};
counter=0;
step_num=l;
if length(weights)>1
    for k=1:l
        for i=1:l_d_s
            counter=counter+1;
            if counter==1
                j=1;
                in=[(connections(j,1)-1)*siz+(1:siz),(connections(j,2)-1)*siz+(1:siz),anc_s];
                pa=[{phi,weights(j)*t*d_seq(i)/step_num}];
                if j~=length(weights)
                    gate_steps={gate_steps{:},init};
                else
                    gate_steps={gate_steps{:},symit};
                end
                index_steps={index_steps{:},in};
                param_steps={param_steps{:},pa};
            end
            for j=2:length(weights)
                in=[(connections(j,1)-1)*siz+(1:siz),(connections(j,2)-1)*siz+(1:siz),anc_s];
                pa=[{phi,weights(j)*t*d_seq(i)/step_num}];
                if j~=length(weights)
                    gate_steps={gate_steps{:},init};
                else
                    gate_steps={gate_steps{:},symit};
                end
                index_steps={index_steps{:},in};
                param_steps={param_steps{:},pa};
            end
            for j=length(weights)-1:-1:2
                in=[(connections(j,1)-1)*siz+(1:siz),(connections(j,2)-1)*siz+(1:siz),anc_s];
                pa=[{phi,weights(j)*t*d_seq(i)/step_num}];
                gate_steps={gate_steps{:},outit};
                index_steps={index_steps{:},in};
                param_steps={param_steps{:},pa};
            end
            if counter~=step_num*l_d_s
                j=1;
                in=[(connections(j,1)-1)*siz+(1:siz),(connections(j,2)-1)*siz+(1:siz),anc_s];
                i_1=i+1;
                if i_1>l_d_s
                    i_1=1;
                end
                pa=[{phi,weights(j)*t*d_seq(i)/step_num}];
                pa2=[{phi,weights(j)*t*d_seq(i_1)/step_num}];
                pa3=[{phi,weights(j)*t*(d_seq(i)+d_seq(i_1))/step_num}];
                if j~=length(weights)
                    gate_steps={gate_steps{:},init2,con,outit2};
                    index_steps={index_steps{:},in,in,in};
                    param_steps={param_steps{:},pa,pa3,pa2};
                elseif counter>1
                    gate_steps={gate_steps{:},symit};
                    index_steps={index_steps{:},in};
                    param_steps={param_steps{:},pa};
                end

            else
                j=1;
                in=[(connections(j,1)-1)*siz+(1:siz),(connections(j,2)-1)*siz+(1:siz),anc_s];
                pa=[{phi,weights(j)*t*d_seq(i)/step_num}];
                if j~=length(weights)
                    gate_steps={gate_steps{:},outit};
                    index_steps={index_steps{:},in};
                    param_steps={param_steps{:},pa};
                elseif counter>1
                    gate_steps={gate_steps{:},symit};
                    index_steps={index_steps{:},in};
                    param_steps={param_steps{:},pa};
                end
            end
        end
    end
else %2 mode case -> only uneven mode interaction exists
    j=1;
    counter=0;
    step_num=l;
    for k=1:l
        for i=1:l_d_s %Elements of the higher order scheme
            counter=counter+1;
            if counter==1 %Initial
                in=[(1:siz),siz+(1:siz),anc_s];
                pa=[{phi,weights(j)*t*d_seq(i)/step_num}];
                gate_steps={gate_steps{:},init,symit};
                index_steps={index_steps{:},in,in};
                param_steps={param_steps{:},pa,pa};
            end
            if counter~=step_num*l_d_s %Subsequent
                i_1=i+1; %Next element
                if i_1>l_d_s
                    i_1=1;
                end
                in=[(1:siz),siz+(1:siz),anc_s];
                pa=[{phi,weights(j)*t*(d_seq(i)+d_seq(i_1))/step_num}];
                pa2=[{phi,weights(j)*t*d_seq(i_1)/step_num}];
                gate_steps={gate_steps{:},con,symit};
                index_steps={index_steps{:},in,in};
                param_steps={param_steps{:},pa,pa2};
                
                if j~=length(weights)
                    gate_steps={gate_steps{:},outit2,con,init2};
                    index_steps={index_steps{:},in,in,in};
                    param_steps={param_steps{:},pa,pa3,pa2};
                elseif counter>1
                    gate_steps={gate_steps{:},symit};
                    index_steps={index_steps{:},in};
                    param_steps={param_steps{:},pa};
                end
            else %Final
                in=[(1:siz),siz+(1:siz),anc_s];
                pa=[{phi,weights(j)*t*d_seq(i)/step_num}];
                gate_steps={gate_steps{:},outit};
                index_steps={index_steps{:},in};
                param_steps={param_steps{:},pa};
            end
        end
    end
end

if length(name)==0
    name=['H_{N=' num2str(WG) ',n=' num2str(sizle) '}^{(' num2str(length(step_num)) ')}(' latex(sym(t)) ')'];
end
if length(circ_string)==0
    circ_string=['H_{N=' num2str(WG) ',n=' num2str(sizle) '}^{(' num2str(length(step_num)) ')}(' latex(sym(t)) ')'];
end
gate=Create_Comp_Gate(name,s,anc_sizle,gate_steps,index_steps,param_steps);
gate=Generate_Gate_Circuit(gate,circ_string,1:s,[],[]);
end

