function [ cost_maps,extrap_maps,trott_steps_maps,s_q_maps ] = Fastest_Fidelity( fid , F , s_q , costs , F_predictors, s_q_predictors)
%FASTEST_FIDELITY How fast can one achieve a desired fidelity (fid) with 
%certain parameters. Finds the algorithm with the lowest cost. 
%(costs are given for single and two qubit gates) 
% F  ={F_1,F_2,...}
% s_q={s_q_1,s_q_2,...}
%      -> The dimensions (of F{i} and s_q{i} are [#WG's,#Photons 
%
% costs=[single qubit cost, double qubit cost]
% ranges=[N_m(1),N_m(end),n_m(1),n_m(end),steps];


%How many datasets are compared
how_many=length(F);
if length(s_q)~=length(F)
    error('F and s_q need to be of equal lengths');
end
%Determine size of simulation data
sizes=size(F{1});
for i=2:how_many
    if sizes~=size(F{i})
        error('Datasets (F{i}) must have equal sizes.')
    end
    for j=1:length(sizes)
    if sizes(j)~=size(s_q{i},j)
        error('Datasets (s_q{i}) and (F{i}) must have equal sizes, allthough s_q has an additional size parameter, for single and double qubits.')
    end
    end

f_goal=1-fid;   %Error corresponding to fidelity rate

%% Go through all parameters, to find the number of trotter steps needed to get below a certain error rate/above a certain fidelity
cost_map_z=zeros(sizes([1 2]));
cost_maps=cell(how_many,1);
trott_steps_map_z=zeros(sizes([1 2]));
trott_steps_maps=cell(how_many,1);
s_q_map_z=zeros([sizes([1 2]) 2]);
s_q_maps=cell(how_many,1);
extrap_maps_z=zeros(sizes([1 2]));
extrap_maps=cell(how_many,1);
for i=1:how_many
    cost_mapper=cost_map_z; %Inititalize with zeros
    trott_steps_mapper=trott_steps_map_z;
    s_q_mapper=s_q_map_z;
    extrap_mapper=extrap_maps_z;
    F_i=F{i};
    F_coeff=F_predictors{1,i};
    F_fn=F_predictors{2,i};
    F_n=F_predictors{3,i};
    s_q_n=s_q_predictors{i};
    
    s_q_i=s_q{i};
    for j=1:sizes(1) %#of WG's
        for k=1:sizes(2) %#of Photons
            found=0;                    %Has the lowest number of Trotter steps for this fidelity been found?
            for l=1:sizes(3) %Trotter steps
                if F_i(j,k,l)>=fid
                    if found==0
                        found=1; %-> found at l trotter steps
                        cost_mapper(j,k)=s_q_i(j,k,l,1)*costs(1)+s_q_i(j,k,l,2)*costs(2);
                        trott_steps_mapper(j,k)=l;
                        s_q_mapper(j,k,:)=[s_q_i(j,k,l,1) s_q_i(j,k,l,2)];
                        extrap_mapper(j,k)=0;
                    end
                end
            end
            if found==0 %Still not found? -> Extrapolate from data
                %Extrapolate convergence order -> predict error rate
                c=F_coeff(j,k);
                fn=F_fn(j,k);
                n=F_n(j,k);
                m_goal=ceil(n*(fn/f_goal)^(1/c));
                s_q_now=[s_q_i(j,k,end,1)*s_q_n(j,k,1)*(m_goal-sizes(3)),s_q_i(j,k,end,2)*s_q_n(j,k,2)*(m_goal-sizes(3))];
                cost_mapper(j,k)=s_q_now(1)*costs(1)+s_q_now(2)*costs(2);
                trott_steps_mapper(j,k)=m_goal;
                s_q_mapper(j,k,:)=s_q_now;
                extrap_mapper(j,k)=m_goal-sizes(3);
            end
        end
    end
    cost_maps{i}=cost_mapper;
    trott_steps_maps{i}=trott_steps_mapper;
    s_q_maps{i}=s_q_mapper;
    extrap_maps{i}=extrap_mapper;
end
end

