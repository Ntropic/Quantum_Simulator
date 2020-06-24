function [rho sizes anc_sizes] = Noisy_Expansion2Matrix(rho,errors,gate,plotting)
%Noisy_Expansion2Matrix creates a matrix out of the elementary gate 
%expansion defined in exp_steps of the "gate" using a markovian error model
%Check if gates have been prepared?!
%       ->  Prepare Gates via "Gates_Tables_Prepare"
%       ->  current_param.twirling '0','','cz'
%       ->                      errors=[t_step,T1,T2 , phi,delta,E1] 
%       ->          the first 3 parameters in errrors refer to ''
%       ->          the 4rd,5th and 6th refer to errors of type 'cz'
%Errors are not assuming parallelizability of operations yet!
if nargin==3
    plotting=1;
end
if nargin<3
    error('more parameters needed')
end
fac=500; %How often do we have to Print the remaining time, progress and elapsed time

sizes=gate.size;
anc_sizes=gate.anc_size;
N=sizes+anc_sizes;


index=[gate.exp_steps.index];
matrix=[gate.exp_steps.matrix];
par=[gate.exp_steps.param];
param=unique([gate.exp_steps.param{:}]);
global_phase=gate.exp_steps.global_phase;
gate_twirling=gate.exp_steps.gate_twirling;

t_step=errors(1);
T1=errors(2);
T2=errors(3);
phi=errors(4);
delta=errors(5);
E1=errors(6);

K=Efficient_Kraus_Twirling( N,t_step,T1,T2 );
[K2,p2]=CZ_Twirling(phi,delta,E1);

%% Main Procedure
    if length(param)>0
        error('Not all parameters have been substituted')
    end
    %Embed Gates iteratively
    if isa(matrix{1},'function_handle')==0
        if length(global_phase)>0
            anc_matrix=Embed_Gate(matrix{1},index{1},N)*global_phase;
        else
            anc_matrix=Embed_Gate(matrix{1},index{1},N);
        end
    else
        if length(global_phase)>0
            anc_matrix=Embed_Gate(sym(matrix{1}),index{1},N)*global_phase;
        else
            anc_matrix=Embed_Gate(sym(matrix{1}),index{1},N);
        end
    end
    rho=anc_matrix*rho*anc_matrix';
    %Determine twirling type and twirl apropriately
    if any(findstr(gate_twirling{1},'cz'))
        rho=CZ_Twirl_Step(rho,K2,p2,index{1});
        rho=Kraus_Twirl_Step(rho,N,K);
    elseif length(gate_twirling{1})==0
        rho=Kraus_Twirl_Step(rho,N,K);
    end
    
    
    nums=floor(length(matrix)/fac);
    lens=length(matrix);
    if nums<1
        for i=2:length(matrix)
            if isa(matrix{i},'function_handle')==0
                anc_matrix=Embed_Gate(matrix{i},index{i},N);
            else
                anc_matrix=Embed_Gate(sym(matrix{i}),index{i},N);
            end
           	%rho=anc_matrix*rho*anc_matrix';
            %Determine twirling type and twirl apropriately
            if any(findstr(gate_twirling{1},'cz'))
                rho=CZ_Twirl_Step(rho,K2,p2,index{i},anc_matrix);
                rho=Kraus_Twirl_Step(rho,N,K,anc_matrix);
            elseif length(gate_twirling{i})==0
                rho=Kraus_Twirl_Step(rho,N,K,anc_matrix);
            end
        end
    
    else
        timer=tic();
        for i=2:fac
            if isa(matrix{i},'function_handle')==0
                anc_matrix=Embed_Gate(matrix{i},index{i},N);
            else
                anc_matrix=Embed_Gate(sym(matrix{i}),index{i},N);
            end
            %rho=anc_matrix*rho*anc_matrix';
            %Determine twirling type and twirl apropriately
            if any(findstr(gate_twirling{i},'cz'))
                rho=CZ_Twirl_Step(rho,K2,p2,index{i},anc_matrix);
                rho=Kraus_Twirl_Step(rho,N,K);
            elseif length(gate_twirling{i})==0
                rho=Kraus_Twirl_Step(rho,N,K,anc_matrix);
            end
        end
        t=toc(timer);
        if plotting==1
            fprintf(['    -> Remaining time: ' num2str(lens/fac*t-t) 's, Elapsed: ' num2str(t) 's, ' num2str(fac) '/' num2str(lens) '\n']);
        end
        for j=2:nums
            for i=fac*(j-1)+1:fac*j
                if isa(matrix{i},'function_handle')==0
                    anc_matrix=Embed_Gate(matrix{i},index{i},N);
                else
                    anc_matrix=Embed_Gate(sym(matrix{i}),index{i},N);
                end
               	rho=anc_matrix*rho*anc_matrix';
                %Determine twirling type and twirl apropriately
                if any(findstr(gate_twirling{i},'cz'))
                    rho=CZ_Twirl_Step(rho,K2,p2,index{i});
                    rho=Kraus_Twirl_Step(rho,N,K);
                elseif length(gate_twirling{i})==0
                    rho=Kraus_Twirl_Step(rho,N,K);
                end
            end
            t=toc(timer);
            if plotting==1
                fprintf(['    -> Remaining time: ' num2str(lens/(fac*j)*t-t) 's, Elapsed: ' num2str(t) 's, ' num2str(fac*j) '/' num2str(lens) '\n']);
            end
        end
        for i=fac*nums+1:length(matrix)
            if isa(matrix{i},'function_handle')==0
                anc_matrix=Embed_Gate(matrix{i},index{i},N);
            else
                anc_matrix=Embed_Gate(sym(matrix{i}),index{i},N);
            end
            rho=anc_matrix*rho*anc_matrix';
            %Determine twirling type and twirl apropriately
            if any(findstr(gate_twirling{i},'cz'))
                rho=CZ_Twirl_Step(rho,K2,p2,index{i});
                rho=Kraus_Twirl_Step(rho,N,K);
            elseif length(gate_twirling{i})==0
                rho=Kraus_Twirl_Step(rho,N,K);
            end
        end
        t=toc(timer);
        if plotting==1
            fprintf(['    -> Remaining time: 0s, Elapsed: ' num2str(t) 's, ' num2str(lens) '/' num2str(lens) '\n']);
        end
    end
    
end