function [phi,t2] = Expansion2Matrix_Return_Wave(phi0,gate,single_error,double_error)
%Expansion2Matrix_Return_Wave creates a matrix out of the elementary gate 
%expansion defined in exp_steps of the "gate"
%Check if gates have been prepared?!
%       ->  Prepare Gates via "Gates_Tables_Prepare"
if nargin==2
    single_error=1;
    double_error=10;
end

fac=500;

sizes=gate.size;
anc_sizes=gate.anc_size;
N=sizes+anc_sizes;


index=[gate.exp_steps.index];
matrix=[gate.exp_steps.matrix];
par=[gate.exp_steps.param];
param=unique([gate.exp_steps.param{:}]);
global_phase=gate.exp_steps.global_phase;

phi=zeros(length(phi0),length(matrix)+1);
phi(:,1)=phi0;
phi_n=phi0;
t2=zeros(1,length(matrix)+1);
for i=1:length(matrix)
    if length(matrix{i})==2
        t2(i+1)=t2(i)+single_error;
    else
        t2(i+1)=t2(i)+double_error;
    end
end

    if length(param)>0
        warning('Not all parameters have been substituted')
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
    phi(:,2)=anc_matrix*phi(:,1);
    
    nums=floor(length(matrix)/fac);
    lens=length(matrix);
    if nums<1
        for i=2:length(matrix)
            if isa(matrix{i},'function_handle')==0
                anc_matrix=Embed_Gate(matrix{i},index{i},N);
            else
                anc_matrix=Embed_Gate(sym(matrix{i}),index{i},N);
            end
            phi(:,i+1)=anc_matrix*phi(:,i);
        end
    else
        timer=tic();
        for i=2:fac
            if isa(matrix{i},'function_handle')==0
                anc_matrix=Embed_Gate(matrix{i},index{i},N);
            else
                anc_matrix=Embed_Gate(sym(matrix{i}),index{i},N);
            end
            phi(:,i+1)=anc_matrix*phi(:,i);
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
                phi(:,i+1)=anc_matrix*phi(:,i);
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
            phi(:,i+1)=anc_matrix*phi(:,i);
        end
        t=toc(timer);
        if plotting==1
            fprintf(['    -> Remaining time: 0s, Elapsed: ' num2str(t) 's, ' num2str(lens) '/' num2str(lens) '\n']);
        end
    end
    
end