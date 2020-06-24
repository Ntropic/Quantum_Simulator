function [phi,t2] = DepthExpansion2Matrix_Return_Wave(phi0,gate,depth,single_error,double_error,plotting)
%Expansion2Matrix_Return_Wave creates a matrix out of the elementary gate 
%expansion defined in exp_steps of the "gate"
%Check if gates have been prepared?!
%       ->  Prepare Gates via "Gates_Tables_Prepare"
if nargin==2
    depth=-1;
end
if nargin<=3
    single_error=1;
    double_error=10;
    plotting=1;
end
if nargin==5
    plotting=1;
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

%% Find the gate groups up to depth
len=length(gate.exp_steps.depth);
grouping=zeros(len,depth);
for i=1:len
    mini=min([depth,length(gate.exp_steps.depth{i})]);
    grouping(i,1:mini)=gate.exp_steps.depth{i}(1:mini);
end
g=diff(grouping,1,1);
h=[any(g,2);1];
inds=find(h);

phi=zeros(length(phi0),length(inds)+1);
phi(:,1)=phi0;
phi2=phi0;
t2=zeros(1,length(inds)+1);
indis=[1;inds];
for j=1:length(inds)
    t2(j+1)=t2(j);
    if j==1
        if length(matrix{i})==2
            t2(j+1)=t2(j+1)+single_error;
        else
            t2(j+1)=t2(j+1)+double_error;
        end
    end
    for i=indis(j)+1:indis(j+1)
        if length(matrix{i})==2
            t2(j+1)=t2(j+1)+single_error;
        else
            t2(j+1)=t2(j+1)+double_error;
        end
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
    phi2=anc_matrix*phi2;
    if any(1==inds)
        a=find(1==inds);
        phi(:,a+1)=phi2;
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
            phi2=anc_matrix*phi2;
            if any(i==inds)
                a=find(i==inds);
                phi(:,a+1)=phi2;
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
            phi2=anc_matrix*phi2;
            if any(i==inds)
                a=find(i==inds);
                phi(:,a+1)=phi2;
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
                phi2=anc_matrix*phi2;
                if any(i==inds)
                    a=find(i==inds);
                    phi(:,a+1)=phi2;
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
            phi2=anc_matrix*phi2;
            if any(i==inds)
                a=find(i==inds);
                phi(:,a+1)=phi2;
            end
        end
        t=toc(timer);
        if plotting==1
            fprintf(['    -> Remaining time: 0s, Elapsed: ' num2str(t) 's, ' num2str(lens) '/' num2str(lens) '\n']);
        end
    end
    
end