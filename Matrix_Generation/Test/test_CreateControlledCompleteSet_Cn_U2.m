%test_CreateControlledCompleteSet_Cn_U2.m 
clc;
clear all;
close all;

p=pwd;
elem=strsplit(p,'\');
shortened=strjoin(elem((1:length(elem)-2)),'\');
addpath(genpath(shortened))

%Load gates
load('../../Gates_Table/elem_gates.mat','-mat')
load('../../Gates_Table/comp_gates.mat','-mat')

%Test CreateControlledCompleteSet_Cn_U2 for n
for n=1:5
    stringing=0;    %Strings of matrices?

    global alpha beta delta theta subs_comp
    alpha=sym('alpha');
    beta=sym('beta');
    delta=sym('delta');
    theta=sym('theta');
    assume(alpha,'real');
    assume(beta,'real');
    assume(delta,'real');
    assume(theta,'real');
    subs_comp={(exp(-(alpha*1i)/2)*exp(-(beta*1i)/2)*exp(delta*1i)*exp(-(theta*1i)/2))/2 + (exp(-(alpha*1i)/2)*exp(-(beta*1i)/2)*exp(delta*1i)*exp((theta*1i)/2))/2,- (exp(-(alpha*1i)/2)*exp((beta*1i)/2)*exp(delta*1i)*exp(-(theta*1i)/2)*1i)/2 + (exp(-(alpha*1i)/2)*exp((beta*1i)/2)*exp(delta*1i)*exp((theta*1i)/2)*1i)/2;
               (exp((alpha*1i)/2)*exp(-(beta*1i)/2)*exp(delta*1i)*exp(-(theta*1i)/2)*1i)/2 - (exp((alpha*1i)/2)*exp(-(beta*1i)/2)*exp(delta*1i)*exp((theta*1i)/2)*1i)/2,(exp((alpha*1i)/2)*exp((beta*1i)/2)*exp(delta*1i)*exp(-(theta*1i)/2))/2 + (exp((alpha*1i)/2)*exp((beta*1i)/2)*exp(delta*1i)*exp((theta*1i)/2))/2};


    %Generate Gates
    t=tic();
    [PU_gates,Cn]=CreateCompleteSet_Cn_U2(n,elem_gates,comp_gates,[1,1,0],1);
    t2=toc(t);
    fprintf('Created %1d gates in %1.2f s.\n',size(PU_gates,1)*size(PU_gates,2),t2)
    %Print Gates
    num_i=length(sprintf('%01d',size(PU_gates,1)));
    num_j=length(sprintf('%01d',size(PU_gates,2)));

    if stringing==1
        for i=1:size(Cn,2)
            [M index M_mode]=Sparse2Square(Cn(i).matrix,{'fock unit full dense',1,n+1},'substitution_list',{'U11','U12';'U21','U22'},'substitution_compare',subs_comp);
            fprintf(['\n' Cn(i).names ':\n'])
            fprintf(M_mode)
        end
        fprintf('---------------------------------------------------')

        stringer=['\nP matrix n=' num2str(n) ', i=%0' num2str(num_i) 'd, j=%0' num2str(num_j) 'd\n'];
        for i=1:size(PU_gates,1)
            for j=1:size(PU_gates,2)
                [M index M_mode]=Sparse2Square(PU_gates(i,j).matrix,{'fock unit full dense',1,n+1},'substitution_list',{'U11','U12';'U21','U22'},'substitution_compare',subs_comp);
                fprintf(stringer,i,j)
                fprintf([PU_gates(i,j).names ':\n'])
                fprintf(M_mode)
            end
        end
    end

    fprintf('Index list:\n')
    index_list=[];
    for i=1:size(PU_gates,1)
        for j=1:size(PU_gates,2)
            fprintf([PU_gates(i,j).names '  ' num2str(PU_gates(i,j).fock_index') '\n'])
            index_list=[index_list;PU_gates(i,j).fock_index'];
        end
    end
    index_list=sortrows(index_list,2);
    phi=linspace(0,pi,50);
    max_r=max(abs((index_list(:,1)-index_list(:,2))/2));
    cmap=jet(max_r*2);
    for i=1:length(index_list)
        b=index_list(i,:);
        m=mean(b);
        r=(max(b)-min(b))/2;
        x=m+r*cos(phi);
        y=r*sin(phi);
        plot(x,y,'Color',cmap(r*2,:));
        hold on;
    end
    axis equal;
    axis([0 max(index_list(:))+1 0 (max(index_list(:)))/4+0.5]);
    pause(0.1);
end