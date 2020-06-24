%const_n_permutations_trot_steps_higher_distribution.m
%Permutates the order of the Elementary (Interaction) Operations in Trotter Expansion
clear all;
close all;
clc;

p=pwd;
if any(strfind(p,'\'));
    elem=strsplit(p,'\');
else
    elem=strsplit(p,'/');
end
clear p;
shortened=fullfile(elem{1:end-3});
addpath(genpath(shortened));

%Photon number
n_max=7;
how_many=100; %The first how many permutations (maximum)
n_bins=100;
t=pi/4;
plotting=1;
trott_type=0;   %0=trotter, 1=strang, 2,else=higher order
rec_order=0;

for n=2:n_max
    %% Create Gray ordering
    s=ceil(log2(n+1));

    %Create Hamiltonian and seperate into elementary interaction parts
    tree=bintree(s,1);
    gray=all_num_algorithm(tree);

    index=zeros(n,2);
    
    %% Simple Ordering
    k=1;
    prefactor=[];
    for l=n%1:n
        for i=1:l
            j=l-i+1;
            %prefactor=[prefactor;1+0*(i*j)];  %Strength of interaction 
            prefactor=[prefactor;sqrt(i*j)];

            gray_num=gray([i+1,j]);
            graybin=dec2bin(gray_num-1,s)';
            graybin=graybin(:)';
            graybin=graybin-'0';

            gray_num2=gray([i,j+1]);
            graybin2=dec2bin(gray_num2-1,s)';
            graybin2=graybin2(:)';
            graybin2=graybin2-'0';
            index(k,:)=[bin2dec(num2str(graybin)) bin2dec(num2str(graybin2))]+1;
            n_i(k)=i-1+j; %Total particles
            m_i(k,:)=[i-1,j,i,j-1]; %Before -> after (each with occupation of both modes)
            k=k+1;
        end
    end

    %% Even uneven ordering
    p_e_u=[1:2:n 2:2:n];

    si=size(index,1);
    %Create elementary Hamiltonians
    H=cell(1,si);
    Ui=H;
    H0=sparse(2^(2*s),2^(2*s));
    Hij=H0;
    for i=1:si
        Hi=H0;
        in=index(i,:);
        Hi(in(1),in(2))=prefactor(i);
        Hi(in(2),in(1))=prefactor(i);
        Hij=Hij+Hi;
        H{i}=Hi;
    end
    
    Uij=expm(-1i*Hij*t);
    indexes=Gray_Indexes_Const_N(n);
    if trott_type==0
        for i=1:si
            Ui{i}=expm(-1i*H{i}*t);
        end

        %% Even uneven ordering
        U=Ui{p_e_u(1)};
        for j=2:si
            U=U*Ui{p_e_u(j)};
        end
        [~,F_e_u]=Average_Fidelity(U,Uij,indexes);
    elseif trott_type==1
        for i=1:si
            Ui2{i}=expm(-1i*H{i}*t/2);
        end
        U1=Ui2{p_e_u(1)};
        U2=U1;
        for j=2:si
            U1=U1*Ui2{p_e_u(j)};
            U2=Ui2{p_e_u(j)}*U2;
        end
        U=(U1*U2);
        [~,F_e_u]=Average_Fidelity(U,Uij,indexes);
    else
        U=U2MWG_Higher_Trotter(H(p_e_u),t,1,rec_order,2);
        [~,F_e_u]=Average_Fidelity(U,Uij,indexes);
    end
    
        %Calculate all permutations
    %p=perms(1:si);
    p=multi_randperm(si,min([factorial(si),how_many]));
    pi=size(p,1);
    F=zeros(min([factorial(si),how_many]),1);
    indexes=Gray_Indexes_Const_N(n);
    if trott_type==0
        for i=1:si
            Ui{i}=expm(-1i*H{i}*t);
        end

        %% Even uneven ordering
        for i=1:min([factorial(si),how_many])
            U=Ui{p(i,1)};
            for j=2:si
                U=U*Ui{p(i,j)};
            end
            [~,F(i)]=Average_Fidelity(U,Uij,indexes);
        end
    elseif trott_type==1
        for i=1:si
            Ui2{i}=expm(-1i*H{i}*t/2);
        end
        
        for i=1:min([factorial(si),how_many])
             U1=Ui2{p(i,1)};
            U2=U1;
            for j=2:si
                U1=U1*Ui2{p(i,j)};
                U2=Ui2{p(i,j)}*U2;
            end
            U=(U1*U2);
            [~,F(i)]=Average_Fidelity(U,Uij,indexes);
        end
    else
        for i=1:min([factorial(si),how_many])
            U=U2MWG_Higher_Trotter(H(p(i,:)),t,1,rec_order,2);
            [~,F(i)]=Average_Fidelity(U,Uij,indexes);
        end
    end
    
    
    if plotting==1
        title_str=['$N_\text{max}=' num2str(n) '$'];
%         if n~=n_max
%             pretty_histogram(F,F_e_u,[],n_bins,0,title_str);
%         else
%             pretty_histogram(F,F_e_u,[],n_bins,1,title_str);
%         end
        if n~=n_max
            [n_F_u F_u]=pretty_histogram(F,F_e_u,[],n_bins,[0 1],title_str);
        else
            [n_F_u F_u]=pretty_histogram(F,F_e_u,[],n_bins,[0 1],title_str);
        end
        %matlab2tikz(['Histograms\pretty_histogram_const_n_' num2str(n) '_tt_' num2str(trott_type) '_r_' num2str(rec_order,'%i') '.tex'],'parseStrings',false,'width','6.2cm','height','4.5cm')
        drawnow()
    end
    
%     fprintf(['n=' num2str(n) ' ------------------\n'])
%     [F,order]=sort(F);
%     p_order=p(order,:);
%     F_u=uniquetol(F,10^-14);
%     for i=1:length(F_u)
%         fprintf(['  F_n=' num2str(F_u(i)) ' -----\n'])
%         inds=find(abs(F-F_u(i))<10^-15);
%         p_order(inds,:)
%     end
    
        %Performance of Even-Uneven ordering
    fprintf(['N=' num2str(n) ' -------------------------------------\n'])
    F_e_u;
    len_F(n)=length(F);
    Feu(n)=F_e_u
    mean_F(n)=mean(F)
    len_better=sum(F>F_e_u);
    bet(n)=len_better/len_F(n)
    len_bins(n)=length(uniquetol(F,10^-16*9))
end
sum(bet(2:7).*(2:7))/sum(2:7)

Feu(1)=1;
E_eu=1-Feu;
figure()
plot(E_eu)