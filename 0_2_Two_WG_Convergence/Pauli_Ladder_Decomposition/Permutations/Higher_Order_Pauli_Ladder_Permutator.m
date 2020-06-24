%Pauli_Ladder_Permutator.m
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
how_many=10000; %The first how many permutations (maximum)
n_bins=100;
t=pi/4;
plotting=1;

for n=6%1:n_max
    %% Create Gray ordering
    s=ceil(log2(n+1));

    %Create Hamiltonian and seperate into elementary interaction parts
    tree=bintree(s,1);
    gray=all_num_algorithm(tree);

    index=zeros((n+1)*n/2,2);


    i=1;
    prefactor2=[];
    for num1=1:n+1 %This generates lists for non symmetric mode!
        for num2=1:n+2-num1
            gray_num=gray([num1,num2]);
            graybin=dec2bin(gray_num-1,s)';
            graybin=graybin(:)';
            graybin=graybin-'0';
            if num1>1 && num2<2^n
                gray_num2=gray([num1-1,num2+1]);
                graybin2=dec2bin(gray_num2-1,s)';
                graybin2=graybin2(:)';
                graybin2=graybin2-'0';
                prefactor2=[prefactor2;sqrt(num2)*sqrt(num1-1)];  %Strength of interaction 
                index2(i,:)=[bin2dec(num2str(graybin)) bin2dec(num2str(graybin2))]+1;
                i=i+1;
            end
        end
    end

    si=size(index2,1);
    %Create elementary Hamiltonians
    H2=cell(1,si);
    Ui2=H2;
    H0=sparse(2^(2*s),2^(2*s));
    Hij=H0;
    for i=1:si
        Hi2=H0;
        in=index2(i,:);
        Hi2(in(1),in(2))=prefactor2(i);
        Hi2(in(2),in(1))=prefactor2(i);
        Hij=Hij+Hi2;
        H2{i}=Hi2;
    end
    Uij=expm(-1i*Hij*t);
    for i=1:si
        Ui2{i}=expm(-1i*H2{i}*t);
    end
    U=Ui2{1};
    for j=2:si
        U=U*Ui2{j};
    end
    indexes=Gray_Indexes(n);
    [~,F_g]=Average_Fidelity(U,Uij,indexes);
    %F_g=Quick_Average_Fidelity(Uij,U);

    %% Simple Ordering
    k=1;
    prefactor=[];
    for l=1:n
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
    p_e_u=[];
    for l=1:n
        pl=sum(1:l-1);
        for i=1:2:l
            p_e_u=[p_e_u pl+i];
        end
        for i=2:2:l
            p_e_u=[p_e_u pl+i];
        end
    end

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
    for i=1:si
        Ui{i}=expm(-1i*H{i}*t);
    end

    %Calculate all permutations
    %p=perms(1:si);
    p=multi_randperm(si,min([factorial(si),how_many]));
    pi=size(p,1);
    F=zeros(1,min([factorial(si),how_many]));
    indexes=Gray_Indexes(n);
    for i=1:min([factorial(si),how_many])
        U=Ui{p(i,1)};
        for j=2:si
            U=U*Ui{p(i,j)};
        end

        [~,F(i)]=Average_Fidelity(U,Uij,indexes);
        %F(i)=Quick_Average_Fidelity(Uij,U);
    end

    %% Even uneven ordering
    U=Ui{p_e_u(1)};
    for j=2:si
        U=U*Ui{p_e_u(j)};
    end
    indexes=Gray_Indexes(n);
    [~,F_e_u]=Average_Fidelity(U,Uij,indexes);
    %F_e_u=Quick_Average_Fidelity(Uij,U);

    if plotting==1
        title_str=['$N_\text{max}=' num2str(n) '$'];
        if n~=n_max
            pretty_histogram(F,F_e_u,[],n_bins,0,title_str);
        else
            pretty_histogram(F,F_e_u,[],n_bins,1,title_str);
        end
        matlab2tikz(['Histograms\pretty_histogram_' num2str(n) '.tex'],'parseStrings',false,'width','6.2cm','height','4.5cm')
        drawnow()
    end

    %% Observe ordering correlations
    % [F,order]=sort(F);
    % n_order=n_i(p);
    % n_order=n_order(order,:);
    % m_order=m_i(p);
    % m_order=m_order(order,:);
    % for i=1:size(n_order,1)
    %     n_now=n_order(i,:);
    %     m_now=m_order(i,:);
    %     [al,be]=sort(n_now);
    %     n_order(i,:)=al;
    %     m_order(i,:)=m_now(be);
    % end
    % %a=[n_order F' p_order];
    % a=[n_order F' m_order];
    % array2table(a)
end