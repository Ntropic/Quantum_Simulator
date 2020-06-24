%all_number_Pauli_Ladder_Permutator_comparison.m
%Permutates the exchange operations for different total numbers of
%particles
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
mode=1;

figure('Position',[100 100 800 400])
for n_n=2:n_max
    s=ceil(log2(n_n+1));

    %Create Hamiltonian and seperate into elementary interaction parts
    tree=bintree(s,1);
    gray=all_num_algorithm(tree);


    %% Create Gray ordering
    i=1;
    prefactor2=[];
    for num1=1:n_n+1 %This generates lists for non symmetric mode!
        for num2=1:n_n+2-num1
            gray_num=gray([num1,num2]);
            graybin=dec2bin(gray_num-1,s)';
            graybin=graybin(:)';
            graybin=graybin-'0';
            if num1>1 && num2<2^n_n
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
    indexes=Gray_Indexes(n_n);
    [~,F_g(n_n-1)]=Average_Fidelity(U,Uij,indexes);

    %% Simple Ordering
    k=1;
    prefactor=[];
    index=[];

    for n=1:n_n
        for i=1:n
            j=n-i+1;
            if mode==0 
                prefactor=[prefactor;1+0*(i*j)];  %Strength of interaction 
            else
                prefactor=[prefactor;sqrt(i*j)];
            end
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
    for l=1:n_n
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
    indexes=Gray_Indexes(n);
    p=multi_randperm(si,min([factorial(si),how_many]));
    pi=size(p,1);
    F=zeros(1,min([factorial(si),how_many]));
    if min([factorial(si),how_many])>5000
        fprintf([' - > 0/' num2str(min([factorial(si),how_many])) '\n']);
        timer=tic();
        for i=1:min([factorial(si),how_many])
            U=Ui{p(i,1)};
            for j=2:si
                U=U*Ui{p(i,j)};
            end
            [~,F(i)]=Average_Fidelity(U,Uij,indexes);
            if mod(i,500)==0
                fprintf([' - > ' num2str(i) '/' num2str(min([factorial(si),how_many])) ' - t=' num2str(toc(timer)) '\n']);
            end
        end
    else
        for i=1:min([factorial(si),how_many])
            U=Ui{p(i,1)};
            for j=2:si
                U=U*Ui{p(i,j)};
            end
            [~,F(i)]=Average_Fidelity(U,Uij,indexes);
        end
    end

    %% Even uneven ordering
    U=Ui{p_e_u(1)};
    for j=2:si
        U=U*Ui{p_e_u(j)};
    end
    [~,F_e_u(n_n-1)]=Average_Fidelity(U,Uij,indexes);
    %F_e_u=Quick_Average_Fidelity(Uij,U);

    h=histogram(F,n_bins,'BinLimits',[0 1],'EdgeColor','none','Normalization','probability');
    hi=h.Values;
    binE=h.BinEdges;
    binE=binE(1:end-1)+(binE(2)-binE(1))/2;
    A(n_n-1,:)=hi/sum(hi);

%         [F,order]=sort(F);
%         p_order=p(order,:);
%         inds=find(abs(F-F(end))<10^-15);
%         if mode==0
%             n_F(n)=n_F_u;
%             minmax_F(n,:)=[min(F_u) max(F_u)];
%             n_F_best(n)=length(inds);
%         else
%             n_F2(n)=n_F_u;
%             minmax_F2(n,:)=[min(F_u) max(F_u)];
%             n_F_best2(n)=length(inds);
%         end
end
imagesc(2:n_n,binE,A');
colormap(parula(256))
set(gca,'YDir','normal')
set(gca,'XTick',[2:n]);
xlabel('$N_\text{max}$');
hold on;
xi=[];
yi=[];
for i=1:n-1
    xi=[xi i+0.5 i+1.5];
    yi=[yi [1 1]*F_e_u(i)];
end
plot(xi,yi,'r--')
xi2=[];
yi2=[];
for i=1:n-1
    xi2=[xi2 i+0.5 i+1.5];
    yi2=[yi2 [1 1]*F_g(i)];
end
plot(xi2,yi2,'b--')

cbar=colorbar();
ylabel(cbar,'$p(F_1,N_\text{max})$');
drawnow()

matlab2tikz('Histograms/overview_all_n.tex','standalone',false,'parseStrings',false,'width','12cm','height','4.5cm');

% fprintf('Ising\n')
% n_F  % Ising -> seems to be A001405 	n_F(n) = binomial(n, floor(n/2)). -> tested till n=8
% minmax_F
% n_F_best
% fprintf('Boson\n')
% n_F2 
% minmax_F2
% n_F_best2