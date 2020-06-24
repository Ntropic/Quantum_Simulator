%const_number_Pauli_Ladder_Permutator.m
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
how_many=100000; %The first how many permutations (maximum)
n_bins=100;
t=pi/4;
plotting=1;
printing=0;
trott_type=0;

for mode=1 %Constant and bosonic strength of interaction
    for n=2:n_max
        s=ceil(log2(n+1));

        %Create Hamiltonian and seperate into elementary interaction parts
        tree=bintree(s,1);
        gray=all_num_algorithm(tree);

        k=1;
        prefactor=[];
        index=[];
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

        %% Even uneven ordering
        p_e_u=[];
        for i=1:2:n
            p_e_u=[p_e_u i];
        end
        for i=2:2:n
            p_e_u=[p_e_u i];
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
        indexes=Gray_Indexes_Const_N(n);
        p=multi_randperm(si,min([factorial(si),how_many]));
        pi=size(p,1);
        F=zeros(1,min([factorial(si),how_many]));
        if min([factorial(si),how_many])>5000
            if printing==1
                fprintf([' - > 0/' num2str(min([factorial(si),how_many])) '\n']);
            end
            timer=tic();
            for i=1:min([factorial(si),how_many])
                U=Ui{p(i,1)};
                for j=2:si
                    U=U*Ui{p(i,j)};
                end
                [~,F(i)]=Average_Fidelity(U,Uij,indexes);
                if mod(i,500)==0
                    if printing==1
                        fprintf([' - > ' num2str(i) '/' num2str(min([factorial(si),how_many])) ' - t=' num2str(toc(timer)) '\n']);
                    end
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
        indexes=Gray_Indexes_Const_N(n);
        [~,F_e_u]=Average_Fidelity(U,Uij,indexes);
        %F_e_u=Quick_Average_Fidelity(Uij,U);
    
        if plotting==1
            title_str=['$n=' num2str(n) '$'];
            if n~=n_max
                [n_F_u F_u]=pretty_histogram(F,F_e_u,[],n_bins,[0 1],title_str);
            else
                [n_F_u F_u]=pretty_histogram(F,F_e_u,[],n_bins,[0 1],title_str);
            end
            drawnow();
            if mode==0
                matlab2tikz(['Histograms\const_n_pretty_histogram_' num2str(n) '_equal_strength.tex'],'parseStrings',false,'width','6.2cm','height','4.5cm'); %'standalone',true,
            else
                matlab2tikz(['Histograms\const_n_pretty_histogram_' num2str(n) '.tex'],'parseStrings',false,'width','6.2cm','height','4.5cm'); %'standalone',true,
            end

            [F,order]=sort(F);
            p_order=p(order,:);
            inds=find(abs(F-F(end))<10^-15);
            if mode==0
                n_F(n)=n_F_u;
                minmax_F(n,:)=[min(F_u) max(F_u)];
                n_F_best(n)=length(inds);
            else
                n_F2(n)=n_F_u;
                minmax_F2(n,:)=[min(F_u) max(F_u)];
                n_F_best2(n)=length(inds);
            end
            n_F_u;
        end
        %Find the blocks
        un_ique=uniquetol(F,10^-15);
        ind_exes=cell(length(un_ique),1);
        for i=1:length(un_ique)
            ind_exes{i}=find(abs(F-un_ique(i))<10^-15);
            len_ind(i)=length(ind_exes{i});
        end
        fprintf(['n=' num2str(n) '  -------------------------------\n'])
        len_ind/2
        p_red=cell(length(un_ique),1);
        for i=1:length(un_ique)
            p_now=p(ind_exes{i},:);
            j=1;
            while j<size(p_now,1)
                p_n=p_now(j,:);
                [tf in]=ismember(p_n(end:-1:1),p_now,'rows');
                if tf==1
                    p_now(in,:)=[];
                end
                j=j+1;
            end
            p_red{i}=p_now;
        end
        for i=1:length(un_ique)
            p_red{i}
        end
    end
end

% if plotting==1
%     fprintf('Ising\n')
%     n_F  % Ising -> seems to be A001405 	n_F(n) = binomial(n, floor(n/2)). -> tested till n=8
%     minmax_F
%     n_F_best
%     fprintf('Boson\n')
%     n_F2 
%     minmax_F2
%     n_F_best2
% end