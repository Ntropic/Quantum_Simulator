%Total_Increase_All_GateNumbers2Datasets.m
%How many more Trotter steps are needed when adding photons  or modes
%(compared to 2 photons and one mode)
clear all;
close all;
clc;

p=pwd;
if any(strfind(p,'\'));
    elem=strsplit(p,'\');
else
    elem=strsplit(p,'/');
end
shortened=fullfile(elem{1:end-3});
addpath(genpath(shortened));

%How many Trotter
st=2;%2:10;
%Which Type of Simulation?
method2=2:3; %Gray?

%Which algorithm?
time_step2=[1 2 4 5]; %Trotter?


for j=1:2
    fprintf([num2str(j) 'th dataset: \n'])
    if j==1
        %Photon number
        n_m=1:7;

        %Trotter Iteration number
        steps=10;

        N_m=2:4;
    elseif j==2
        %Photon number
        n_m=1:3;

        %Trotter Iteration number
        steps=10;

        N_m=4:6;
    end

    r_h=0:3;

    %Iteration time
    t=pi/4;

    %% Fidelities and gate numbers (gn)
    %% ------------------------------------------------------------------------------------------------------------------------
    exits={'_pauli_ladder','_pauli_gray','_gj'};
    exit_Folders={'Pauli_Ladder_Decomposition','Pauli_Gray_Decomposition','Gauss_Jordan_Decomposition'};

    trot_name=['mwg_trot_compare_fid_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_N_' num2str(max(N_m))];
    trot_name_gn=['mwg_trot_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_N_' num2str(max(N_m))];

    sym_name=['mwg_sym_compare_fid_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_N_' num2str(max(N_m))];
    sym_name_gn=['mwg_sym_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_N_' num2str(max(N_m))];

    fprintf('-> Check if simulation results are located in the designated folders:\n\n')

    %% -------------------------------------------------------------------------------------------------------------------------
    for i=1:length(exits) %Only search for gate numbers
        fprintf(['- Searching for calculations of type ' exit_Folders{i} '\n'])

        if exist(['../' exit_Folders{i} '/Simple_Trotter/' trot_name exits{i} '.mat'])
            a=load(['../' exit_Folders{i} '/Simple_Trotter/' trot_name exits{i} '.mat']);
            F_trot=a.F_av_ind;
            fprintf(' - Found Simple Trotter calculations\n')
        else
            error('Not found (trott fid)')
        end
        if exist(['../' exit_Folders{i} '/Simple_Trotter/' trot_name_gn exits{i} '.mat'])
            a=load(['../' exit_Folders{i} '/Simple_Trotter/' trot_name_gn exits{i} '.mat']);
            s_q_trot=a.s_q;
            fprintf(' - Found Simple Trotter (gate number) calculations\n')
        else
            error('Not found (trott gate numbers)')
        end

        %% -------------------------------------------------------------------------------------------------------------------------
        if exist(['../' exit_Folders{i} '/Symmetric_Trotter/' sym_name exits{i} '.mat'])
            a=load(['../' exit_Folders{i} '/Symmetric_Trotter/' sym_name exits{i} '.mat']);
            F_sym=a.F_av_ind;
            fprintf(' - Found Simple Trotter calculations\n')
        else
            error('Not found (sym fid)')
        end
        if exist(['../' exit_Folders{i} '/Symmetric_Trotter/' sym_name_gn exits{i} '.mat'])
            a=load(['../' exit_Folders{i} '/Symmetric_Trotter/' sym_name_gn exits{i} '.mat']);
            s_q_sym=a.s_q;
            fprintf(' - Found Symetric Trotter (gate number) calculations\n')
        else
            error('Not found (sym gate numbers)');
        end

        %% -------------------------------------------------------------------------------------------------------------------------
        counter=0;
        for r=r_h
            counter=counter+1;
            higher_name_r0=['mwg_higher_compare_fid_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_r_' num2str(r) '_N_' num2str(max(N_m))];
            higher_name_r0_gn=['mwg_higher_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_r_' num2str(r) '_N_' num2str(max(N_m))];

            if exist(['../' exit_Folders{i} '/Higher_Order_Trotter/' higher_name_r0 exits{i} '.mat'])
                a=load(['../' exit_Folders{i} '/Higher_Order_Trotter/' higher_name_r0 exits{i} '.mat']);
                F_high_r0{counter}=a.F_av_ind;
                fprintf([' - Found Higher Order Trotter [r=' num2str(r) '] calculations\n'])
            else
                error('Not found (higher fid)')
            end
            if exist(['../' exit_Folders{i} '/Higher_Order_Trotter/' higher_name_r0_gn exits{i} '.mat'])
                a=load(['../' exit_Folders{i} '/Higher_Order_Trotter/' higher_name_r0_gn exits{i} '.mat']);
                s_q_high_r0{counter}=a.s_q;
                fprintf([' - Found Higher_Order_Trotter [r=' num2str(r) '] (gate number) calculations\n'])
            else
                error(['Not found (Higher_Order_Trotter) [r=' num2str(r) ']'])
            end
        end
        
        if j==1
            if i==1
                F={};
                s_q={};
            end
            F={F{:},F_trot,F_sym,F_high_r0{:}};
            s_q={s_q{:},s_q_trot,s_q_sym,s_q_high_r0{:}};
            if i==1
                len_sq=length(s_q);
            end
        else
            if i==1
                F2={};
                s_q2={};
            end
            F2={F2{:},F_trot,F_sym,F_high_r0{:}};
            s_q2={s_q2{:},s_q_trot,s_q_sym,s_q_high_r0{:}};
        end
    end
end

F=reshape(F,len_sq,length(exits));
F=F';
s_q=reshape(s_q,len_sq,length(exits));
s_q=s_q';

F2=reshape(F2,len_sq,length(exits));
F2=F2';
s_q2=reshape(s_q2,len_sq,length(exits));
s_q2=s_q2';

fprintf('\n--> All files have been found, all variables extracted. <--\n')
fprintf('         --> !Proceeding with data analysis! <--\n\n')

%% Analyze Data---------------------------------------------------------------------------------------------------------
st2=st;

for method=method2
    for time_step=time_step2
        for st=st2

            %Get data
            F_av_ind=F{method,time_step};
            F_av_ind2=F2{method,time_step};
            F_av_ind=F_av_ind(:,:,1:st);
            F_av_ind2=F_av_ind2(:,:,1:st);
            F_av_ind(1,1)=1;

            % Compute order of convergence for every point
            for i=1:size(F_av_ind,1)
                for j=1:size(F_av_ind,2)
                    ooc=Order_Of_Convergence(1-F_av_ind(i,j,:));
                    c(i,j)=round(ooc(end));
                    e(i,j)=1-F_av_ind(i,j,st);
                end
            end
            for i=1:size(F_av_ind2,1)
                for j=1:size(F_av_ind2,2)
                    ooc=Order_Of_Convergence(1-F_av_ind2(i,j,:));
                    c2(i,j)=round(ooc(end));
                    e2(i,j)=1-F_av_ind2(i,j,st);
                end
            end

            for i=1:size(F_av_ind,1) %Waveguides
                for j=1:size(F_av_ind,2)-1
                    cn=c(i,j);
                    rat(i,j)=(e(i,j+1)/e(i,j))^(1/cn);
                end
            end
            for i=1:size(F_av_ind2,1) %Waveguides
                for j=1:size(F_av_ind2,2)-1
                    cn=c2(i,j);
                    rat2(i,j)=(e2(i,j+1)/e2(i,j))^(1/cn);
                end
            end
            cn=round(mean(c2(:)));

            %% Efficiency
            if method~=3
                eff=zeros(size(rat)+[0 1]);
                for i=1:size(rat,1)
                    if i==1
                        eff(1,1)=NaN;
                        eff(1,2)=1; %Normalization
                        for j=3:size(rat,2)+1
                            eff(i,j)=eff(i,j-1)*rat(i,j-1);
                        end
                    else
                        if i==2 %First calculate the effort between (1,2) and (i,1)
                            eff(2,1)=(e(2,1)/e(1,2))^(1/cn);
                        else    %Calculate effort between (i,1) and (i-1,1)
                            eff(i,1)=(e(i,1)/e(1,2))^(1/cn);
                        end
                        for j=2:size(rat,2)+1
                            eff(i,j)=eff(i,j-1)*rat(i,j-1);
                        end
                    end
                end
                eff=eff/min(eff(:));

                eff2=zeros(size(rat2)+[0 1]);
                for i=1:size(rat2,1)
                    if i==1
                        eff2(1,1)=NaN;
                        eff2(1,2)=1; %Normalization
                        for j=3:size(rat2,2)+1
                            eff2(i,j)=eff2(i,j-1)*rat2(i,j-1);
                        end
                    else
                        if i==2 %First calculate the effort between (1,2) and (i,1)
                            eff2(2,1)=(e2(2,1)/e2(1,2))^(1/cn);
                        else    %Calculate effort between (i,1) and (i-1,1)
                            eff2(i,1)=(e2(i,1)/e2(1,2))^(1/cn);
                        end
                        for j=2:size(rat2,2)+1
                            eff2(i,j)=eff2(i,j-1)*rat2(i,j-1);
                        end
                    end
                end
                eff2=eff2/min(eff2(:));
            else %GJ approach (remove lowest order
                eff=zeros(size(rat)+[0 1]);
                for i=2:size(rat,1)
                    if i==2
                        eff(1,:)=NaN;
                        eff(2,1)=NaN;
                        eff(2,2)=1; %Normalization
                        for j=3:size(rat,2)+1
                            eff(i,j)=eff(i,j-1)*rat(i,j-1);
                        end
                    else
                        if i==3 %First calculate the effort between (1,2) and (i,1)
                            eff(3,1)=(e(3,1)/e(2,2))^(1/cn);
                        else    %Calculate effort between (i,1) and (i-1,1)
                            eff(i,1)=(e(i,1)/e(2,2))^(1/cn);
                        end
                        for j=2:size(rat,2)+1
                            eff(i,j)=eff(i,j-1)*rat(i,j-1);
                        end
                    end
                end
                eff=eff/min(eff(:));

                eff2=zeros(size(rat2)+[0 1]);
                for i=2:size(rat2,1)
                    if i==2
                        eff2(1,:)=NaN;
                        eff2(2,1)=NaN;
                        eff2(2,2)=1; %Normalization
                        for j=3:size(rat2,2)+1
                            eff2(i,j)=eff2(i,j-1)*rat2(i,j-1);
                        end
                    else
                        if i==3 %First calculate the effort between (1,2) and (i,1)
                            eff2(3,1)=(e2(3,1)/e2(2,2))^(1/cn);
                        else    %Calculate effort between (i,1) and (i-1,1)
                            eff2(i,1)=(e2(i,1)/e2(2,2))^(1/cn);
                        end
                        for j=2:size(rat2,2)+1
                            eff2(i,j)=eff2(i,j-1)*rat2(i,j-1);
                        end
                    end
                end
                eff2=eff2/min(eff2(:));
            end

        end
        
        figure()
        %% Fit photon number increase in effort 
        for i=2:size(eff,1)
            eff_now=eff(i,:);
            photon_num=1:size(eff,2);           
            [xData, yData] = prepareCurveData(photon_num,eff_now);
            ft=fittype( 'poly1' );
            [fitresult,gof]=fit(xData,yData,ft);
            
            fun=@(p1,p2,x) p1*x+p2;
            p1(i)=fitresult.p1;
            p2(i)=fitresult.p2;
            bounds=confint(fitresult);
            p1_min_95(i)=bounds(1,1);
            p1_max_95(i)=bounds(2,1);
            p2_min_95(i)=bounds(1,2);
            p2_max_95(i)=bounds(2,2);
            bounds=confint(fitresult,0.68);
            p1_min_68(i)=bounds(1,1);
            p1_max_68(i)=bounds(2,1);
            p2_min_68(i)=bounds(1,2);
            p2_max_68(i)=bounds(2,2);
            rsq(i)=gof.adjrsquare;
            
            phot_num{i}=xData;
            eff_num{i}=yData;
            %Plot fit result
        end
%         for i=2:size(eff,1)
%             patch([phot_num{i}',phot_num{i}(end:-1:1)'],[fun(p1_max_95(i),p2_max_95(i),phot_num{i}'),fun(p1_min_95(i),p2_min_95(i),phot_num{i}(end:-1:1)')],0.9*[1 1 1],'EdgeColor','none')
%             hold on;
%         end
%         for i=2:size(eff,1)
%             patch([phot_num{i}',phot_num{i}(end:-1:1)'],[fun(p1_max_68(i),p2_max_68(i),phot_num{i}'),fun(p1_min_68(i),p2_min_68(i),phot_num{i}(end:-1:1)')],0.8*[1 1 1],'EdgeColor','none')
%         end
        for i=2:size(eff,1)
            mar={'b','r','g'};
            plot(phot_num{i}',fun(p1(i),p2(i),phot_num{i}'),mar{i})
            hold on;
        end
        for i=2:size(eff,1)
            mar2={'bo','ro','go'};
            plot(phot_num{i}',eff_num{i}',mar2{i})
        end
        matlab2tikz(['ph_inc' exits{method} '_' num2str(time_step) '.tex'],'standalone',false,'parsestrings',false,'width','6cm','height','4cm')

    end
end