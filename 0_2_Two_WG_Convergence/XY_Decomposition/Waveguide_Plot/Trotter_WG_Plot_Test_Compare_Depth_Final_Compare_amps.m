%Trotter_WG_Plot_Test_Compare_Depth_Final_Compare_amps.m
%Show final distribution and compare with ideal
clc;
clear all;
close all;

p=pwd;
if any(strfind(p,'\'));
    elem=strsplit(p,'\');
else
    elem=strsplit(p,'/');
end
shortened=fullfile(elem{1:end-3});
addpath(genpath(shortened));
%Load gates
load('../../../Gates_Table/elem_gates.mat','-mat')
load('../../../Gates_Table/comp_gates.mat','-mat')
load('XY_quantum_gates.mat','-mat')

plotting=[1 1];
%Waveguides
n=9;
N=n-1;  %Photons? classically

depth=2;
s_e=1;
d_e=1;
how_often=50;
horz_lines=1;

%Trotter Iteration number
max_steps=5;

%Iteration time
t=pi/4;
posx=[0 t/2 t];
stringx={'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$'};
with_header=0;

%Create initial pure state
m=n;
theta=0;    %relative phase
phi0=NOON_Wave(1,n,m,theta);

%% Exact evolution
[Hij]=XY_Exchange_Hamiltonian(n);
%FockPrint(Hij)
s=50;
phi_i=zeros(length(phi0),s);
phi_i(:,1)=phi0;
ts=linspace(0,t,s);
for i=2:s
    phi_i(:,i)=expm(-1i*Hij*ts(i))*phi0;
end

indexes=2.^(0:n-1)+1;
phi_i_red=phi_i(indexes,:);
  
        mode2={'fock',n-1,2};
        l=length(dec2bin(length(indexes)-1));
        for i=1:length(indexes)
            d=[i-1,length(indexes)-i+1];
            indexes2(i)=bin2dec([dec2bin(d(1),l) dec2bin(d(2),l)]);
        end
        %ket_string2=ket_fock_str(indexes2,mode2);
        ket_string2=ket_fock_tex(indexes2,mode2);
        for i=1:length(ket_string2)
            ketter=ket_string2{i};
            ket_string2{i}=['$' ketter(2:end) '$'];
        end
        [surf_mat xi ti dx surf_mat_real surf_mat_imag]=WG_Plot(phi_i_red,ts,1,10,1);
       	id_r=surf_mat_real(:,2); %Save the real imag and abs values for future reference (difference) plots
        id_i=surf_mat_imag(:,2);
        id_a=surf_mat(:,end);
        
if plotting(1)==1
    %% Plot exact evolution
    f=figure('Position', [10 10 1300 400]);

    ax=subplot(1,10,1:2);
    plot(surf_mat_real(:,1),xi,'b--');
    hold on;
    plot(surf_mat_imag(:,1),xi,'r--');
    plot(surf_mat(:,1),xi,'k');
    plot([0,0],[0.5,length(ket_string2)+0.5],'k')

    x_lim=[-1 1];
    for i=1:length(ket_string2)-1
        a=patch([x_lim x_lim(end:-1:1)],[[1 1]*(i+dx) [1 1]*(i+1-dx)], 'k');
        a.FaceAlpha = 0.25;
        a.EdgeColor='None';
        plot(x_lim,[1 1]*(i-dx),'Color',[1 1 1]*0.5);
        plot(x_lim,[1 1]*(i+dx),'Color',[1 1 1]*0.5);
    end
    a=patch([x_lim x_lim(end:-1:1)],[[1 1]*(0) [1 1]*(1-dx)], 'k');
    a.FaceAlpha = 0.25;
    a.EdgeColor='None';
    a=patch([x_lim x_lim(end:-1:1)],[[1 1]*(i+1+dx) [1 1]*(i+1.5)], 'k');
    a.FaceAlpha = 0.25;
    a.EdgeColor='None';
    axis([-1 1 0.5 length(ket_string2)+0.5])
    ax.YTick=[1:length(ket_string2)];
    ax.YTickLabel=ket_string2;
    ax.XTick=[0 1];
    ax.XTickLabel=[0 1];

    ax1=subplot(1,10,3:8);
    imagesc(ti,xi,surf_mat);
    set(gca,'YDir','normal')
    hold on;
    x_lim=xlim(gca);
    for i=1:length(ket_string2)-1
        a=patch([x_lim x_lim(end:-1:1)],[[1 1]*(i+dx) [1 1]*(i+1-dx)], 'k');
        a.FaceAlpha = 0.25;
        a.EdgeColor='None';
        plot(x_lim,[1 1]*(i-dx),'Color',[1 1 1]*0.5);
        plot(x_lim,[1 1]*(i+dx),'Color',[1 1 1]*0.5);
    end
    a=patch([x_lim x_lim(end:-1:1)],[[1 1]*(0) [1 1]*(1-dx)], 'k');
    a.FaceAlpha = 0.25;
    a.EdgeColor='None';
    a=patch([x_lim x_lim(end:-1:1)],[[1 1]*(i+1+dx) [1 1]*(i+1.5)], 'k');
    a.FaceAlpha = 0.25;
    a.EdgeColor='None';
    ax1.YTick=[];
    ax1.YTickLabel=[];
    ax1.XTick=posx;
    ax1.XTickLabel=stringx;
    %xlabel('t')

    ax2=subplot(1,10,9:10);
    plot(surf_mat_real(:,2),xi,'b--');
    hold on;
    plot(surf_mat_imag(:,2),xi,'r--');
    plot(surf_mat(:,end),xi,'k');
    plot([0,0],[0.5,length(ket_string2)+0.5],'k')

    x_lim=[-1 1];
    for i=1:length(ket_string2)-1
        a=patch([x_lim x_lim(end:-1:1)],[[1 1]*(i+dx) [1 1]*(i+1-dx)], 'k');
        a.FaceAlpha = 0.25;
        a.EdgeColor='None';
        plot(x_lim,[1 1]*(i-dx),'Color',[1 1 1]*0.5);
        plot(x_lim,[1 1]*(i+dx),'Color',[1 1 1]*0.5);
    end
    a=patch([x_lim x_lim(end:-1:1)],[[1 1]*(0) [1 1]*(1-dx)], 'k');
    a.FaceAlpha = 0.25;
    a.EdgeColor='None';
    a=patch([x_lim x_lim(end:-1:1)],[[1 1]*(i+1+dx) [1 1]*(i+1.5)], 'k');
    a.FaceAlpha = 0.25;
    a.EdgeColor='None';
    axis([-1 1 0.5 length(ket_string2)+0.5])
    ax2.YTick=[];
    ax2.YTickLabel=[];
    ax2.XTick=[0 1];
    ax2.XTickLabel=[0 1];

    drawnow()
    if with_header==1
        matlab2tikz(['WG_Trotter/exact_n_' num2str(n) '_t_' num2str(t) '_comp.tex'],'standalone',true,'extraCode','\usepackage{braket}','parseStrings',false,'width','14cm','height','3cm')
    else
        matlab2tikz(['WG_Trotter/exact_n_' num2str(n) '_t_' num2str(t) '_comp.tex'],'parseStrings',false,'width','14cm','height','3cm')
    end
end

if plotting(2)==1
    %% Calculate Trotter steps
    for steps=1:max_steps
        %% Approximate decomposition
        %% Prepare gates
        Ni=N:-1:1;
        in=1:N;
        Jx=sqrt(in).*sqrt(Ni);

        %Gate A
        sizle=floor(n/2)*2;
        anc_size=0;
        nameA=['A_n_' num2str(n)];
        circ_string=['A_{XY}^{n=' num2str(n) '}'];
        gates_indexes={};
        param_steps={};
        gates_names={};
        for k=1:2:2*floor(n/2)-1
            gates_indexes={gates_indexes{:},[k,k+1]};
            param_steps={param_steps{:},[{sym('phi'),sym('phi')*Jx(k)}]};
        end
        [gates_names{1:length(gates_indexes)}]=deal('XX_YY');
        A=Create_Comp_Gate(nameA,sizle,anc_size,gates_names,gates_indexes,param_steps);
        %A=Generate_Gate_Circuit(A,circ_string,1:sizle,[],[]);

        %Gate B
        sizle=floor((n-1)/2)*2;
        anc_size=0;
        nameB=['B_n_' num2str(n)];
        circ_string=['B_{XY}^{n=' num2str(n) '}'];
        gates_indexes={};
        param_steps={};
        gates_names={};
        for k=1:2:2*floor((n-1)/2)-1
            gates_indexes={gates_indexes{:},[k,k+1]};
            param_steps={param_steps{:},[{sym('phi'),sym('phi')*Jx(k+1)}]};
        end
        [gates_names{1:length(gates_indexes)}]=deal('XX_YY');
        B=Create_Comp_Gate(nameB,sizle,anc_size,gates_names,gates_indexes,param_steps);
        %B=Generate_Gate_Circuit(B,circ_string,1:sizle,[],[]);

        %Create Trotter steps
        name=['XY_n_' num2str(n)];
        circ_string=['XY_{\text{approx}}^{n=' num2str(n) '}'];
        circ_string_step=['XY_{\text{step}}^{n=' num2str(n) '}'];
        gates_names={nameA,nameB};
        gates_indexes={[1:A.size],1+[1:B.size]};
        other_gates1=Comp_Gate_Merger(comp_gates,xy_gates,A,B);
        %This command creates a Trotter gate for all the steps -> Great for
        [gate,step_gate trotter_name trotter_step_name]=Trotter(t,steps,name,circ_string,circ_string_step,n,0,gates_names,gates_indexes);

        %Do Trotter steps
        other_gates=Comp_Gate_Merger(other_gates1,gate,step_gate);
        other_gates=Gates_Tables_Prepare(other_gates,elem_gates);

        ind=Gate_Index_by_Name(trotter_name,elem_gates,other_gates);
        [phi_a,ta]=DepthExpansion2Matrix_Return_Wave(phi0,other_gates(ind(2)),depth,s_e,d_e);
        phi_a_red=phi_a(indexes,:);


        %% Plotting Trotter steps
        %Plot waveguide development of ideal waveguide
        f=figure('Position', [10 10 1300 400]);
        ax2=subplot(1,6,1:2);
        [surf_mat xi ti dx surf_mat_real surf_mat_imag]=WG_Plot(phi_a_red,ta,2,how_often,1); %2 smooth interpolation
        %h=pcolor(ti,xi,surf_mat);
        %set(h, 'EdgeColor', 'none');
        imagesc(ti,xi,surf_mat);
        set(gca,'YDir','normal')
        hold on;
        x_lim=xlim(gca);
        for i=1:length(ket_string2)-1
            a=patch([x_lim x_lim(end:-1:1)],[[1 1]*(i+dx) [1 1]*(i+1-dx)], 'k');
            a.FaceAlpha = 0.4;
            a.EdgeColor='None';
            plot(x_lim,[1 1]*(i-dx),'Color',[1 1 1]*0.5);
            plot(x_lim,[1 1]*(i+dx),'Color',[1 1 1]*0.5);
        end
        a=patch([x_lim x_lim(end:-1:1)],[[1 1]*(0) [1 1]*(1-dx)], 'k');
        a.FaceAlpha=0.25;
        a.EdgeColor='None';
        a=patch([x_lim x_lim(end:-1:1)],[[1 1]*(i+1+dx) [1 1]*(i+1.5)], 'k');
        a.FaceAlpha=0.25;
        a.EdgeColor='None';

        if horz_lines==1
            y_lim=ylim(gca);
            for i=2:length(ta)-1
                plot([1 1]*ta(i),ylim,'Color',[1 1 1]*0.5);
            end
        end

        ax2.YTick=[];%[1:length(ket_string2)];
        ax2.YTickLabel=[];%ket_string2;
        
        ax2.XTick=[ta(1:2:end)];
        ax2.XTickLabel=[0:length(ta(1:2:end))];
        %xlabel('Steps')
        set(gcf,'color',0.98*[1 1 1])
        drawnow;

        ax2=subplot(1,6,3:4);
        plot(surf_mat_real(:,2),xi,'b--');
        hold on;
        plot(surf_mat_imag(:,2),xi,'r--');
        plot(surf_mat(:,end),xi,'k');
        plot([0,0],[0.5,length(ket_string2)+0.5],'k')

        x_lim=[-1 1];
        for i=1:length(ket_string2)-1
            a=patch([x_lim x_lim(end:-1:1)],[[1 1]*(i+dx) [1 1]*(i+1-dx)], 'k');
            a.FaceAlpha = 0.25;
            a.EdgeColor='None';
            plot(x_lim,[1 1]*(i-dx),'Color',[1 1 1]*0.5);
            plot(x_lim,[1 1]*(i+dx),'Color',[1 1 1]*0.5);
        end
        a=patch([x_lim x_lim(end:-1:1)],[[1 1]*(0) [1 1]*(1-dx)], 'k');
        a.FaceAlpha = 0.25;
        a.EdgeColor='None';
        a=patch([x_lim x_lim(end:-1:1)],[[1 1]*(i+1+dx) [1 1]*(i+1.5)], 'k');
        a.FaceAlpha = 0.25;
        a.EdgeColor='None';
        axis([-1 1 0.5 length(ket_string2)+0.5])
        ax2.YTick=[];
        ax2.YTickLabel=[];
        ax2.XTick=[0 1];
        ax2.XTickLabel=[0 1];

        ax3=subplot(1,6,5:6);
        plot(surf_mat_real(:,2)-id_r,xi,'b--');
        hold on;
        plot(surf_mat_imag(:,2)-id_i,xi,'r--');
        plot(surf_mat(:,end)-id_a,xi,'k');
        plot([0,0],[0.5,length(ket_string2)+0.5],'k')

        x_lim=[-2 2];
        for i=1:length(ket_string2)-1
            a=patch([x_lim x_lim(end:-1:1)],[[1 1]*(i+dx) [1 1]*(i+1-dx)], 'k');
            a.FaceAlpha = 0.25;
            a.EdgeColor='None';
            plot(x_lim,[1 1]*(i-dx),'Color',[1 1 1]*0.5);
            plot(x_lim,[1 1]*(i+dx),'Color',[1 1 1]*0.5);
        end
        a=patch([x_lim x_lim(end:-1:1)],[[1 1]*(0) [1 1]*(1-dx)], 'k');
        a.FaceAlpha = 0.25;
        a.EdgeColor='None';
        a=patch([x_lim x_lim(end:-1:1)],[[1 1]*(i+1+dx) [1 1]*(i+1.5)], 'k');
        a.FaceAlpha = 0.25;
        a.EdgeColor='None';
        axis([-2 2 0.5 length(ket_string2)+0.5])
        ax3.YTick=[];
        ax3.YTickLabel=[];
        ax3.XTick=[0 2];
        ax3.XTickLabel=[0 2];
        
        if with_header==1
            matlab2tikz(['WG_Trotter/trotter_n_' num2str(n) '_s_' num2str(steps)  '_t_' num2str(t) '_comp.tex'],'standalone',true,'extraCode','\usepackage{braket}','parseStrings',false,'width','7cm','height','3cm')
        else
            matlab2tikz(['WG_Trotter/trotter_n_' num2str(n) '_s_' num2str(steps)  '_t_' num2str(t) '_comp.tex'],'parseStrings',false,'width','7cm','height','3cm')
        end
    end
end