%Trotter_WG_Plot_Test_Compare_Depth_Final_Compare.m
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

%Waveguides
n=6;
N=n-1;  %Photons? classically

depth=2;
s_e=1;
d_e=1;
how_often=50;
horz_lines=1;

%Trotter Iteration number
max_steps=10;

%Iteration time
t=pi/4;

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
  
%% Plot exact evolution
f=figure('Position', [10 10 1300 400]);

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
    
ax1=subplot(1,1,1);
[surf_mat xi ti dx]=WG_Plot(phi_i_red,ts,1,10,1);
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
a.FaceAlpha = 0.25;
a.EdgeColor='None';
a=patch([x_lim x_lim(end:-1:1)],[[1 1]*(i+1+dx) [1 1]*(i+1.5)], 'k');
a.FaceAlpha = 0.25;
a.EdgeColor='None';

ax1.YTick=[1:length(ket_string2)];
ax1.YTickLabel=ket_string2;
drawnow;
matlab2tikz(['WG/exact_n_' num2str(n) '_t_' num2str(t) '.tex'],'standalone',true,'extraCode','\usepackage{braket}','parseStrings',false,'width','6cm','height','4cm')

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
    ax2=subplot(1,1,1);
    [surf_mat xi ti dx]=WG_Plot(phi_a_red,ta,2,how_often,1); %2 smooth interpolation
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
    a.FaceAlpha = 0.25;
    a.EdgeColor='None';
    a=patch([x_lim x_lim(end:-1:1)],[[1 1]*(i+1+dx) [1 1]*(i+1.5)], 'k');
    a.FaceAlpha = 0.25;
    a.EdgeColor='None';

    if horz_lines==1
        y_lim=ylim(gca);
        for i=2:length(ta)-1
            plot([1 1]*ta(i),ylim,'Color',[1 1 1]*0.5);
        end
    end
    if length(ket_string2)==0
        ax2.YTick=[];
    else
        ax2.YTick=[1:length(ket_string2)];
        ax2.YTickLabel=ket_string2;
    end
    %cb=colorbar();
    %ylabel(cb,'p')
    set(gcf,'color',0.98*[1 1 1])
    drawnow;
    
    matlab2tikz(['WG/trotter_n_' num2str(n) '_s_' num2str(steps)  '_t_' num2str(t) '.tex'],'standalone',true,'extraCode','\usepackage{braket}','parseStrings',false,'width','6cm','height','4cm')
end