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
%Load gates
load('../../../Gates_Table/elem_gates.mat','-mat')
load('../../../Gates_Table/comp_gates.mat','-mat')
load('../XY_quantum_gates.mat','-mat')

mode=1;

n=6;
N=n-1;%Photon number

%Trotter step number
steps=10;

%Iteration time
t=pi/4;

bins=200;
how_many=10^5;

Hij=XY_Exchange_Hamiltonian(n);
U_exact=expm(-1i*t*Hij);

%Subspace for analysis of algorithm
index=2.^(0:n-1)+1;

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
for i=1:2:2*floor(n/2)-1
    gates_indexes={gates_indexes{:},[i,i+1]};
    param_steps={param_steps{:},[{sym('phi'),sym('phi')*Jx(i)}]};
end
[gates_names{1:length(gates_indexes)}]=deal('XX_YY');
A=Create_Comp_Gate(nameA,sizle,anc_size,gates_names,gates_indexes,param_steps);
A=Generate_Gate_Circuit(A,circ_string,1:sizle,[],[]);

%Gate B
sizle=floor((n-1)/2)*2;
anc_size=0;
nameB=['B_n_' num2str(n)];
circ_string=['B_{XY}^{n=' num2str(n) '}'];
gates_indexes={};
param_steps={};
gates_names={};
for i=1:2:2*floor((n-1)/2)-1
    gates_indexes={gates_indexes{:},[i,i+1]};
    param_steps={param_steps{:},[{sym('phi'),sym('phi')*Jx(i+1)}]};
end
[gates_names{1:length(gates_indexes)}]=deal('XX_YY');
B=Create_Comp_Gate(nameB,sizle,anc_size,gates_names,gates_indexes,param_steps);
B=Generate_Gate_Circuit(B,circ_string,1:sizle,[],[]);
    
    
for l=1:steps
    fprintf([num2str(l) ' Trotter steps\n'])
    timer=tic();

    %Create Trotter steps
    name=['XY_n_' num2str(n) '_sym'];
    circ_string=['XY_{\text{approx, sym}}^{n=' num2str(n) '}'];
    circ_string_step=['XY_{\text{step, sym}}^{n=' num2str(n) '}'];
    gates_names={nameA,nameB};
    gates_indexes={[1:A.size],1+[1:B.size]};
    other_gates1=Comp_Gate_Merger(comp_gates,xy_gates,A,B);
    %This command creates a Trotter gate for all the steps -> Great for
    %creating a gate circuit
    [gate,step_gate,trotter_name,trotter_step_name]=Strang_Trotter(t,l,name,circ_string,circ_string_step,n,0,gates_names,gates_indexes);
    if l==1
        Gates2Circ2Preview(gate,elem_gates,other_gates1,'xy_gate','gate=expansion',2,30);
    end

    %Do Trotter steps
    other_gates=Comp_Gate_Merger(other_gates1,gate,step_gate);
    other_gates=Gates_Tables_Prepare(other_gates,elem_gates);
    ind=Gate_Index_by_Name(trotter_name,elem_gates,other_gates);
    U_approx=Expansion2Matrix(other_gates(ind(2)));
    

    %average_fidelity
    [F_min(l) F_av(l)]=Min_Av_Fidelity(U_approx,U_exact,index);
    U_ap{l}=U_approx;
    %average_fidelity
    time_l(l)=toc(timer);
end 

bins_int=[min(1-F_av) max(1-F_min)];
bins_int=log10([bins_int(2),bins_int(1)]);
xs=10.^(linspace(bins_int(1),bins_int(2),bins));
x=xs;

phi_noon=NOON_Wave_XY(N);
phi=RandWave(n,how_many);
for l=1:steps
    F_bins(l,:)=Fidelity(U_ap{l},U_exact,index,phi);
    F_noon(l)=Fidelity(U_ap{l},U_exact,index,phi_noon(index));
    F_hist(l,:)=hist(1-F_bins(l,:),xs(end:-1:1));
end
F_hist=F_hist(:,end:-1:1);


%Calculate Order of Convergence
f1_bins=1-F_bins;
for i=1:how_many
    order(:,i)=Order_Of_Convergence(f1_bins(:,i));
end
xo_min=min(order(:));
xo_max=max(order(:));
bins_o=[xo_min xo_max];
xo=linspace(bins_o(1),bins_o(2),bins);
for l=1:steps-1
    order_hist(l,:)=hist(order(l,:),xo);
end

f1_av=1-F_av;
f1_min=1-F_min;
f1_bins=F_hist'/max(F_hist(:));
order_hist=order_hist'/max(order_hist(:));
f1_noon=1-F_noon;

Order_f1_av=Order_Of_Convergence(f1_av);
Order_f1_min=Order_Of_Convergence(f1_min);
Order_f1_noon=Order_Of_Convergence(f1_noon);

%% 1-Fidelity Plot
figure('Position',[100 100 1200 800])
ax1=subplot(1,2,1);
lin_x=linspace(min(x),max(x),length(x));
dx=lin_x(2)-lin_x(1);
lin_x=lin_x+dx/2-(x(end-1)-x(end));
imagesc(ax1,1:steps,lin_x,f1_bins(end:-1:1,:),'HandleVisibility','off');
je=gray(255);
je=je(end:-1:1,:);
colormap(je)
caxis([0 min(max(f1_bins))*4])
set(ax1,'YDir','normal')
set(ax1, 'YScale', 'log')
hold on;
plot3(ax1,1:steps,f1_av,ones(size(f1_av)),'k')
plot3(ax1,1:steps,f1_min,ones(size(f1_min)),'--k')
plot3(ax1,1:steps,f1_noon,ones(size(f1_noon)),'--r')
xlabel(ax1,'Symmetric Trotter steps')
ylabel(ax1,'1-Fidelity')
legend(ax1,{'Average Fidelity','Minimum Fidelity','NOON Fidelity'},'Location','northeast')
axis(ax1,[1 steps min(x) max(x) 0 1]);
view(ax1,0,90)
hold off;
drawnow;

%% Order of Convergence of 1-Fidelity Plot
ax2=subplot(1,2,2);
lin_xo=linspace(min(xo),max(xo),length(xo));
dxo=lin_xo(2)-lin_xo(1);
lin_xo=lin_xo-dxo/2;
imagesc(ax2,1:(steps-1),lin_xo(end:-1:1),order_hist(end:-1:1,:),'HandleVisibility','off');
je=gray(255);
je=je(end:-1:1,:);
colormap(je)
caxis([0 max(max(order_hist))*2])
set(ax2,'YDir','normal')
hold on;
plot3(ax2,1:(steps-1),Order_f1_av,ones(size(Order_f1_av)),'k')
plot3(ax2,1:(steps-1),Order_f1_min,ones(size(Order_f1_min)),'--k')
plot3(ax2,1:(steps-1),Order_f1_noon,ones(size(Order_f1_noon)),'--r')
xlabel(ax2,'Symmetric Trotter steps')
ylabel(ax2,'1-Fidelity')
legend(ax2,{'Average Fidelity','Minimum Fidelity','NOON Fidelity'},'Location','northeast')
axis(ax2,[1 steps-1 min(xo) max(xo) 0 1]);
view(ax2,0,90)
hold off;
drawnow;

%matlab2tikz(['Fidelity_Density/sym_fidelity_density_N_' num2str(N) '.tex'],'standalone',true)