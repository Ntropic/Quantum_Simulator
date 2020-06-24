%test_Plot_Gauss_Jordan_Decomposition.m
clear all;
close all;
clc;

make_mov=1;

p=pwd;
if any(strfind(p,'\'));
    elem=strsplit(p,'\');
else
    elem=strsplit(p,'/');
end
shortened=fullfile(elem{1:end-4});
addpath(genpath(shortened));

%Number of qubits per mode
n=7;
s=ceil(log2(n+1));
t=pi/4;
file_name_fast=['fast_gj_' num2str(n) ];
file_name=['gj_' num2str(n) ];
%Create/Load Controlled Unitaries for different sizes

H=Gray_Exchange_Hamiltonian_Particles(n);
matrix=expm(-1i*H*t);

%% Fast
mat2=Plot_Fast_Gauss_Jordan_Decomposition(matrix,1); 

h2=figure();
ax2=subplot(1,1,1);

if make_mov==1
    v=VideoWriter([file_name_fast '.avi'],'Uncompressed AVI');
    v.FrameRate=10;
    open(v);
end
for i=1:length(mat2)
    Add_PColorMat(ax2,full(mat2{i}),{'fock_gray_phase',n,2})
    title(['Step #' num2str(i)]) 
    pause(0.05);
    if make_mov==1
        frame=getframe(h2);
        writeVideo(v,frame);
    end
end
close(v);

%% Default
mat=Plot_Gauss_Jordan_Decomposition(matrix,1); 

h=figure();
ax=subplot(1,1,1);

if make_mov==1
    v=VideoWriter([file_name '.avi'],'Uncompressed AVI');
    v.FrameRate=10;
    open(v);
end
for i=1:length(mat)
    Add_PColorMat(ax,mat{i},{'fock_gray_phase',n,2})
    title(['Step #' num2str(i)]) 
    pause(0.05);
    if make_mov==1
        frame=getframe(h);
        writeVideo(v,frame);
    end
end
close(v);
