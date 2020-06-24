%test_Plot_Gauss_Jordan_Decomposition2TikzMovie.m
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
n=3;
s=ceil(log2(n+1));
t=pi/4;

file_name=['symm_compare_both_gj_' num2str(n) ];
%Create/Load Controlled Unitaries for different sizes

H=Gray_Exchange_Hamiltonian_Particles(n);
matrix=expm(-1i*H*t);

%% Fast
mat2=Plot_Fast_Gauss_Jordan_Decomposition(matrix,1); 
mat=Plot_Gauss_Jordan_Decomposition(matrix,1); 

h=figure('Position',[100 100 600 1000]);
ax1=subplot(2,1,1);
ax2=subplot(2,1,2);
if make_mov==1
    v=VideoWriter([file_name '.avi'],'Uncompressed AVI');
    v.FrameRate=5;
    open(v);
end
for i=1:max([length(mat) length(mat2)])
    if i<=length(mat)
        [bra_str ket_str indexes]=Add_PColorMat(ax1,full(mat{i}) ,{'fock_gray_tex_phase',n,2});
        title(ax1,['$N_{max}=' num2str(n) '$'])
    end
    if i<=length(mat2)
        [bra_str ket_str indexes]=Add_PColorMat(ax2,full(mat2{i}) ,{'fock_gray_tex_phase',n,2});
    end
    drawnow();
    B=PColor2Tikz2Pdf2Png(bra_str,ket_str,indexes,'Aux_Files','frame',1200);
    C(:,:,1)=[B(1:600,:,1);B(700:end,:,1)];
    C(:,:,2)=[B(1:600,:,2);B(700:end,:,2)];
    C(:,:,3)=[B(1:600,:,3);B(700:end,:,3)];
    A{i}=C;
    clear C;
    if i>1
        if any(size(A{i-1})~=size(A{i}))
            C=A{i-1};
            B=A{i};
            C(1:min([size(A{i-1},1),size(A{i},1)]),1:min([size(A{i-1},2),size(A{i},2)]),1:min([size(A{i-1},3),size(A{i},3)]))=B(1:min([size(A{i-1},1),size(A{i},1)]),1:min([size(A{i-1},2),size(A{i},2)]),1:min([size(A{i-1},3),size(A{i},3)]));
            A{i}=C;
            clear C;
        end
    end
    if make_mov==1
        writeVideo(v,A{i});
    end
end
for i=1:10
    writeVideo(v,A{end});
end
c=length(A)+1;
for i=1:max([length(mat) length(mat2)])
    if i<=length(mat)
        [bra_str ket_str indexes]=Add_PColorMat(ax1,full(mat{end-i+1}) ,{'fock_gray_tex_phase',n,2});
        title(ax1,['$N_{max}=' num2str(n) '$'])
    end
    if i<=length(mat2)
        [bra_str ket_str indexes]=Add_PColorMat(ax2,full(mat2{end-i+1}) ,{'fock_gray_tex_phase',n,2});
    end
    drawnow();
    B=PColor2Tikz2Pdf2Png(bra_str,ket_str,indexes,'Aux_Files','frame',1200);
    C(:,:,1)=[B(1:600,:,1);B(700:end,:,1)];
    C(:,:,2)=[B(1:600,:,2);B(700:end,:,2)];
    C(:,:,3)=[B(1:600,:,3);B(700:end,:,3)];
    A{c}=C;
    clear C;
        if any(size(A{c-1})~=size(A{c}))
            C=A{c-1};
            B=A{c};
            C(1:min([size(A{c-1},1),size(A{c},1)]),1:min([size(A{c-1},2),size(A{c},2)]),1:min([size(A{c-1},3),size(A{c},3)]))=B(1:min([size(A{c-1},1),size(A{c},1)]),1:min([size(A{c-1},2),size(A{c},2)]),1:min([size(A{c-1},3),size(A{c},3)]));
            A{c}=C;
            clear C;
        end
    if make_mov==1
        writeVideo(v,A{c});
    end
    c=c+1;
end
for i=1:10
    writeVideo(v,A{end});
end
close(v);