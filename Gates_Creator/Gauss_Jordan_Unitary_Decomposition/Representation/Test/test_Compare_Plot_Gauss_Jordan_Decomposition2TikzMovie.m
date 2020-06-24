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

file_name=['compare_norm_gj_' num2str(n) ];
file_name_fast=['compare_fast_norm_gj_' num2str(n) ];
%Create/Load Controlled Unitaries for different sizes

H=Gray_Exchange_Hamiltonian_Particles(n);
matrix=expm(-1i*H*t);

%% Fast
mat2=Plot_Fast_Gauss_Jordan_Decomposition(matrix,1); 
mat=Plot_Gauss_Jordan_Decomposition(matrix,1); 

h=figure('Position',[100 100 1200 800]);
ax1=subplot(1,1,1);
if make_mov==1
    v=VideoWriter([file_name '.avi'],'Uncompressed AVI');
    v.FrameRate=5;
    open(v);
end
for i=1:length(mat)
    [bra_str ket_str indexes]=Add_PColorMat(ax1,full(mat{i}) ,{'fock_gray_tex_phase',n,2});
    title(['$N_{max}=' num2str(n) '$, Step \#' num2str(i)])
    
    A{i}=PColor2Tikz2Pdf2Png(bra_str,ket_str,indexes,'Aux_Files','frame',800);
    if i>1
        if any(size(A{i-1})~=size(A{i}))
            C=A{i-1};
            B=A{i};
            C(1:min([size(A{i-1},1),size(A{i},1)]),1:min([size(A{i-1},2),size(A{i},2)]),1:min([size(A{i-1},3),size(A{i},3)]))=B(1:min([size(A{i-1},1),size(A{i},1)]),1:min([size(A{i-1},2),size(A{i},2)]),1:min([size(A{i-1},3),size(A{i},3)]));
            A{i}=C;
        end
    end
    if make_mov==1
        writeVideo(v,A{i});
    end
end
for i=1:10
    writeVideo(v,A{end});
end
for i=length(mat):-1:1
    writeVideo(v,A{i});
end
close(v);



if make_mov==1
    v=VideoWriter([file_name_fast '.avi'],'Uncompressed AVI');
    v.FrameRate=5;
    open(v);
end
for i=1:length(mat2)
    [bra_str ket_str indexes]=Add_PColorMat(ax1,full(mat2{i}) ,{'fock_gray_tex_phase',n,2});
    title(['$N_{max}=' num2str(n) '$, Step \#' num2str(i)])
    
    A{i}=PColor2Tikz2Pdf2Png(bra_str,ket_str,indexes,'Aux_Files','frame',800);
    
    if i>1
        if any(size(A{i-1})~=size(A{i}))
            C=A{i-1};
            B=A{i};
            C(1:min([size(A{i-1},1),size(A{i},1)]),1:min([size(A{i-1},2),size(A{i},2)]),1:min([size(A{i-1},3),size(A{i},3)]))=B(1:min([size(A{i-1},1),size(A{i},1)]),1:min([size(A{i-1},2),size(A{i},2)]),1:min([size(A{i-1},3),size(A{i},3)]));
            A{i}=C;
        end
    end
    if make_mov==1
        writeVideo(v,A{i});
    end
end
for i=1:10
    writeVideo(v,A{length(mat2)});
end
for i=length(mat2):-1:1
    writeVideo(v,A{i});
end
close(v);
