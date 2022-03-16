% for each shape, the parameters need to change are: n, True_V,
close all;clear;clc;
addpath ./my_core
%% generate the point cloud in designed shape
radius=1; % the ridus of the sphere
neg_radius=(-1)*radius;
ax=neg_radius;bx=radius;
ay=neg_radius;by=radius;
l1=2*radius;l2=2*radius;
n=10; % partition number
x=linspace(ax,bx,l1*n+1);
y=linspace(ay,by,l2*n+1);
[X,Y]=meshgrid(x,y);
Z=-1.*sqrt(1-X.^2-Y.^2);
Z=real(Z); 
%% calculate the 1st-order volume 
[M,N]=size(Z);
s=l1*l2/((M-1)*(N-1));
V=zeros(M,N);
for i=1:M-1
    for j=1:N-1
        f1=Z(i,j);
        f2=Z(i+1,j);
        f3=Z(i,j+1);
        f4=Z(i+1,j+1);
        average_h=(f1+f2+f3+f4)/4;
        V(i,j)=s*average_h;
    end
end
Volume=sum(sum(V));

%% the second-order volume estimation
% filter 6*6
H=[0 0 4 4 0 0;0 3 -7 -7 3 0;4 -7 3 3 -7 4;4 -7 3 3 -7 4;0 3 -7 -7 3 0;0 0 4 4 0 0];
% figure;mesh(H);
Hf=filter2(H,Z,'same');
SS=l1*l2;  % the projected area
Hf_full=filter2(H,Z,'full');
mf=3; % the size of the filter
[M_full,N_full]=size(Hf_full);
Hf_mf=Hf_full(mf:M_full-mf+1,mf:N_full-mf+1);

dividend=4*12*4;  %four direction, four points, 12 is the common dividend
V_Hf=sum(sum(Hf))*s/dividend;
% for the middle area
mD=ceil(6/2);
V_Hf_middle=sum(sum(Hf(mD:M-mD,mD:N-mD)))*s/dividend;
V_Hf_middle_scale=sum(sum(Hf(mD:M-mD,mD:N-mD)))*M*N*s/((M-2*mD)*(N-2*mD)*dividend);
% in fact, below is the true 2th-order estimation
averageH=sum(sum(Hf))/((M-1)*(N-1));
V_Hf_averageH=averageH*SS/dividend;

nf=6;
Hf_nf=Hf_full(nf:M_full-nf+1,nf:N_full-nf+1);
[M_nf,N_nf]=size(Hf_nf);
averageH_nf=sum(sum(Hf_nf))/(M_nf*N_nf);
V_nf=averageH_nf*SS/dividend;


[M_mf,N_mf]=size(Hf_mf);
averageH_mf=sum(sum(Hf_mf))/(M_mf*N_mf);
V_mf=averageH_mf*SS/dividend;

averageH_full=sum(sum(Hf_full))/(M_full*N_full);
V_full=averageH_full*SS/dividend;

V_all_average2=(V_mf+V_full)/2;

figure;subplot(1,2,1);surf(X,Y,Hf,'EdgeColor','none');
axis tight;axis off;
% axis([-1,1,-1,1]);
subplot(1,2,2);imagesc([-1;1],[-1;1],Hf);
axis tight;axis off;
save('XYHfplot_bathtub_001.mat','X','Y','Hf');



Hf1=reshape(Hf,[],1);
Q3=prctile(Hf1,75);
Q1=prctile(Hf1,25);
Q2=prctile(Hf1,50);
% IQR=Q3-Q1;
% Upper_adjacent=Q3+1.5*IQR;
% Lower_adjacent=Q1-1.5*IQR;
% Hf_a=Hf1(Hf1<Upper_adjacent);
% Hf_a=Hf_a(Hf_a>Lower_adjacent);
% Mean_without_outliers=sum(Hf_a)/(M*N);
% V_Hf_mean=Mean_without_outliers*SS/dividend;

V_Hf_Q2=Q2*SS/dividend;
V_all_average3=(V_Hf_Q2+V_mf+V_full)/3;

%% volume calculation and error estimation
True_V=13*pi/6;
V2=Volume-V_all_average3;
V_2021=Volume-V_Hf_middle;
V_2021_scale=Volume-V_Hf_middle_scale;
V_Hf_averageV=Volume-V_Hf_averageH;
V_nf2=Volume-V_nf;

error1=True_V+Volume;
error2=True_V+Volume-V_all_average3;

%% the second-order estimation
% SS=l1*l2;  % the projected area
% True_V=2*pi/3;
% [V_vector,volume_vector,E_vector,Hf]=second_order(Z,s,SS,True_V,Volume);
% 
% figure;surf(X,Y,Z+1.1,'EdgeColor','none');
% hold on
% imagesc([-1;1],[-1;1],Z+1.1);
% xlabel('X');ylabel('Y');zlabel('Z');
% axis tight;axis off;
%% evaluate the 3-D alpha shape volume and Simpson2d
% [alpha3D_V,error_alpha3D]=alpha_3D_V(X,Y,Z,True_V);
% [volume_simp2,error_simp,coe_equal]=my_simpson2(Z,ax,bx,ay,by,True_V);

%% save the figure and point cloud
% f=getframe(gca);
% imwrite(f.cdata,['./figures/ball' num2str(n) '.tif'],...
%     'Resolution',[900,900]); 
% % imwrite(f.cdata,['./figures/sphere' num2str(n) '.png'],...
% %     'XResolution',600,'YResolution',600);
% save(['./figures/ball_' num2str(n) '_Hf.mat'],'X','Y','Hf');
% save(['./point clouds/XYZ_ball_' num2str(n) '.mat'],'X','Y','Z');
% save(['./volume and error results/ball_V_E_' num2str(n) '.txt'],'V_vector',...
%     'volume_vector','E_vector','alpha3D_V','error_alpha3D','volume_simp2',...
%     'error_simp','-ascii');
% % save for the civil 3D evaluation
% Px1=reshape(X,[],1);
% Py1=reshape(Y,[],1);
% Pz1=reshape(Z,[],1);
% comparison_surface=[Px1,Py1,Pz1];
% % figure;plot3(P1(:,1),P1(:,2),P1(:,3),'r*');
% % for the base surface points saving
% nnZ=size(Pz1,1);
% Z2=zeros(nnZ,1);
% base_surface=[Px1,Py1,Z2];
% 
% base_surface_100=base_surface.*100;
% comparison_surface_100=comparison_surface.*100;
% save(['./point clouds/base_surface_ball_100_' num2str(n) '.txt'],'base_surface_100','-ascii');
% save(['./point clouds/comparison_surface_ball_100_' num2str(n) '.txt'],'comparison_surface_100','-ascii');
% % 
% % comparison_cloud=pointCloud(comparison_surface);
% % B_C_cloud=pointCloud(B_C_surface);
% % pcwrite(comparison_cloud,'comparison_cloud_01','PLYFormat','binary');
% % pcwrite(B_C_cloud,'B_C_cloud_01','PLYFormat','binary');


