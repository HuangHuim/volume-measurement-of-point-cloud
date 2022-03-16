close all;
clear;clc;
r=1; % the ridus of the cylinder
neg_r=(-1)*r;
a=1;b=3; % the sections of y
n=10; % partition number
% Equidistant partition
l1=abs(b-a);l2=abs(r-neg_r);
x=linspace(neg_r,r,l2*n+1);
y=linspace(a,b,l1*n+1);
[Y,X]=meshgrid(y,x);
Z=-1.*sqrt(1-X.^2);
Z=real(Z);
%% add two sections
% x=linspace(neg_r,r,l2*100+1);
M=size(Z,1);
M_centre=ceil(M/2);
Z_a=Z(1:M_centre-1,1);
[Xa,Za]=meshgrid(x,Z_a);
zz=-1.*sqrt(1-x.^2);

[ii,jj]=size(Xa);
for i=1:ii
    for j=1:jj
        if Xa(i,j)^2+Za(i,j)^2>1+1e-8
            Xa(i,j)=NaN;
            Za(i,j)=NaN;
        end
    end
end
Ya=ones(ii,jj);
Yb=3*ones(ii,jj);
% ii-1 is the number of patches, which equals to the number of partitions
% along Z axis
X_patch=zeros(4,ii-1);
Y_patch1=ones(4,ii-1);
Z_patch=zeros(4,ii-1);

for i=1:ii-1
    X_patch(1,i)=Xa(i,i);X_patch(2,i)=Xa(i+1,i+1);X_patch(3,i)=Xa(i+1,jj-i);X_patch(4,i)=Xa(i,jj-i+1);
    Z_patch(1,i)=Za(i,i);Z_patch(2,i)=Za(i+1,i+1);Z_patch(3,i)=Za(i+1,jj-i);Z_patch(4,i)=Za(i,jj-i+1);
end
% set the color for the cylinder
C=Z;
C(M_centre-1:M-1,:)=Z(M_centre:M,:);
figure;surf(X,Y,Z+1.1,C,'EdgeColor','none');
xlabel('X');ylabel('Y');zlabel('Z');
% set the color for the patch at the ends
C_patch=Z(1:ii-1,1);
hold on;
patch(X_patch,Y_patch1,Z_patch+1.1,C_patch,'EdgeColor','none');

Y_patch2=3*ones(4,ii-1);
hold on;
patch(X_patch,Y_patch2,Z_patch+1.1,C_patch,'EdgeColor','none');
hold on;
imagesc([-1;1],[1;3],Z');
axis tight;
save('XYZplot_cylinder_001.mat','X','Y','Z','C','X_patch','Y_patch1','Z_patch','Y_patch2','C_patch');
% colorbar;
% save the point cloud
Px1=reshape(X,[],1);
Py1=reshape(Y,[],1);
Pz1=reshape(Z,[],1);
P1=[Px1,Py1,Pz1];
% figure;plot3(P1(:,1),P1(:,2),P1(:,3),'r*');
% for the base surface points saving
nnZ=size(Pz1,1);
Z2=zeros(nnZ,1);
P1_base=[Px1,Py1,Z2];

Sxa=reshape(Xa,[],1);
Sya=reshape(Ya,[],1);
Sza=reshape(Za,[],1);
Psa=[Sxa,Sya,Sza];
indexA=~isnan(Sxa);
Psa=Psa(indexA,:);
% hold on;plot3(Psa(:,1),Psa(:,2),Psa(:,3),'r*');

Syb=reshape(Yb,[],1);
Psb=[Sxa,Syb,Sza];
Psb=Psb(indexA,:);
% hold on;plot3(Psb(:,1),Psb(:,2),Psb(:,3),'r*');

P_whole=[P1;Psa;Psb];
figure;plot3(P_whole(:,1),P_whole(:,2),P_whole(:,3),'b.');
% save('cylinder_001.mat','P_whole');
%% volume calculate

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
% filter 6*6
% H=[0 0 4 4 0 0;0 3 -7 -7 3 0;4 -7 3 3 -7 4;4 -7 3 3 -7 4;0 3 -7 -7 3 0;0 0 4 4 0 0];
% % figure;mesh(H);
% Hf=filter2(H,Z,'same');
% SS=l1*l2;  % the projected area
% Hf_full=filter2(H,Z,'full');
% mf=3; % the size of the filter
% [M_full,N_full]=size(Hf_full);
% Hf_mf=Hf_full(mf:M_full-mf+1,mf:N_full-mf+1);
% 
% dividend=4*12*4;  %four direction, four points, 12 is the common dividend
% V_Hf=sum(sum(Hf))*s/dividend;
% % in fact, below is the true 2th-order estimation
% averageH=sum(sum(Hf))/(M*N);
% V_Hf_averageH=averageH*SS/dividend;
% 
% nf=6;
% Hf_nf=Hf_full(nf:M_full-nf+1,nf:N_full-nf+1);
% [M_nf,N_nf]=size(Hf_nf);
% averageH_nf=sum(sum(Hf_nf))/(M_nf*N_nf);
% V_nf=averageH_nf*SS/dividend;
% 
% 
% [M_mf,N_mf]=size(Hf_mf);
% averageH_mf=sum(sum(Hf_mf))/(M_mf*N_mf);
% V_mf=averageH_mf*SS/dividend;
% 
% averageH_full=sum(sum(Hf_full))/(M_full*N_full);
% V_full=averageH_full*SS/dividend;
% 
% V_all_average2=(V_mf+V_full)/2;
% 
% figure;subplot(1,2,1);surf(X,Y,Hf,'EdgeColor','none');
% axis tight;axis off;
% % axis([-1,1,-1,1]);
% subplot(1,2,2);imagesc([-1;1],[-1;1],Hf);
% axis tight;axis off;
% save('XYHfplot_cylinder_001.mat','X','Y','Hf');
% 
% Hf1=reshape(Hf,[],1);
% Q3=prctile(Hf1,75);
% Q1=prctile(Hf1,25);
% Q2=prctile(Hf1,50);
% IQR=Q3-Q1;
% Upper_adjacent=Q3+1.5*IQR;
% Lower_adjacent=Q1-1.5*IQR;
% Hf_a=Hf1(Hf1<Upper_adjacent);
% Hf_a=Hf_a(Hf_a>Lower_adjacent);
% Mean_without_outliers=sum(Hf_a)/(M*N);
% V_Hf_mean=Mean_without_outliers*SS/dividend;
% 
% V_Hf_Q2=Q2*SS/dividend;
% V_all_average3=(V_Hf_Q2+V_mf+V_full)/3;
% V_Hf_averageV=Volume-V_Hf_averageH;
% 
% %% the volume calculation
% True_V=pi;
% base_surface=P1_base;
% comparison_surface_ends=P_whole;
% comparison_surface_noends=P1;
% 
% 
% B_C_surface=[base_surface;comparison_surface_ends];
% shp=alphaShape(B_C_surface);
% figure;plot(shp);
% alpha3D_V=volume(shp);
% error_alpha3D=abs(True_V-alpha3D_V);
% error_V1=abs(True_V+Volume);
% error_V2=abs(True_V+Volume-V_all_average3);
% V2=Volume-V_all_average3;
% 
% base_surface_100=base_surface.*1000;
% comparison_surface_100=comparison_surface_ends.*1000;
% save('base_surface_cylinder_1000_001.txt','base_surface_100','-ascii');
% save('comparison_surface_cylinder_1000_001.txt','comparison_surface_100','-ascii');
% 
% comparison_cloud=pointCloud(comparison_surface_ends);
% B_C_cloud=pointCloud(B_C_surface);
% pcwrite(comparison_cloud,'cylinder_comparison_cloud_001','PLYFormat','binary');
% pcwrite(B_C_cloud,'cylinder_B_C_cloud_001','PLYFormat','binary');






