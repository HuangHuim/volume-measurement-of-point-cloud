% the shape of the potholl is the playground shape
close all;clear;clc;
%% generate the point cloud
R1=2;R2=1; % the ridus of two balls overlapping
r=1/2;
t1=2+sqrt(15)/2;
t2=t1+sqrt(3)/2;
a=0;b=t2+1;c=0;d=4;
l1=abs(b-a);l2=abs(d-c);
n=100; % partition number
% Equidistant partition
x=linspace(c,d,l2*n+1);
y=linspace(a,b,l1*n+1);
[X,Y]=meshgrid(x,y);
[M,N]=size(X);
Z=zeros(M,N);
for i=1:M  % y label
    for j=1:N
        if Y(i,j)<t1
            Z(i,j)=(-1)*sqrt(4-(X(i,j)-2)^2-(Y(i,j)-2)^2);
        else
            Z(i,j)=(-1)*sqrt(1-(X(i,j)-2)^2-(Y(i,j)-t2)^2);
        end
    end
end


% Z=-1.*sqrt(1-X.^2-Y.^2);
Z=real(Z);
% figure;surf(X,Y,Z,'EdgeColor','none');
% colorbar;
% xlabel('X');ylabel('Y');zlabel('Z');
% ball_01=[X,Y,Z];
%% the first-order estimation
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
% figure;surf(Hf,'EdgeColor','none');
% colorbar;
% xlabel('X');ylabel('Y');zlabel('Z');
% 
% figure;subplot(1,2,1);surf(X,Y,Hf,'EdgeColor','none');
% axis tight;axis off;
% % axis([-1,1,-1,1]);
% subplot(1,2,2);imagesc([-1;1],[-1;1],Hf);
% axis tight;axis off;
% save('XYZplot_gourd_001.mat','X','Y','Z');
% save('XYHfplot_gourd_001.mat','X','Y','Hf');
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
% 
% %% the error estimation and save the data for ploting and calculation
% True_V=3*pi+sqrt(15)*pi+sqrt(3)*pi/4-5*sqrt(15)*pi/16-sqrt(3)*pi/16;
% V2=Volume-V_all_average3;
% 
% V_Hf_averageV=Volume-V_Hf_averageH;
% 
% 
% error1=True_V+Volume;
% error2=True_V+Volume-V_all_average3;
% 
% X1=reshape(X,[],1);
% Y1=reshape(Y,[],1);
% Z1=reshape(Z,[],1);
% nnZ=size(Z1,1);
% Z2=zeros(nnZ,1);
% base_surface=[X1,Y1,Z2];
% comparison_surface=[X1,Y1,Z1];
% B_C_surface=[base_surface;comparison_surface];
% shp=alphaShape(B_C_surface);
% figure;plot(shp);
% alpha3D_V=volume(shp);
% error_alpha3D=abs(True_V-alpha3D_V);
% 
% % save('base_surface_ball_01.txt','base_surface','-ascii');
% % save('comparison_surface_ball_01.txt','comparison_surface','-ascii');
% 
% base_surface_100=base_surface.*100;
% comparison_surface_100=comparison_surface.*100;
% save('base_surface_gourd_100_001.txt','base_surface_100','-ascii');
% save('comparison_surface_gourd_100_001.txt','comparison_surface_100','-ascii');
% 
% comparison_cloud=pointCloud(comparison_surface);
% B_C_cloud=pointCloud(B_C_surface);
% pcwrite(comparison_cloud,'gourd_comparison_cloud_001','PLYFormat','binary');
% pcwrite(B_C_cloud,'gourd_B_C_cloud_001','PLYFormat','binary');








