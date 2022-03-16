% a convex ball with radius 3 and a concanvo ball with radius 2
close all;clear;clc;
%% generate the point cloud
R1=3;R2=2; % the ridus of two balls overlapping
h1=2*sqrt(2);
h2=sqrt(3);
t1=3+h1;
t2=3+h1+h2;
a=0;b=t2+2;c=0;d=6;
l1=abs(b-a);l2=abs(d-c);
n=10; % partition number
% Equidistant partition
x=linspace(c,d,l2*n+1);
y=linspace(a,b,l1*n+1);
[X,Y]=meshgrid(x,y);
[M,N]=size(X);
Z=zeros(M,N);
for i=1:M  % y label
    for j=1:N
        t=(X(i,j)-3)^2+(Y(i,j)-3)^2;
        if Y(i,j)<t1
            if t<1
                Z(i,j)=sqrt(1-t)-h1;
            else
                Z(i,j)=-sqrt(9-t);
            end
        else
            Z(i,j)=(-1)*sqrt(4-(X(i,j)-3)^2-(Y(i,j)-t2)^2);
        end
    end
end
Z=real(Z);

% ball_01=[X,Y,Z];
%% the second-order estimation
% delt_r=1/n;
% s=delt_r*delt_r;
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
% % filter 6*6
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
% V_Hf_averageV=Volume-V_Hf_averageH;
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
% figure;surf(Hf_mf,'EdgeColor','none');
% colorbar;
% xlabel('X');ylabel('Y');zlabel('Z');
% title('Hf_mf with the full filter');
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
% figure;surf(X,Y,Z+3,'EdgeColor','none');
% hold on
% imagesc([c;d],[a;b],Z+3);
% % axis([-1,1,-1,1]);
% xlabel('X');ylabel('Y');zlabel('Z');
% axis tight;axis off;
% % title('Hf_mf with the full filter');
% 
% figure;subplot(1,2,1);surf(X,Y,Hf,'EdgeColor','none');
% axis tight;axis off;
% % axis([-1,1,-1,1]);
% subplot(1,2,2);imagesc([c;d],[a;b],Hf);
% axis tight;axis off;
% 
% save('XYZplot_CC_gourd_001.mat','X','Y','Z');
% save('XYHfplot_CC_gourd_001.mat','X','Y','Hf');
% %% the 3D alpha volume estimation and errors calculation
% True_V=19*sqrt(2)*pi+3*sqrt(3)*pi/2-7*pi;
% V2=Volume-V_all_average3;
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
% save('base_surface_CC_gourd_100_001.txt','base_surface_100','-ascii');
% save('comparison_surface_CC_gourd_100_001.txt','comparison_surface_100','-ascii');
% 
% comparison_cloud=pointCloud(comparison_surface);
% B_C_cloud=pointCloud(B_C_surface);
% pcwrite(comparison_cloud,'CC_gourd_comparison_cloud_001','PLYFormat','binary');
% pcwrite(B_C_cloud,'CC_gourd_B_C_cloud_001','PLYFormat','binary');






