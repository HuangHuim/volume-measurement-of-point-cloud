function [V2,V_2021,V_2021_scale]=second_order(Z,s,SS,True_V,Volume)

%% the second-order volume estimation
% filter 6*6
H=[0 0 4 4 0 0;0 3 -7 -7 3 0;4 -7 3 3 -7 4;4 -7 3 3 -7 4;0 3 -7 -7 3 0;0 0 4 4 0 0];
% figure;mesh(H);
Hf=filter2(H,Z,'same');
Hf_full=filter2(H,Z,'full');
mf=3; % the size of the filter
[M,N]=size(Z);
[M_full,N_full]=size(Hf_full);
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

% figure;subplot(1,2,1);surf(X,Y,Hf,'EdgeColor','none');
% axis tight;axis off;
% % axis([-1,1,-1,1]);
% subplot(1,2,2);imagesc([-1;1],[-1;1],Hf);
% axis tight;axis off;
% save('XYHfplot_bathtub_001.mat','X','Y','Hf');



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

% % filter matrix with size of 6*6
% H=[0 0 4 4 0 0;0 3 -7 -7 3 0;4 -7 3 3 -7 4;4 -7 3 3 -7 4;0 3 -7 -7 3 0;0 0 4 4 0 0];
% % figure;surf(H);
% Hf=filter2(H,Z,'same');
% Hf_full=filter2(H,Z,'full');
% dividend=4*12*4;  %four direction, four points, 12 is the common dividend
% 
% [M,N]=size(Z);
% [M_full,N_full]=size(Hf_full);
% 
% V_Hf=sum(sum(Hf))*s/dividend;
% % in fact, below is the true 2th-order estimation
% averageH=sum(sum(Hf))/((M-1)*(N-1));
% V_Hf_averageH=averageH*SS/dividend;
% 
% nf=6;
% Hf_nf=Hf_full(nf:M_full-nf+1,nf:N_full-nf+1);
% [M_nf,N_nf]=size(Hf_nf);
% averageH_nf=sum(sum(Hf_nf))/((M_nf-1)*(N_nf-1));
% V_nf=averageH_nf*SS/dividend;
% 
% averageH_full=sum(sum(Hf_full))/((M_full-1)*(N_full-1));
% V_full=averageH_full*SS/dividend;
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
% V_all_average1=(V_Hf+V_full+V_nf)/3;
% V_all_average2=(V_Hf_Q2+V_Hf)/2;
% V_all_average3=(V_Hf_Q2+V_Hf+V_full)/3;
% V_all_average4=(V_Hf_Q2+V_all_average1)/2;
% 
% volume_Hf=Volume-V_Hf;
% volume_full=Volume-V_full;
% volume_Hf_averageH=Volume-V_Hf_averageH;
% volume_nf=Volume-V_nf;
% volume_Hf_Q2=Volume-V_Hf_Q2;
% volume_Hf_mean=Volume-V_Hf_mean;
% volume_average1=Volume-V_all_average1;
% volume_average2=Volume-V_all_average2;
% volume_average3=Volume-V_all_average3;
% volume_average4=Volume-V_all_average4;
% 
% error1=True_V+Volume;
% error_Hf=True_V+Volume-V_Hf;
% error_full=True_V+Volume-V_full;
% error_Hf_averageH=True_V+Volume-V_Hf_averageH;
% error_nf=True_V+Volume-V_nf;
% error_Hf_Q2=True_V+Volume-V_Hf_Q2;
% error_Hf_mean=True_V+Volume-V_Hf_mean;
% error_average1=True_V+Volume-V_all_average1;
% error_average2=True_V+Volume-V_all_average2;
% error_average3=True_V+Volume-V_all_average3;
% error_average4=True_V+Volume-V_all_average4;
% 
% V_vector=[True_V,V_Hf,V_full,V_Hf_averageH,V_nf,V_Hf_Q2,V_Hf_mean,V_all_average1,...
%     V_all_average2,V_all_average3,V_all_average4];
% volume_vector=[Volume,volume_Hf,volume_full,volume_Hf_averageH,volume_nf,...
%     volume_Hf_Q2,volume_Hf_mean,volume_average1,volume_average2,...
%     volume_average4,volume_average3];
% E_vector=[error1,error_Hf,error_full,error_Hf_averageH,error_nf,error_Hf_Q2,...
%     error_Hf_mean,error_average1,error_average2,error_average3,error_average4];
% % figure;subplot(1,2,1);surf(X,Y,Hf,'EdgeColor','none');
% % axis tight;axis off;
% % % axis([-1,1,-1,1]);
% % subplot(1,2,2);imagesc([-1;1],[-1;1],Hf);
% % axis tight;axis off;

