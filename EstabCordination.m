% function [Anorm, Bnorm, Nnorm]=EstabCordination
%% 加载路面参数信息，将坑洼点云投影到路面平面
close all;clear;clc;
load roadsurface.mat
%% 投影坐标测试部分
pavement=model1.Parameters;
n1=pavement(1);
n2=pavement(2);
n3=pavement(3);
a1t=1;
a2t=-n1/(n2+n3);
a3t=a2t;
b1t=1;
b2t=(n3-n1*a3t)/(n2*a3t-n3*a2t);
b3t=(n1*a2t-n2)/(n2*a3t-n3*a2t);
A1=[n1;n2;n3];
A2=[a1t;a2t;a3t];
A3=[b1t;b2t;b3t];
A1_params=pavement./norm(A1);
A2_params=[A2;0]./norm(A2);
A3_params=[A3;0]./norm(A3);
figure;quiver3(0,0,0,A1_params(1),A1_params(2),A1_params(3));
hold on; quiver3(0,0,0,A2_params(1),A2_params(2),A2_params(3));
hold on; quiver3(0,0,0,A3_params(1),A3_params(2),A3_params(3));
title('road coordinate system');
points=load('potehole.mat');
Points=points.potehole;
% pca_normals=pca(Points);
add_cor=single(ones(size(Points,1),1));
Points=[Points,add_cor];
Z=Points*A1_params';
X=Points*A2_params;
Y=Points*A3_params;
Pro_points=[X,Y,Z];
Pro_points=pointCloud(Pro_points);
figure;pcshow(Pro_points);
X_grid=linspace(min(X),max(X),150);
Y_grid=linspace(min(Y),max(Y),170);
[X_v,Y_v]=meshgrid(X_grid,Y_grid);
X_v=double(X_v);
Y_v=double(Y_v);
l1=max(X)-min(X);
l2=max(Y)-min(Y);
X=double(X);
Y=double(Y);
Z=double(Z);
Z_grid=griddata(X,Y,Z,X_v,Y_v,'nearest');

Z=Z_grid;
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
%% the second-order estimation
SS=l1*l2;  % the projected area
True_V=Volume;
[V2,V_2021,V_2021_scale]=second_order(Z,s,SS,True_V,Volume);
% pcwrite(Pro_points,'Pro_potehole','PLYFormat','binary');
% three vectors products are zeros
% ProBN=Bnorm*pavement(1:3)';
% ProAB=Bnorm*Anorm';
% N=pavement(1:3);
% Nnorm=N./norm(N);