% Input 坑洼点以及路面平面参数向量
% Output 各个坑洼的体积及表面积
% 到这里坑洼点云周围还是有噪声点，基于对初始点云观察可知：
% 这是道路大坑洼周围的不规则坑洼造成的
% load hole1.mat;
% H1=pointCloud(hole1);
% figure;pcshow(H1);
close all;
clear;clc;
load potehole.mat;
H2=pointCloud(potehole);
figure;pcshow(H2);
% 通过后面的三角划分显示发现：hole2包含一些噪声点
% 也就是说需要对一些外点进行剔除
hole2=potehole;
Hole2Mean=mean(hole2);
Hole2Std=std(hole2);   % std为标准差，var为方差
n=size(hole2,1);
nhole2=hole2-repmat(Hole2Mean,n,1);
threshd=3.*Hole2Std;
Hole2=hole2(nhole2(:,1)<threshd(1) & nhole2(:,2)<threshd(2) & nhole2(:,3)<threshd(3),:);
H3=pointCloud(Hole2);
figure;pcshow(H3);
% hole2=double(hole2);
% establish coordination system
[Anorm, Bnorm, Nnorm]=EstabCordination;
axisP=[Anorm',Bnorm'];
CoorP=Hole2*axisP;
CoorP=double(CoorP);
TRI=delaunay(CoorP(:,1),CoorP(:,2));
figure;
triplot(TRI,CoorP(:,1),CoorP(:,2));

% To calculate the projected area, the edge of the pothole must be found
% the edge of the pothole 
HighP=Hole2*Nnorm';    % the height of the pothole
% 要想找到坑洼的边界点需结合投影坐标和深度坐标同时进行判断
% 最小深度点作为起点，利用投影坐标找到邻近点，
% 找到临近点中最高的点（除了之前已作为边界的点以外）作为边界点第二点
% 已经作为边界的点移出备选区域（不能移出，不然最后不知道是否已经找完所有的点），
% 依次类推，当找到的最后的点两个最高的点都在边界中时停止搜索。




% TRI=delaunayTriangulation(hole2(:,1),hole2(:,2),hole2(:,3));
% figure;tetramesh(TRI);
% 理论上应该写一个函数，计算三个点到一个面的体积
% 将三个点投影到路面上的三角形就是三柱体的底面，所以总的底面积就是
% 路面对应的空隙大小，而高就是坑洼点云的距离路面的‘加权’平均距离，
