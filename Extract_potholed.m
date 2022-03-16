close all;
clear; clc;
load('point3Ddisp.mat');
% 对点云进行降采样
X1=downsample(point3Ddisp,7);

X=pointCloud(X1);
% 不进行去燥可以减少误差，而且误差主要集中在路面上，
% 所以对坑洼体积计算没有太大影响
% X2=pointCloud(X1);
% X=pcdenoise(X2,'Threshold',0.02);
% pcshow(X2);title('origin downsample point cloud');
% figure;pcshow(X);title('denoised downsample point cloud');
maxDistance = 0.02;
referenceVector = [0,0,1];
maxAngularDistance = 5;
[model1,inlierIndices,outlierIndices] = pcfitplane(X,...
            maxDistance,referenceVector,maxAngularDistance);
plane1 = select(X,inlierIndices);
remainPtCloud = select(X,outlierIndices);
figure;pcshow(plane1);
figure;pcshow(remainPtCloud);
idx=kmeans(remainPtCloud.Location,3);
figure;plot3(remainPtCloud.Location(idx==1,1),remainPtCloud.Location(idx==1,2),...
    remainPtCloud.Location(idx==1,3),'r*');
figure;plot3(remainPtCloud.Location(idx==2,1),remainPtCloud.Location(idx==2,2),...
    remainPtCloud.Location(idx==2,3),'y*');
figure;plot3(remainPtCloud.Location(idx==3,1),remainPtCloud.Location(idx==3,2),...
    remainPtCloud.Location(idx==3,3),'b*');
% 在进行聚类后需要对点云进行外点剔除，
% 由于噪声主要集中在路面上方，所以只需要找出在平面下方的点
% 计算点云的均值与方差，位于中心点正负3sigma之间的数据为非离群点
% 可以直接利用点与点群的中心点的距离来进行剔除，大于平均距离的就算是离群点
pavement=model1.Parameters;
firstPart=remainPtCloud.Location(idx==1,:);
secondPart=remainPtCloud.Location(idx==2,:);
thirdPart=remainPtCloud.Location(idx==3,:);

%% 选取路面下的点云
% 对于第一个点云
% 将点云坐标转化为齐次坐标
n1 = size(firstPart,1);
homones=ones(n1,1);
FirP=[firstPart,homones];

% 平面齐次方程系数向量 pavement
% 与平面系数相乘小于0的就是坑洼内的点
dis1=FirP*pavement';
hole1=firstPart(dis1<=0,:);

% 对于第二个点云
% 将点云坐标转化为齐次坐标
n2 = size(secondPart,1);
homones2=ones(n2,1);
SecP=[secondPart,homones2];

dis2=SecP*pavement';
hole2=secondPart(dis2>=0,:);

% 对于第二个点云
% 将点云坐标转化为齐次坐标
n3 = size(thirdPart,1);
homones3=ones(n3,1);
ThiP=[thirdPart,homones3];

dis3=ThiP*pavement';
hole3=thirdPart(dis3<=0,:);

save('hole1.mat','hole1');
save('hole2.mat','hole2');
save('plantParams.mat','pavement');
