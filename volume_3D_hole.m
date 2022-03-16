
close all; clear; clc;
ptCloud = pcread('hole000.ply');
point=ptCloud.Location;
Points1=point(point(:,1)>0.27,:);
Points2=Points1(Points1(:,2)>-0.164,:);
figure;pcshow(Points2);
xlabel('X axis');ylabel('Y axis');zlabel('Z axis');

%% find the approximate direction of the road
% p1p3=[0.005 0.16004 0.0208]';
% p1p2=[0.146 -0.0048 -0.0207]';
% referenceVector=-cross(p1p3,p1p2);
referenceVector = [0.13,-0.13,0.99];
referenceVector=referenceVector/norm(referenceVector);

X=pointCloud(Points2);
maxDistance = 0.0005;
% the reference vector is the approximate direction we defined in advance
maxAngularDistance = 0.2;
[model1,inlierIndices,outlierIndices] = pcfitplane(X,...
    maxDistance,referenceVector);%,maxAngularDistance);
plane1 = select(X,inlierIndices);
remainPtCloud = select(X,outlierIndices);
figure;pcshow(plane1);
figure;pcshow(remainPtCloud);
xlabel('X');ylabel('Y');zlabel('Z');
remainPoints=remainPtCloud.Location;
potehole=remainPoints(remainPoints(:,1)<0.433,:);
potehole=potehole(potehole(:,2)<0.03,:);
poteholePt=pointCloud(potehole);
figure;pcshow(poteholePt);
% save('roadsurface.mat','model1');
save('potehole.mat','potehole');
% pcwrite(poteholePt,'poteholeB','PLYFormat','binary');
%% alpha shape volume
% point=potehole;
% x=point(:,1);y=point(:,2);z=point(:,3);%get point out
% x = double(x); y=double(y); z= double(z);
%     %获取点云坐标
% alp = 5;region = 0.75;%hole = 1; region = 0.75;
% shp = alphaShape(x,y,z,alp);
%     %生产点云的包络数据
% %ref:http://cn.mathworks.com/help/matlab/ref/alphashape.html
% figure;plot(shp)
%         %显示点云包络
% v= volume(shp);