close all;
clear; clc;
load('point3Ddisp.mat');
% �Ե��ƽ��н�����
X1=downsample(point3Ddisp,7);

X=pointCloud(X1);
% ������ȥ����Լ��������������Ҫ������·���ϣ�
% ���ԶԿ����������û��̫��Ӱ��
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
% �ڽ��о������Ҫ�Ե��ƽ�������޳���
% ����������Ҫ������·���Ϸ�������ֻ��Ҫ�ҳ���ƽ���·��ĵ�
% ������Ƶľ�ֵ�뷽�λ�����ĵ�����3sigma֮�������Ϊ����Ⱥ��
% ����ֱ�����õ����Ⱥ�����ĵ�ľ����������޳�������ƽ������ľ�������Ⱥ��
pavement=model1.Parameters;
firstPart=remainPtCloud.Location(idx==1,:);
secondPart=remainPtCloud.Location(idx==2,:);
thirdPart=remainPtCloud.Location(idx==3,:);

%% ѡȡ·���µĵ���
% ���ڵ�һ������
% ����������ת��Ϊ�������
n1 = size(firstPart,1);
homones=ones(n1,1);
FirP=[firstPart,homones];

% ƽ����η���ϵ������ pavement
% ��ƽ��ϵ�����С��0�ľ��ǿ����ڵĵ�
dis1=FirP*pavement';
hole1=firstPart(dis1<=0,:);

% ���ڵڶ�������
% ����������ת��Ϊ�������
n2 = size(secondPart,1);
homones2=ones(n2,1);
SecP=[secondPart,homones2];

dis2=SecP*pavement';
hole2=secondPart(dis2>=0,:);

% ���ڵڶ�������
% ����������ת��Ϊ�������
n3 = size(thirdPart,1);
homones3=ones(n3,1);
ThiP=[thirdPart,homones3];

dis3=ThiP*pavement';
hole3=thirdPart(dis3<=0,:);

save('hole1.mat','hole1');
save('hole2.mat','hole2');
save('plantParams.mat','pavement');
