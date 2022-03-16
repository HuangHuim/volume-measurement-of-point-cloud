% Input ���ݵ��Լ�·��ƽ���������
% Output �������ݵ�����������
% ��������ݵ�����Χ�����������㣬���ڶԳ�ʼ���ƹ۲��֪��
% ���ǵ�·�������Χ�Ĳ����������ɵ�
% load hole1.mat;
% H1=pointCloud(hole1);
% figure;pcshow(H1);
close all;
clear;clc;
load potehole.mat;
H2=pointCloud(potehole);
figure;pcshow(H2);
% ͨ����������ǻ�����ʾ���֣�hole2����һЩ������
% Ҳ����˵��Ҫ��һЩ�������޳�
hole2=potehole;
Hole2Mean=mean(hole2);
Hole2Std=std(hole2);   % stdΪ��׼�varΪ����
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
% Ҫ���ҵ����ݵı߽������ͶӰ������������ͬʱ�����ж�
% ��С��ȵ���Ϊ��㣬����ͶӰ�����ҵ��ڽ��㣬
% �ҵ��ٽ�������ߵĵ㣨����֮ǰ����Ϊ�߽�ĵ����⣩��Ϊ�߽��ڶ���
% �Ѿ���Ϊ�߽�ĵ��Ƴ���ѡ���򣨲����Ƴ�����Ȼ���֪���Ƿ��Ѿ��������еĵ㣩��
% �������ƣ����ҵ������ĵ�������ߵĵ㶼�ڱ߽���ʱֹͣ������




% TRI=delaunayTriangulation(hole2(:,1),hole2(:,2),hole2(:,3));
% figure;tetramesh(TRI);
% ������Ӧ��дһ�����������������㵽һ��������
% ��������ͶӰ��·���ϵ������ξ���������ĵ��棬�����ܵĵ��������
% ·���Ӧ�Ŀ�϶��С�����߾��ǿ��ݵ��Ƶľ���·��ġ���Ȩ��ƽ�����룬
