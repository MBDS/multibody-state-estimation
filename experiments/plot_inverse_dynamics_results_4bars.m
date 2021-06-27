close all;
clc;
clear;

%forceIdx = 1 + 2; % y1
%forceIdx = 1 + 3; % x2
forceIdx = 1 + 5; % theta

dir='.';
q=load(sprintf('%s/q.txt',dir));
Q=load(sprintf('%s/Q_forces.txt',dir));
refTrajectory=load('../config/trajectories/fourbars1-with-rel-angle-trajectory.txt');
MSC=load('InverseDynamics_4bars_MSC.tab');
MSC=MSC(1:(end-1),:);

afigure();
plot(q(:,1),q(:,6),'r','Marker','.','MarkerSize',9,'MarkerIndices',1:12:size(q,1));
hold on;
plot(MSC(:,1),pi/180*MSC(:,3),'b','Marker','o','MarkerSize',9,'MarkerIndices',4:12:size(q,1));
plot(q(:,1),refTrajectory(:,2),'k','Marker','s','MarkerSize',9,'MarkerIndices',8:12:size(q,1));
axis tight;
grid minor;
legend('FG','MSC','Reference');
xlabel('t [s]');
ylabel('\theta [rad]');

afigure();
plot(q(:,1),q(:,6)-refTrajectory(:,2),'r','Marker','.','MarkerSize',9,'MarkerIndices',6:12:size(q,1));
hold on;
plot(MSC(:,1),pi/180*MSC(:,3)-refTrajectory(:,2),'b','Marker','o','MarkerSize',9,'MarkerIndices',1:12:size(MSC,1));
axis tight;
legend('FG','MSC');
xlabel('t [s]');
ylabel('\theta - \theta_{ref} [rad]');
grid minor;

afigure();
plot(Q(:,1),Q(:,forceIdx),'r-','Marker','.','MarkerSize',9,'MarkerIndices',6:12:size(q,1));
hold on;
plot(MSC(:,1),-MSC(:,2),'b','Marker','o','MarkerSize',9,'MarkerIndices',1:12:size(MSC,1));
axis tight;
legend('FG','MSC');
xlabel('t [s]');
ylabel('Torque [N m]');
grid minor;



