close all;
clc;
clear;

forceIdx = 1 + 5; % theta index in vector Q[·]

prefix='./invdyn_4bars_';
q=load(sprintf('%sq.txt',prefix));
Q=load(sprintf('%sQ_forces.txt',prefix));
refTrajectory=load('../config/trajectories/fourbars1-with-rel-angle-trajectory.txt');
MSC=load('InverseDynamics_4bars_MSC.tab');
MSC=MSC(1:(end-1),:);

if (0)
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
end

if (0)
afigure();
plot(q(:,1),q(:,6)-refTrajectory(:,2),'r','Marker','.','MarkerSize',9,'MarkerIndices',6:12:size(q,1));
hold on;
plot(MSC(:,1),pi/180*MSC(:,2)-refTrajectory(:,2),'b','Marker','o','MarkerSize',9,'MarkerIndices',1:12:size(MSC,1));
axis tight;
legend('FG','MSC');
xlabel('t [s]');
ylabel('\theta - \theta_{ref} [rad]');
grid minor;
end

afigure();
plot(Q(:,1),Q(:,forceIdx),'r-','Marker','.','MarkerSize',9,'MarkerIndices',6:12:size(q,1));
hold on;
plot(MSC(:,1),-MSC(:,2),'b--','Marker','o','MarkerSize',9,'MarkerIndices',1:12:size(MSC,1));
axis tight;
legend('FG','MSC');
xlabel('t [s]');
ylabel('Torque [N m]');
grid minor;



