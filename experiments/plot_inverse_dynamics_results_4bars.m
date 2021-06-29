close all;
clc;
clear;

forceIdx = 1 + 5; % theta index in vector Q[·]

prefix='./invdyn_4bars_';
q=load(sprintf('%sq.txt',prefix));
Q=load(sprintf('%sQ_forces.txt',prefix));
refTrajectory=load('../config/trajectories/fourbars1-with-rel-angle-trajectory.txt');
MSC=load('InverseDynamics_4bars_MSC.tab');
MSCtorque = MSC(:,2);
MSCang = pi/180*MSC(:,3);

afigure();
plot(q(:,1),180/pi*q(:,6),'r','Marker','.','MarkerSize',9,'MarkerIndices',1:30:size(q,1));
hold on;
plot(MSC(:,1),180/pi*MSCang,'b','Marker','o','MarkerSize',9,'MarkerIndices',10:30:size(q,1));
plot(q(:,1),180/pi*refTrajectory(:,2),'k','Marker','s','MarkerSize',9,'MarkerIndices',20:30:size(q,1));
axis tight;
grid minor;
legend('FG','MSC','Reference');
xlabel('t [s]');
ylabel('\theta [deg]');

afigure();
plot(q(:,1),180/pi*(q(:,6)-refTrajectory(:,2)),'r','Marker','.','MarkerSize',9,'MarkerIndices',1:25:size(q,1));
hold on;
plot(MSC(:,1),180/pi*(MSCang-refTrajectory(:,2)),'b','Marker','o','MarkerSize',9,'MarkerIndices',14:25:size(MSC,1));
axis tight;
legend('FG','MSC');
xlabel('t [s]');
ylabel('\theta - \theta_{ref} [deg]');
grid minor;


afigure();
plot(Q(:,1),Q(:,forceIdx),'r-','Marker','.','MarkerSize',9,'MarkerIndices',6:25:size(q,1));
hold on;
plot(MSC(:,1),-MSCtorque,'b--','Marker','o','MarkerSize',9,'MarkerIndices',19:25:size(MSC,1));
axis tight;
legend('FG','MSC');
xlabel('t [s]');
ylabel('Torque [N m]');
grid minor;



