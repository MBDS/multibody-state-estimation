%%
close all;
clc;
clear;

force1Idx = 1 + 21;  % indices in vector Q[·]
force2Idx = 1 + 22;

prefix='./invdyn_pprobot_';
q=load(sprintf('%sq.txt',prefix));
Q=load(sprintf('%sQ_forces.txt',prefix));
refTrajectory=load('../config/trajectories/pick-and-place-robot-trajectory.txt');
MSC=load('InverseDynamics_pprobot_MSC.tab');
MSC=MSC(1:(end-1),:);


%%
afigure();
plot(q(:,1),180/pi*q(:,force1Idx),'r','Marker','.','MarkerSize',9,'MarkerIndices',1:20:size(q,1));
hold on;
plot(q(:,1),180/pi*q(:,force2Idx),'r','Marker','.','MarkerSize',9,'MarkerIndices',1:20:size(q,1));
%plot(MSC(:,1),pi/180*MSC(:,3),'b','Marker','o','MarkerSize',9,'MarkerIndices',4:12:size(q,1));
plot(refTrajectory(:,1),180/pi*refTrajectory(:,2),'k--','Marker','s','MarkerSize',9,'MarkerIndices',10:20:size(q,1));
plot(refTrajectory(:,1),180/pi*refTrajectory(:,3),'k--','Marker','s','MarkerSize',9,'MarkerIndices',10:20:size(q,1));
axis tight;
grid minor;
legend('FG \theta_1','FG \theta_2', 'Reference \theta_1','Reference \theta_2');
xlabel('t [s]');
ylabel('\theta [deg]');

%%
afigure();
plot(q(:,10),q(:,11),'b','Marker','.','MarkerSize',9);
hold on;
plot([-0.5 -0.5 0.5],[-5 -4 -4],'k:','LineWidth',4);
xlim([-1.5 1.5]);
ylim([-6 -4]);
axis equal;
xlabel('x_6');
ylabel('y_6');
grid minor;

%%
afigure();
plot(Q(:,1),Q(:,force1Idx),'r-','Marker','.','MarkerSize',9,'MarkerIndices',1:20:size(q,1));
hold on;
plot(MSC(:,1),MSC(:,2),'b--','Marker','o','MarkerSize',9,'MarkerIndices',10:20:size(MSC,1));
axis tight;
legend('FG \tau_1','MSC \tau_1');
xlabel('t [s]');
ylabel('Torque [N m]');
grid minor;

afigure();
plot(Q(:,1),Q(:,force2Idx),'g-','Marker','.','MarkerSize',9,'MarkerIndices',1:20:size(q,1));
hold on;
plot(MSC(:,1),MSC(:,3),'k--','Marker','o','MarkerSize',9,'MarkerIndices',10:20:size(MSC,1));
axis tight;
legend('FG \tau_2', 'MSC \tau_2');
xlabel('t [s]');
ylabel('Torque [N m]');
grid minor;



