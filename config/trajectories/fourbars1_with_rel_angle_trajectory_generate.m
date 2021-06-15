dt= 5e-3;  % Timestep
N = 2000; % Number of time steps
%n = 5; % Number of dependent coordinates

%q_idx = 5; % Index of the "master" coordinate in "q" (1=first)

w = 2.0; % [rad/s]
M=zeros(N,1+1);
for i=1:N
    t = (i-1)*dt;
    M(i,1) = t;
    M(i,1+1) = 1.2*(1-cos(t*w));
end

save -ascii 'fourbars1-with-rel-angle-trajectory.txt' M;

stem(M(:,1+1));
