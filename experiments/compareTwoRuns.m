% Example usage:
%
% compareTwoRuns('mech_4bars_MSC_1ms.tab', 'forward_experiment1/bars4_dt=0.001_lag=0.002_iters=15', 1)
%
function [rmse_q, rmse_dq] = compareTwoRuns(MSC_FILE,FG_PREFIX, DO_PLOTS, ADAM_DECIMATION)

    % MSC ADAMS
    % "Time" "x1" "x1dot" "x2" "x2dot" "y1"	"y1dot"	"y2"  "y2dot"
    %  1     2     3        4    5      6     7       8    9
	D=load(MSC_FILE);
    D=D(1:ADAM_DECIMATION:end,:);
    % params: loadOneDataset(D,idxX1,idxY1,idxX2,idxY2, idxXd1,idxYd1,idxXd2,idxYd2)
	[ts_a,x1_a,y1_a,x2_a,y2_a, dx1_a,dy1_a,dx2_a,dy2_a]=loadOneDataset(D,2,6,4,8, 3,7,5,9);

    % FG
	D_q=load(sprintf('%s_q.txt',FG_PREFIX));
	D_dq=load(sprintf('%s_dq.txt',FG_PREFIX));
	D_ddq=load(sprintf('%s_ddq.txt',FG_PREFIX));
    n = size(D_q,2);
	[ts_f,x1_f,y1_f,x2_f,y2_f,dx1_f,dy1_f,dx2_f,dy2_f]=loadOneDataset([D_q, D_dq,D_ddq],2,3,4,5, n+2,n+3,n+4,n+5);
    
    if (DO_PLOTS)
        D=10;  % decimation to make PFG figures smaller
        
        afigure();
        set(gcf,'Name','q vector comparison');
        subplot(4,1,1); plot(ts_a(1:D:end),x1_a(1:D:end),'r.', ts_f(1:D:end),x1_f(1:D:end),'b.'); title('x1'); legend('MSC.ADAMS','Ours'); xlabel('time [s]'); ylabel('[m]');
        subplot(4,1,2); plot(ts_a(1:D:end),y1_a(1:D:end),'r.', ts_f(1:D:end),y1_f(1:D:end),'b.'); title('y1'); legend('MSC.ADAMS','Ours'); xlabel('time [s]'); ylabel('[m]');
        subplot(4,1,3); plot(ts_a(1:D:end),x2_a(1:D:end),'r.', ts_f(1:D:end),x2_f(1:D:end),'b.'); title('x2'); legend('MSC.ADAMS','Ours'); xlabel('time [s]'); ylabel('[m]');
        subplot(4,1,4); plot(ts_a(1:D:end),y2_a(1:D:end),'r', ts_f(1:D:end),y2_f(1:D:end),'b.'); title('y2'); legend('MSC.ADAMS','Ours'); xlabel('time [s]'); ylabel('[m]');

        afigure();
        set(gcf,'Name','dq vector comparison');
        subplot(4,1,1); plot(ts_a(1:D:end),dx1_a(1:D:end),'r.', ts_f(1:D:end),dx1_f(1:D:end),'b.'); title('dx1'); legend('MSC.ADAMS','Ours'); xlabel('time [s]'); ylabel('[m/s]');
        subplot(4,1,2); plot(ts_a(1:D:end),dy1_a(1:D:end),'r.', ts_f(1:D:end),dy1_f(1:D:end),'b.'); title('dy1'); legend('MSC.ADAMS','Ours'); xlabel('time [s]'); ylabel('[m/s]');
        subplot(4,1,3); plot(ts_a(1:D:end),dx2_a(1:D:end),'r.', ts_f(1:D:end),dx2_f(1:D:end),'b.'); title('dx2'); legend('MSC.ADAMS','Ours'); xlabel('time [s]'); ylabel('[m/s]');
        subplot(4,1,4); plot(ts_a(1:D:end),dy2_a(1:D:end),'r.', ts_f(1:D:end),dy2_f(1:D:end),'b.'); title('dy2'); legend('MSC.ADAMS','Ours'); xlabel('time [s]'); ylabel('[m/s]');

        % Errors MSC <-> FG
        afigure();
        set(gcf,'Name','Trajectory errors');
        if (0)
            subplot(4,1,1); plot(ts_a,x1_a-x1_f,'r.'); ylim([-0.02 0.02]); title('x1 error'); xlabel('time [s]'); ylabel('[m]');
            subplot(4,1,2); plot(ts_a,y1_a-y1_f,'r.'); ylim([-0.02 0.02]); title('y1 error'); xlabel('time [s]'); ylabel('[m]');
            subplot(4,1,3); plot(ts_a,x2_a-x2_f,'r.'); ylim([-0.02 0.02]); title('x2 error'); xlabel('time [s]'); ylabel('[m]');
            subplot(4,1,4); plot(ts_a,y2_a-y2_f,'r.'); ylim([-0.02 0.02]); title('y2 error'); xlabel('time [s]'); ylabel('[m]');
        else
            clf;
            hold on;
            plot(ts_a(1:D:end),x1_a(1:D:end)-x1_f(1:D:end),'k-'); 
            plot(ts_a(1:D:end),y1_a(1:D:end)-y1_f(1:D:end),'r--');
            plot(ts_a(1:D:end),x2_a(1:D:end)-x2_f(1:D:end),'b:');
            plot(ts_a(1:D:end),y2_a(1:D:end)-y2_f(1:D:end),'g-.');
            legend('x1','y1','x2','y2');
            
            xlabel('time [s]'); ylabel('Error in q [m]');
        end
        
        % vel errors MSC <-> FG
        afigure();
        set(gcf,'Name','Trajectory velocity errors');
        clf;
        hold on;
        plot(ts_a(1:D:end),dx1_a(1:D:end)-dx1_f(1:D:end),'k-'); 
        plot(ts_a(1:D:end),dy1_a(1:D:end)-dy1_f(1:D:end),'r--');
        plot(ts_a(1:D:end),dx2_a(1:D:end)-dx2_f(1:D:end),'b:');
        plot(ts_a(1:D:end),dy2_a(1:D:end)-dy2_f(1:D:end),'g-.');
        legend('dx1','dy1','dx2','dy2');

        xlabel('time [s]'); ylabel('Error in dq [m/s]');

    end
    
    errs_q = [...
        x1_a-x1_f; ...
        y1_a-y1_f; ...
        x2_a-x2_f; ...
        y2_a-y2_f; ...
        ];
    rmse_q=sqrt(sum(errs_q.^2)/length(errs_q));

    errs_dq = [...
        dx1_a-dx1_f; ...
        dy1_a-dy1_f; ...
        dx2_a-dx2_f; ...
        dy2_a-dy2_f; ...
        ];
    rmse_dq=sqrt(sum(errs_dq.^2)/length(errs_dq));
    
    if (DO_PLOTS)
        fprintf('RMSE total in q  = %g\n',rmse_q);
        fprintf('RMSE total in dq = %g\n',rmse_dq);
    end
end

function [ts,x1,y1,x2,y2, dx1,dy1,dx2,dy2]=loadOneDataset(D,idxX1,idxY1,idxX2,idxY2, idxXd1,idxYd1,idxXd2,idxYd2)
	ts=D(:,1);
	x1=D(:,idxX1);
	y1=D(:,idxY1);
	x2=D(:,idxX2);
	y2=D(:,idxY2);
	
    dx1=D(:,idxXd1);
	dy1=D(:,idxYd1);
	dx2=D(:,idxXd2);
	dy2=D(:,idxYd2);

	
    if (0)
        figure('Name',sprintf('Trajectory-%s',name));
        plot(x1,y1,'b.'); hold on;
        plot(x2,y2,'r.');
        axis equal;
        legend('Point 1', 'Point 2');

        
        figure('Name',sprintf('Coordinates-%s',name));

        idx=1; d=x1; t='x1';
        subplot(2,2, idx); plot(ts,d,'.'); title(t);

        idx=2; d=y1; t='y1';
        subplot(2,2, idx); plot(ts,d,'.'); title(t);

        idx=3; d=x2; t='x2';
        subplot(2,2, idx); plot(ts,d,'.'); title(t);

        idx=4; d=y2; t='y2';
        subplot(2,2, idx); plot(ts,d,'.'); title(t);
    end
end

