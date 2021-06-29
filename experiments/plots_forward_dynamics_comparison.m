% Experiment: eval the effect of the lag for a fixed no. of iterations:
% *** See README.md first ***

function [] = plots_forward_dynamics_comparison(basename)

if (1)
    lags=[0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.010];
    smootheriterations=15;

    rmse_q=zeros(1,length(lags));
    rmse_dq=zeros(1,length(lags));

    for i=1:length(lags)
        lag = lags(i);
        [rmse_q(i),rmse_dq(i)]=compareTwoRuns(...
            'mech_4bars_MSC_0.1ms.tab', ...
            sprintf('%s_dt=0.001_lag=%.03f_iters=%i',basename, lag, smootheriterations), ...
            1,...  % PLOTS
            10 ...  % adams decimation
            );
    end
    format long
    rmse_q
    rmse_dq

    afigure(); abar(lags*1e3, rmse_q); title('rmse in q'); xlabel('Window length [ms]');
    afigure(); abar(lags*1e3, rmse_dq); title('rmse in dq'); xlabel('Window length [ms]');
end

% Experiment 2: eval the effect of the no. of iterations for fixed lag
%
if (0)
    smootheriterationss=[3 4 5 6 7 8 9 10];
    lag=0.005;

    rmse_q=zeros(1,length(smootheriterationss));
    rmse_dq=zeros(1,length(smootheriterationss));

    for i=1:length(smootheriterationss)
        smootheriterations = smootheriterationss(i);
        [rmse_q(i),rmse_dq(i)]=compareTwoRuns(...
            'mech_4bars_MSC_0.1ms.tab', ...
            sprintf('forward_experiment2/bars4_dt=0.001_lag=%.03f_iters=%i',lag, smootheriterations), ...
            0,...  % PLOTS
            10 ...  % adams decimation
            );
    end

    afigure(); abar(smootheriterationss*1e3, rmse_q); title('rmse in q'); xlabel('Max. iterations');
    afigure(); abar(smootheriterationss*1e3, rmse_dq); title('rmse in dq'); xlabel('Max. iterations');
end

% Like experiment 1, for indep. coord graph:
if (0)
    lags=[0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.010]; % 0.011 0.012 0.013 0.014 0.015 0.016 0.017 0.018 0.019 0.020];
    smootheriterations=15;

    rmse_q=zeros(1,length(lags));
    rmse_dq=zeros(1,length(lags));

    for i=1:length(lags)
        lag = lags(i);
        [rmse_q(i),rmse_dq(i)]=compareTwoRuns(...
            'mech_4bars_MSC_0.1ms.tab', ...
            sprintf('forward_icoords_experiment1/bars4_dt=0.001_lag=%.03f_iters=%i',lag, smootheriterations), ...
            0,...  % PLOTS
            10 ...  % adams decimation
            );
        %pause;
    end
    rmse_q

    afigure(); abar(lags*1e3, rmse_q); title('rmse in q'); xlabel('Window length [ms]');
    afigure(); abar(lags*1e3, rmse_dq); title('rmse in dq'); xlabel('Window length [ms]');
end
