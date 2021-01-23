time_los.t0      = time.list(length(time.list));
time_los.LOS     = 1.8 * 30 * 24 * 60 * 60; %1.6ヶ月
time_los.simDt   = 60;
time_los.stepNum = ceil(time_los.LOS/time_los.simDt);
time_los.list = linspace(time_los.t0,time_los.t0+time_los.simDt * time_los.stepNum,time_los.stepNum+1);

% 推定値初期値
scEstBySc_los.X = scEstByScEkf.X;
scEstByGs_los.X = scEstByGsEkf.X;
scEstBySc_los.P = scEstByScEkf.P;
scEstByGs_los.P = scEstByGsEkf.P;
scEstBySc_los.t = time_los.list;
scEstByGs_los.t = time_los.list;
% 真値
scTrue_los = Spacecraft(time_los);
scTrue_los.state = scTrue.state(:,size(scTrue.state,2));

% 真値の伝搬
scTrue_los.calcOrbitTwoBody(constant.sunMu,error.dynamics)

for i_los = 1:time_los.stepNum+1
    % 記録
    scEstByGs_los.clockError(i_los) = error.clock0- scEstByGs_los.X(1); %残りの誤差
    scEstByGs_los.state(:,i_los)    = scEstByGs_los.X(2:7);
    scEstByGs_los.P_list(:,:,i_los) = scEstByGs_los.P;
    scEstBySc_los.clockError(i_los) = error.clock0- scEstBySc_los.X(1); %残りの誤差
    scEstBySc_los.state(:,i_los)    = scEstBySc_los.X(2:7);
    scEstBySc_los.P_list(:,:,i_los) = scEstBySc_los.P;    
    % 時間伝搬
    [scEstBySc_los.X, scEstBySc_los.P] = Spacecraft.timeUpdateEkf(scEstBySc_los.X, scEstBySc_los.P, constant, time_los.simDt,time_los.simDt,error);
    [scEstByGs_los.X, scEstByGs_los.P] = Spacecraft.timeUpdateEkf(scEstByGs_los.X, scEstByGs_los.P, constant, time_los.simDt,time_los.simDt,error);
end

showResult(scTrue_los,scEstBySc_los,scEstByGs_los,error,n*3-2,[],gs,resultPath);

% 
