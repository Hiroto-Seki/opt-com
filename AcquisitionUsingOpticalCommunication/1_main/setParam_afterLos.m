% 長時間の通信途絶後の軌道決定

function [time_afterLos,gsTrue_afterLos,earth_afterLos,scTrue_afterLos,scEstByScSeq_afterLos,scEstByGsSeq_afterLos]...
    = setParam_afterLos(time,time_los,gs,scTrue_los,scEstBySc_los,scEstByGs_los,scEstByScEkf,scEstByGsEkf,constant,error)

time_afterLos    = time;
time_afterLos.t0 = time_los.list(end);
time_afterLos.stepNum = 10000;
time_afterLos.list = linspace(time_afterLos.t0,time_afterLos.t0+time_afterLos.simDt * time_afterLos.stepNum,time_afterLos.stepNum+1);

% 地上局
gsTrue_afterLos = GroundStation(gs,constant,time_afterLos);
% 地球
earth_afterLos  = CelestialBody(time_afterLos,"Earth");
earth_afterLos.getEphem(time_afterLos);
% 宇宙機(真値)
scTrue_afterLos = Spacecraft(time_afterLos);
scTrue_afterLos.state = scTrue_los.state(:,end);


%% ここから
scEstByScSeq_afterLos            = Spacecraft(time_afterLos);
scEstByScSeq_afterLos.state      = scEstBySc_los.state(:,end);
scEstByScSeq_afterLos.clockError = scEstBySc_los.clockError(end); 
% 地上局がEKFで推定した値
scEstByGsSeq_afterLos            = Spacecraft(time_afterLos);
scEstByGsSeq_afterLos.state      = scEstByGs_los.state(:,end);
scEstByGsSeq_afterLos.clockError = scEstByGs_los.clockError(end); 

%% 推定値 (clockのオフセット + 宇宙機の位置・速度)
scEstByScSeq_afterLos.X             = scEstBySc_los.X;
scEstByScSeq_afterLos.P             = scEstBySc_los.P;
% scEstByScSeq_afterLos.P             = [error.clockSigma^2,                                             zeros(1,6);
%                                        zeros(3,1), 1/3 * error.scPosSigma^2 * eye(3),                  zeros(3,3);
%                                                                     zeros(3,4), 1/3* error.scVelSigma^2 * eye(3)];

scEstByScSeq_afterLos.P_list        = zeros(7,7,length(time_afterLos.list));
scEstByScSeq_afterLos.P_list(:,:,1) = scEstByScSeq_afterLos.P;
%%  ここ
scEstByGsSeq_afterLos.X             = scEstByGs_los.X;
scEstByGsSeq_afterLos.P             = scEstByGs_los.P;
% scEstByGsSeq_afterLos.P             = [error.clockSigma^2,                                             zeros(1,6);
%                                        zeros(3,1), 1/3 * error.scPosSigma^2 * eye(3),                  zeros(3,3);
%                                                                     zeros(3,4), 1/3* error.scVelSigma^2 * eye(3)];

scEstByGsSeq_afterLos.P_list        = zeros(7,7,length(time_afterLos.list));
scEstByGsSeq_afterLos.P_list(:,:,1) = scEstByGsSeq_afterLos.P;


% 使う観測の設定0or1で記述する. 0は使用しない. 1は使用する
scEstByScSeq_afterLos.useObs = scEstByScEkf.useObs;
scEstByGsSeq_afterLos.useObs = scEstByGsEkf.useObs;



scEstByScSeq_afterLos.R = scEstByScEkf.R;
scEstByGsSeq_afterLos.R = scEstByGsEkf.R;




% 送受信した回数の初期化
gsTrue_afterLos.ut_counter=0;
gsTrue_afterLos.dr_counter=0;
gsTrue_afterLos.ut2w_counter = 0;
scTrue_afterLos.ur_counter=0;
scTrue_afterLos.dt_counter=0;
scTrue_afterLos.ur2w_counter = 0;
% 初期化
time_afterLos.lastSearch = 0;

end