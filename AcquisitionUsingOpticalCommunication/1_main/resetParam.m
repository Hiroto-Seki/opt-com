function [time, gsTrue, scTrue,scEstBySc,scEstByGs] = resetParam(time, gsTrue, scTrue,scEstBySc,scEstByGs)
    % 結果の出力用に記録
    gsTrue.t_drList = [gsTrue.t_drList, gsTrue.t_dr(1:gsTrue.dr_counter)];
    gsTrue.snr_drList = [gsTrue.snr_drList,gsTrue.snr_dr(1:gsTrue.dr_counter)];
    scTrue.targetError_dtList = [scTrue.targetError_dtList,scTrue.targetError_dt(1:gsTrue.dr_counter) ];
    scTrue.pointError_dtList = [scTrue.pointError_dtList, scTrue.pointingError_dt(1:gsTrue.dr_counter)];

    time.lastSearch = 0;
    gsTrue.ut_counter = 0;
    gsTrue.t_ut = [];
    gsTrue.state_ut = [];
    gsTrue.direction_ut = [];
    gsTrue.directionAccuracy_ut = [];
    gsTrue.pointingError_ut = [];
    gsTrue.ut2w_counter = 0;
    gsTrue.ut2w_counterList = [];
    gsTrue.dr_counter=0;
    gsTrue.dr_observability = [];
    gsTrue.t_dr = [];
    gsTrue.state_dr = [];
    gsTrue.lengthTrue_dr = [];
    gsTrue.lengthObserved_dr = [];
    gsTrue.length2wObserved_dr = [];
    gsTrue.durationAtSc = [];
    gsTrue.directionTrue_dr = [];
    gsTrue.receivedPower_dr = [];
    gsTrue.directionAccuracy_dr = [];
    gsTrue.directionObserved_dr = [];
    gsTrue.scAccel_dr = [];
    gsTrue.scRecAngle_dr = [];
    gsTrue.scRecAngleAccuracy_dr = [];
    gsTrue.transUpAngle_dr = [];
    gsTrue.transUpAngleAccuracy_dr = [];
    scTrue.ur_counter = 0;
    scTrue.ur_observability = []; 
    scTrue.t_ur = [];
    scTrue.state_ur = [];
    scTrue.attStateTrue_ur = [];
    scTrue.attStateObserved_ur = [];
    scTrue.lengthTrue_ur = [];
    scTrue.lengthObserved_ur = [];
    scTrue.ur2w_counter = 0;
    scTrue.length2wObserved_ur = [];
    scTrue.durationAtGs = [];
    scTrue.directionTrue_ur = [];
    scTrue.directionObserved_ur = [];
    scTrue.directionAccuracy_ur = [];
    scTrue.receivedPower_ur = [];
    scTrue.accelTrue_ur = [];
    scTrue.accelObseved_ur = [];       %加速度の観測値
    scTrue.gsT_ur = [];                 %uplinkに載っている，送信時刻
    scTrue.eState_ur    = [];           %uplinkに載っている，送信時刻の地球の状態量   
    scTrue.gsState_ur   = [];           %uplinkに載っている, 送信時刻の地上局の状態量
    scTrue.transDirection_ur = [];       %uplinkに載っている，送信方向(誤差を持つ)
    scTrue.recDownAngle_ur = [];         %uplinkに載っている，地上局の観測量(downlinkの測角)
    scTrue.recDownAngleAccuracy_ur = []; %上の精度
    scTrue.dt_counter  = 0;           % downlinkを送信した回数を数える 
    scTrue.t_dt    = [];                % downlinkを送信した時刻  
    scTrue.state_dt  = [];              % downlinkを送信した時刻の状態量(観測量の計算に用いる)
    scTrue.accel_dt    = [];            % downlinkに載っている，宇宙機の観測量(宇宙機の加速度)
    scTrue.recUpAngle_dt  = [];         % downlinkに載っている，宇宙機の観測量(uplinkの受信角度)
    scTrue.recUpAngleAccuracy_dt= [];   % 上の精度
    scTrue.transUpAngle_dt  = [];       % downlinkに載っている，宇宙機の観測量(uplinkの送信角度)
    scTrue.transUpAngleAccuracy_dt = []; % 上の精度
    scTrue.pointingError_dt = [];       % downlinkの送信方向誤差
    scTrue.tRec_dt  = [];               % downlinkが受信されるはずの時刻 
    % X_dt, P_dt　(scEstBySc)
    scEstBySc.X_dt = [];
    scEstBySc.P_dt = [];
    scEstByGs.X_dt = [];
    scEstByGs.P_dt = [];
    
    scEstBySc.estNoUse =[];
    scEstByGs.estNoUse =[];
    
end