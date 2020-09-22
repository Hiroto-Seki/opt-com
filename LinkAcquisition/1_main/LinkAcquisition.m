%{
---  Contents   ---
end to end link acquisition simulator
---  outline    ---
1. get ephemris date of earth, ground station and spacecraft
2. calculate true orbit of  spacecraft
3. calculate initial guess of spacecraft orbit (not true)
4. ground station Send a signal around the estimated orbit to the spacecraft according to the search pattern.
5. calculate observed value using the send signal
6. estimate spacecraft orbit using observed value by EKF
7. send back a signal to the ground station
8. gruond station observe a signal from spacecraft
---   Date      ---
rev1: 2020/05/29
%}

%% 前処理
% clear all; close all; clc
% add path to SPICE
addpath(genpath('~/Documents/Matlab/SPICE'));
% SPICEのKernel(天体情報)を読み込む
spice_loadkernels();
% 取得した天体情報 + alphaを利用しやすいように構造体へまとめる
SSD = spice_setparams();
% 乱数
rng('default');
rng(1)

%% setting parameter
[constant,time,error,gs,sc,sc_est,sc_estGs] = setparam(SSD);

%% 1 set gorund station, earth 
% 地上局
gsTrue = GroundStation(gs,constant,time);
% 地球
eTrue = Earth(time,constant.sunMu);
eTrue.getEphem();
% 探査機 
% 探査機自身が持つ推定値(観測とともに更新)
scEst = Spacecraft(time,constant.sunMu,"spacecraft");
scEst.state=sc.state0;
scEst.clockError = sc_est.X_hat(1);
scEst.resClockError(1) = error.initialClock - scEst.clockError;
% 真値(誤差を入れる)
scTrue = Spacecraft(time,constant.sunMu,"spacecraft");
scTrue.state = scEst.state + [error.scPos0; error.scVel0];
% 地上局が推定している値
scEstGs = Spacecraft(time,constant.sunMu,"spacecraft");
scEstGs.state = scEst.state;



%% 2. calculate true orbit of spacecraft
% 真値の計算
scTrue.calcOrbitTwoBody(time.simDt,error.dynamics)
%% 3. calculate initial guess of spacecraft orbit (without observation)
scEstGs.calcOrbitTwoBody(time.simDt,0)
scIniGuess = scEstGs.state;

%% ここからは，time stepごとの計算をしていく
gsTransNum = 0;
scReceiveNum = 0;
gsReceiveNum = 0;
% 地上局が2wayを観測できるかどうか
gs2way = 0; % 0: 観測前, 1:観測後
% 探査機が2wayを観測できるようになる時間(初期化)
time.sc2wayGet = time.list(length(time.list));

for i = 1:length(time.list)-1
    disp(i)
    % 探索一回につき観測一回
    if mod(i,time.obsStep)==1 || time.obsStep ==1
        gsTransNum = gsTransNum + 1;
        %% 4.ground station Send a signal around the estimated orbit to the spacecraft according to the search pattern.
        % 地上局が推定している探査機の軌道から，目標方向(stateEstOpn)と到達時刻(tEstOpn)を計算
        opnUpEstTemp = GroundStation.calcTarget(time.list(i),gsTrue.state(:,i),eTrue.state(:,i),scEstGs,time,constant);
        %[new]実際に時刻time.list(i)で向けるべきだった方向を計算する
        opnUpTrueTemp = GroundStation.calcTarget(time.list(i),gsTrue.state(:,i),eTrue.state(:,i),scTrue,time,constant);
        % 推定値周りに探索した時に，探索にかかる時間を求め，tTrans, stateTrans,
        % azmTrans, elvTransを求める
        [gsTrue,eTrue] = GroundStation.search(i,gsTrue,eTrue,gs,time,constant,opnUpEstTemp,opnUpTrueTemp);
        % レーザーを照射した時刻gsTrue.tTrans, とその時刻での
        %% 5. calculate observed value using the send signal
        % tTransにstateTransから照射された光が宇宙機に到達する時刻を求め, その時刻の宇宙機の状態量と観測量を計算する
        % two-way取れるかで場合わけ
        if gs2way == 0
            %one wayの観測
            [scTrue,gsTrue] = Spacecraft.calcObservedValue(scTrue,gsTrue,eTrue,i,constant,time,gs,sc,error);
            scRec.ltdObs(gsTransNum) = scTrue.lengthObserved(i) /constant.lightSpeed;
        else
            %one wayの観測. ただし，地上局は2wayが取れているので，この送信に対して，探査機は2wayが取れる
            [scTrue,gsTrue] = Spacecraft.calcObservedValue2way(scTrue,gsTrue,eTrue,i,constant,time,gs,sc,error,scTrans,gsReceiveNum);    
            scRec.rtltObserved(gsTransNum) = scTrue.rtltObserved(i);
            %地上局がdownlinkを受信する時刻の状態量(DownlinkReceive)uplink信号に含まれる
            scRec.xveDr(:,gsTransNum) = scTrans.eReceive(:,gsReceiveNum); 
            scRec.xvgDr(:,gsTransNum) = scTrans.gsReceive(:,gsReceiveNum);
            scRec.duration(gsTransNum) = scTrue.duration(i);
            % 初めて，探査機が2wayを観測できる時刻を求める
            if gsReceiveNum == 1
                time.sc2wayGet = scTrue.tReceive(i);
            end
        end
        % 1wayでも2wayでも共通のもの
        scRec.t(gsTransNum) = scTrue.tReceive(i); % 受信時刻
        scRec.azmObs(gsTransNum) = scTrue.azmObserved(i);
        scRec.elvObs(gsTransNum) = scTrue.elvObserved(i);
        scRec.acelObs(:,gsTransNum) = scTrue.accelObserved(:,i);
        scRec.tTrans(gsTransNum) = gsTrue.tTrans(i); %光が送信された時刻
        scRec.xve(:,gsTransNum) = eTrue.stateTrans(:,i);
        scRec.xvg(:,gsTransNum) = gsTrue.stateTrans(:,i); 
        scRec.azmTrue(gsTransNum) = scTrue.azmTrue(i);
        scRec.elvTrue(gsTransNum) = scTrue.elvTrue(i);
        scRec.angleError(gsTransNum) = scTrue.angleError(i);
    end
        

    %% 6.estimate spacecraft orbit using observed value by EKF(探査機が自身の軌道を推定する)
    % 2wayが取れる前後で場合わけ
    if time.list(i+1) < time.sc2wayGet
        %%  次の時刻までに観測がなかった時 → リファレンス，誤差共分散を更新
        if time.list(i+1) < scRec.t(scReceiveNum +1)
           sc_est.Dt = time.simDt;
           % リファレンスと誤差共分散行列を更新
           [sc_est.X_hat, sc_est.P] = Spacecraft.updateState(sc_est.X_hat,sc_est.P, sc_est.Dt,scEst.mu);
        else
        % 次の時刻までに観測があった時→dt1とdt2に分けて更新
            sc_est.Dt1 = scRec.t(scReceiveNum +1) - time.list(i);
            % dt1までリファレンスと誤差共分散行列を更新
            [sc_est.X_bar, sc_est.P_bar] = Spacecraft.updateState(sc_est.X_hat,sc_est.P, sc_est.Dt1,scEst.mu);
            % 観測ベクトル
            sc_est.Y = [scRec.ltdObs(scReceiveNum +1);  scRec.azmObs(scReceiveNum +1); scRec.elvObs(scReceiveNum +1); scRec.acelObs(:,scReceiveNum +1)];
            % リファレンスの状態量の時の観測量
            sc_est.Y_bar = Spacecraft.calcG(sc_est.X_bar,scRec.xve(:,scReceiveNum +1),scRec.xvg(:,scReceiveNum +1),constant);
            sc_est.y     = sc_est.Y - sc_est.Y_bar;
            sc_est.H_childa = Spacecraft.delGdelX(sc_est.X_bar,scRec.xve(:,scReceiveNum +1),scRec.xvg(:,scReceiveNum +1),constant);
            % 観測誤差共分散行列の更新
            sc_est.R(2,2) = scRec.angleError(scReceiveNum +1)^2;
            sc_est.R(3,3) = sc_est.R(2,2);
            % カルマンゲインの計算と推定値の更新
            sc_est.K = sc_est.P_bar * sc_est.H_childa.'/(sc_est.H_childa*sc_est.P_bar*sc_est.H_childa.' + sc_est.R);
            sc_est.X_hat = sc_est.X_bar + sc_est.K*sc_est.y;
            sc_est.P = (eye(7) - sc_est.K*sc_est.H_childa)*sc_est.P_bar;
            % dt2の間は普通にリファレンスと誤差共分散行列を更新
            sc_est.Dt2 = time.list(i+1) - scRec.t(scReceiveNum +1);
            [sc_est.X_hat, sc_est.P] = Spacecraft.updateState(sc_est.X_hat,sc_est.P, sc_est.Dt2,scEst.mu);  
        end
        % 探査機の軌道の推定値を記録する
        if i == length(time.list)
        else
        scEst.resClockError(i+1)= error.initialClock - sc_est.X_hat(1);
        scEst.state(:,i+1) = sc_est.X_hat(2:7);
        sc_est.P_list(:,:,i+1) = sc_est.P;
        end
        %2wayの時に使うものだけ抽出
        sc_est.X_hat2w = sc_est.X_hat(2:7);
        sc_est.P2w     = sc_est.P(2:7,2:7);
    else % 2wayが取れた後
        % 次の時刻まで観測がなかった場合
        if time.list(i+1) < scRec.t(scReceiveNum +1)
           sc_est.Dt = time.simDt;
           % リファレンスと誤差共分散行列を更新
           [sc_est.X_hat2w, sc_est.P2w] = Spacecraft.updateState2(sc_est.X_hat2w,sc_est.P2w, sc_est.Dt,scEst.mu);
        else
        % 次の時刻までに観測があった時→dt1とdt2に分けて更新
            sc_est.Dt1 = scRec.t(scReceiveNum +1) - time.list(i);
            % dt1までリファレンスと誤差共分散行列を更新
            [sc_est.X_bar2w, sc_est.P_bar2w] = Spacecraft.updateState2(sc_est.X_hat2w,sc_est.P2w, sc_est.Dt1,scEst.mu);
            % 観測ベクトル
            sc_est.Y2w = [scRec.rtltObserved(scReceiveNum +1);  scRec.azmObs(scReceiveNum +1); scRec.elvObs(scReceiveNum +1); scRec.acelObs(:,scReceiveNum +1)];
            % リファレンスの状態量の時の観測量
            sc_est.Y_bar2w = Spacecraft.calcG2w(sc_est.X_bar2w,scRec.xve(:,scReceiveNum +1),scRec.xvg(:,scReceiveNum +1),scRec.xveDr(:,scReceiveNum +1),scRec.xvgDr(:,scReceiveNum +1),scTrue.duration(scReceiveNum +1),constant);
            sc_est.y2w     = sc_est.Y2w - sc_est.Y_bar2w;
            sc_est.H_childa2w = Spacecraft.delGdelX2w(sc_est.X_bar2w,scRec.xve(:,scReceiveNum +1),scRec.xvg(:,scReceiveNum +1),scRec.xveDr(:,scReceiveNum +1),scRec.xvgDr(:,scReceiveNum +1),scTrue.duration(scReceiveNum +1),constant);
            % 観測誤差共分散行列の更新
            sc_est.R(2,2) = scRec.angleError(scReceiveNum +1)^2;
            sc_est.R(3,3) = sc_est.R(2,2);          
            % カルマンゲインの計算と推定値の更新
            sc_est.K2w = sc_est.P_bar2w * sc_est.H_childa2w.'/(sc_est.H_childa2w*sc_est.P_bar2w*sc_est.H_childa2w.' + sc_est.R);
            sc_est.X_hat2w = sc_est.X_bar2w + sc_est.K2w*sc_est.y2w;
            sc_est.P2w = (eye(6) - sc_est.K2w*sc_est.H_childa2w)*sc_est.P_bar2w;
            % 時計誤差に換算する
            if i ==  length(time.list)
            else
            scEst.resClockError(i+1)= norm(sc_est.X_hat2w(1:3) - scRec.xve(1:3,scReceiveNum +1) - scRec.xvg(1:3,scReceiveNum +1))/constant.lightSpeed ...
                + scRec.tTrans(scReceiveNum +1) - scRec.t(scReceiveNum +1);
            end
            % dt2の間は普通にリファレンスと誤差共分散行列を更新
            sc_est.Dt2 = time.list(i+1) - scRec.t(scReceiveNum +1);
            [sc_est.X_hat2w, sc_est.P2w] = Spacecraft.updateState2(sc_est.X_hat2w,sc_est.P2w, sc_est.Dt2,scEst.mu);  
        end
    end
   % 探査機の軌道の推定値を記録する
   if i == length(time.list)
   else
   scEst.state(:,i+1) = sc_est.X_hat2w;
   sc_est.P_list(1,1,i+1) = sc_est.P_list(1,1,i);
   sc_est.P_list(2:7,2:7,i+1) = sc_est.P2w;
   end
   
%%  7.  ダウンリンク推定方向の計算. (PAAの計算に次時刻の推定値が必要なので，ここで計算している)
    if i ==1
        j = i;
        scEst.calcDownDirection(time.list(j),gsTrue,eTrue,scEst.state(:,j),time,constant,j);
        scTrue.calcDownDirection(time.list(j),gsTrue,eTrue,scTrue.state(:,j),time,constant,j);
        gsTrue.calcObservedValue(scTrue,scEst,j,constant,time,gs,sc,error);
    end
    j = i+1;
    scEst.calcDownDirection(time.list(j),gsTrue,eTrue,scEst.state(:,j),time,constant,j);
    scTrue.calcDownDirection(time.list(j),gsTrue,eTrue,scTrue.state(:,j),time,constant,j);
    gsTrue.calcObservedValue(scTrue,scEst,j,constant,time,gs,sc,error);
      
    % PAAの計算　(今の時刻次の時刻の間にuplink観測があった場合，そのuplinkと現時刻の地上曲推定方向を求める)
    if time.list(i+1) > scRec.t(scReceiveNum+1)
       scTrans.ID(scReceiveNum+1) = i+1;
       scTrans.t(scReceiveNum+1) = time.list(i+1);
       scTrans.stateTrue(:,scReceiveNum+1) = scTrue.state(:,i+1);
       scTrans.stateEstGs(:,scReceiveNum+1) = scEstGs.state(:,i+1);
       scTrans.tReceive(scReceiveNum+1) = scTrue.tDown(i+1);
       scTrans.eReceive(:,scReceiveNum+1) = scTrue.eDown(:,i+1);
       scTrans.gsReceive(:,scReceiveNum+1) = scTrue.gsDown(:,i+1);
       % 地上局で観測される観測量
       scTrans.azmObserved(scReceiveNum+1) = gsTrue.azmObserved(i+1);
       scTrans.elvObserved(scReceiveNum+1) = gsTrue.elvObserved(i+1);
       scTrans.angleError(scReceiveNum+1) = gsTrue.angleError(i+1);
       scReceiveNum = scReceiveNum + 1;
   else    
   end
   
   
   %% 地上局が探査機の軌道を観測値から推定する
   if  time.list(i+1) < scRec.t(1) % scTransが定義できないため場合わけ
       % 時刻time.list(i+1)の推定値を更新する
       sc_estGs.Dt2 = time.simDt;
       % リファレンスと誤差共分散行列を更新
       %% 時計誤差を状態量に含まないので，updateStateを変更する必要あり
       [sc_estGs.X_hat2, sc_estGs.P2] = Spacecraft.updateState2(sc_estGs.X_hat2,sc_estGs.P2, sc_estGs.Dt2,scEst.mu);
   elseif time.list(i+1) < scTrans.tReceive(gsReceiveNum+1) %次の時刻まで観測がない場合はそのまま推定値を更新する
       % 時刻time.list(i+1)の推定値を更新する
       sc_estGs.Dt2 = time.simDt;
       % リファレンスと誤差共分散行列を更新
       [sc_estGs.X_hat2, sc_estGs.P2] = Spacecraft.updateState2(sc_estGs.X_hat2,sc_estGs.P2, sc_estGs.Dt2,scEst.mu);
   else % 次の時刻までに観測があった時．
        %  ダウンリンクを送信した時刻=scTrans.t(gsReceiveNum+1)の状態量を求める．→　time.list(i+1)まで伝搬して，次時刻の状態量を推定する
        if gsReceiveNum == 0 
            sc_estGs.X_bar1 = scTrans.stateEstGs(:,1);
            sc_estGs.P_bar1  = sc_estGs.P1;
            gs2way = 1; %以降地上局が2-wayを取れるようになった
        else
            sc_estGs.Dt1 = scTrans.t(gsReceiveNum+1) - scTrans.t(gsReceiveNum);
            [sc_estGs.X_bar1, sc_estGs.P_bar1] = Spacecraft.updateState2(sc_estGs.X_hat1,sc_estGs.P1, sc_estGs.Dt1,scEst.mu);
        end
        % 観測量 往復にかかった時間，downlinkされた方向の方位角，仰角
        sc_estGs.RTLT = (scTrans.tReceive(gsReceiveNum+1) - scTrans.t(gsReceiveNum+1)) + (scRec.t(gsReceiveNum+1) - scRec.tTrans(gsReceiveNum+1));
        sc_estGs.Y = [sc_estGs.RTLT; scTrans.azmObserved(gsReceiveNum+1);scTrans.elvObserved(gsReceiveNum+1)];
        % リファレンスの状態量の時の観測量の計算
        sc_estGs.Y_bar = Spacecraft.calcG2(sc_estGs.X_bar1,scRec.xve(:,gsReceiveNum+1),scRec.xvg(:,gsReceiveNum+1),scTrans.eReceive(:,gsReceiveNum+1),scTrans.gsReceive(:,gsReceiveNum+1),scTrans.t(gsReceiveNum+1) - scRec.t(gsReceiveNum+1), constant);
        sc_estGs.y     = sc_estGs.Y - sc_estGs.Y_bar;
        sc_estGs.H_childa = Spacecraft.delGdelX2(sc_estGs.X_bar1,scRec.xve(:,gsReceiveNum+1),scRec.xvg(:,gsReceiveNum+1),scTrans.eReceive(:,gsReceiveNum+1),scTrans.gsReceive(:,gsReceiveNum+1),scTrans.t(gsReceiveNum+1) - scRec.t(gsReceiveNum+1), constant);
        % 観測誤差共分散行列の計算
        sc_estGs.R1(2,2) = scTrans.angleError(gsReceiveNum+1)^2;
        sc_estGs.R1(3,3) = sc_estGs.R1(2,2);
        % カルマンゲインの計算
        sc_estGs.K = sc_estGs.P_bar1 * sc_estGs.H_childa.'/(sc_estGs.H_childa*sc_estGs.P_bar1*sc_estGs.H_childa.' + sc_estGs.R1);
        % 状態量の更新
        sc_estGs.X_hat1 = sc_estGs.X_bar1 + sc_estGs.K*sc_estGs.y;
        sc_estGs.P1 = (eye(6) - sc_estGs.K*sc_estGs.H_childa)*sc_estGs.P_bar1;

        % ダウンリンクを送信した時刻の状態量から，ダウンリンクを受信した時刻の状態量を求める
        sc_estGs.X_hat2 = sc_estGs.X_hat1;
        sc_estGs.P2     = sc_estGs.P1;
        for l = 1:round( i+1 - scTrans.ID(gsReceiveNum+1) )
            [sc_estGs.X_hat2, sc_estGs.P2] = Spacecraft.updateState2(sc_estGs.X_hat2,sc_estGs.P2, time.simDt ,scEst.mu);
        end
        gsReceiveNum = gsReceiveNum + 1;
   end
%    探査機の軌道の推定値を記録する
   if i == length(time.list)
   else
   scEstGs.state(:,i+1) = sc_estGs.X_hat2;
   sc_est.P2_list(:,:,i+1) = sc_estGs.P2;
   end
   
   
end

% 時計誤差の時間履歴
figure(1)
title('clock error of spacecraft')
hold on
plot((time.list-time.list(1))/60/60, scEst.resClockError)
xlabel('time [h]')
ylabel('clock error [s]')
hold off

% 軌道決定精度をplotする
% 探査機位置誤差の時間履歴
figure(2)
tiledlayout(3,1)
nexttile
title('position error of estimated value')
hold on
plot((time.list-time.list(1))/60/60, scEst.state(1,:) - scTrue.state(1,:))
plot((time.list-time.list(1))/60/60, scEst.state(2,:) - scTrue.state(2,:))
plot((time.list-time.list(1))/60/60, scEst.state(3,:) - scTrue.state(3,:))
plot((time.list-time.list(1))/60/60,...
    ( (scEst.state(1,:) - scTrue.state(1,:)).^2 + (scEst.state(2,:) - scTrue.state(2,:)).^2 + (scEst.state(3,:) - scTrue.state(3,:)).^2).^0.5)
xlabel('time [h]')
ylabel('position error [km]')
legend('x', 'y', 'z','length')
% ylim([-1000 1000])
hold off
nexttile
title('velocity error of estimated value')
hold on
plot((time.list-time.list(1))/60/60, scEst.state(4,:) - scTrue.state(4,:))
plot((time.list-time.list(1))/60/60, scEst.state(5,:) - scTrue.state(5,:))
plot((time.list-time.list(1))/60/60, scEst.state(6,:) - scTrue.state(6,:))
plot((time.list-time.list(1))/60/60,...
    ( (scEst.state(4,:) - scTrue.state(4,:)).^2 + (scEst.state(5,:) - scTrue.state(5,:)).^2 + (scEst.state(6,:) - scTrue.state(6,:)).^2).^0.5)
xlabel('time [h]')
ylabel('position error [km]')
legend('x', 'y', 'z','speed')
% ylim([-1000 1000])
hold off
nexttile
title('downlink direction error')
hold on
plot((time.list(1:length(time.list)) - time.list(1))/60/60, scEst.azmDown - scTrue.azmDown)
plot((time.list(1:length(time.list)) - time.list(1))/60/60, scEst.elvDown - scTrue.elvDown)
plot((time.list(1:length(time.list)) - time.list(1))/60/60, ((scEst.azmDown - scTrue.azmDown).^2 + (scEst.elvDown - scTrue.elvDown).^2).^0.5);
hold off
xlabel('time [h]')
ylabel('angle [rad]')
legend('azimuth', 'elevation', 'angle')

% 地上局の推定値
figure(3)
tiledlayout(2,1)
nexttile
title('position error of estimated value')
hold on
plot((time.list-time.list(1))/60/60, scEstGs.state(1,:) - scTrue.state(1,:))
plot((time.list-time.list(1))/60/60, scEstGs.state(2,:) - scTrue.state(2,:))
plot((time.list-time.list(1))/60/60, scEstGs.state(3,:) - scTrue.state(3,:))
plot((time.list-time.list(1))/60/60,...
    ( (scEstGs.state(1,:) - scTrue.state(1,:)).^2 + (scEstGs.state(2,:) - scTrue.state(2,:)).^2 + (scEstGs.state(3,:) - scTrue.state(3,:)).^2).^0.5)
xlabel('time [h]')
ylabel('position error [km]')
legend('x', 'y', 'z','length')
nexttile
title('velocity error of estimated value')
hold on
plot((time.list-time.list(1))/60/60, scEstGs.state(4,:) - scTrue.state(4,:))
plot((time.list-time.list(1))/60/60, scEstGs.state(5,:) - scTrue.state(5,:))
plot((time.list-time.list(1))/60/60, scEstGs.state(6,:) - scTrue.state(6,:))
plot((time.list-time.list(1))/60/60,...
    ( (scEstGs.state(4,:) - scTrue.state(4,:)).^2 + (scEstGs.state(5,:) - scTrue.state(5,:)).^2 + (scEstGs.state(6,:) - scTrue.state(6,:)).^2).^0.5)
xlabel('time [h]')
ylabel('velocity error [km]')
legend('x', 'y', 'z','speed')

figure(4)
title('estemated/true position of spacecraft')
hold on
plot3(scTrue.state(1,:),scTrue.state(2,:),scTrue.state(3,:))
plot3(scEst.state(1,:),scEst.state(2,:),scEst.state(3,:))
plot3(scEstGs.state(1,:),scEstGs.state(2,:),scEstGs.state(3,:))
plot3(scIniGuess(1,:),scIniGuess(2,:),scIniGuess(3,:))
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
legend('true', 'estimated(sc)', 'estimated(gs)','initial guess')
hold off




