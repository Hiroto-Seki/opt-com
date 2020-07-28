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
clear all; close all; clc
% add path to SPICE
addpath(genpath('~/Documents/Matlab/SPICE'));
% SPICEのKernel(天体情報)を読み込む
spice_loadkernels();
% 取得した天体情報 + alphaを利用しやすいように構造体へまとめる
SSD = spice_setparams();
% 乱数
rng('default');
rng(3)

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
for i = 1:length(time.list)-1
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
        [scTrue,gsTrue] = Spacecraft.calcObservedValue(scTrue,gsTrue,eTrue,i,constant,time,gs,sc,error);
        scRec.t(gsTransNum) = scTrue.tReceive(i); % 受信時刻
        scRec.ltdObs(gsTransNum) = scTrue.lengthObserved(i) /constant.lightSpeed;
        scRec.azmObs(gsTransNum) = scTrue.azmObserved(i);
        scRec.elvObs(gsTransNum) = scTrue.elvObserved(i);
        scRec.tTrans(gsTransNum) = gsTrue.tTrans(i); %光が送信された時刻
        scRec.xve(:,gsTransNum) = gsTrue.stateTrans(:,i);
        scRec.xvg(:,gsTransNum) = eTrue.stateTrans(:,i); 
        scRec.azmTrue(gsTransNum) = scTrue.azmTrue(i);
        scRec.elvTrue(gsTransNum) = scTrue.elvTrue(i);
    end
        

    %% 6.estimate spacecraft orbit using observed value by EKF(探査機が自身の軌道を推定する)
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
        sc_est.Y = [scRec.ltdObs(scReceiveNum +1);  scRec.azmObs(scReceiveNum +1); scRec.elvObs(scReceiveNum +1)];
        % リファレンスの状態量の時の観測量
        sc_est.Y_bar = Spacecraft.calcG(sc_est.X_bar,scRec.xve(:,scReceiveNum +1),scRec.xvg(:,scReceiveNum +1),constant);
        sc_est.y     = sc_est.Y - sc_est.Y_bar;
        sc_est.H_childa = Spacecraft.delGdelX(sc_est.X_bar,scRec.xve(:,scReceiveNum +1),scRec.xvg(:,scReceiveNum +1),constant);
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
   scEst.clockError(i+1)= sc_est.X_hat(1);
   scEst.state(:,i+1) = sc_est.X_hat(2:7);
   sc_est.P_list(:,:,i+1) = sc_est.P;
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
            sc_estGs.P1_bar  = sc_estGs.P1;
        else
            sc_estGs.Dt1 = scTrans.t(gsReceiveNum+1) - scTrans.t(gsReceiveNum);
            [sc_estGs.X_bar1, sc_estGs.P1_bar] = Spacecraft.updateState2(sc_estGs.X_hat1,sc_estGs.P1, sc_estGs.Dt1,scEst.mu);
        end
        % 観測量 往復にかかった時間，downlinkされた方向の方位角，仰角
        sc_estGs.RTLT = (scTrans.tReceive(gsReceiveNum+1) - scTrans.t(gsReceiveNum+1)) + (scRec.t(gsReceiveNum+1) - scRec.tTrans(gsReceiveNum+1));
        sc_estGs.Y = [sc_estGs.RTLT; scTrans.azmObserved(gsReceiveNum+1);scTrans.elvObserved(gsReceiveNum+1)];
        % リファレンスの状態量の時の観測量の計算
        sc_estGs.Y_bar = Spacecraft.calcG2(sc_estGs.X_bar1,scRec.xve(:,gsReceiveNum+1),scRec.xvg(:,gsReceiveNum+1),scTrans.eReceive(:,gsReceiveNum+1),scTrans.gsReceive(:,gsReceiveNum+1),scTrans.t(gsReceiveNum+1) - scRec.t(gsReceiveNum+1), constant);
        sc_estGs.y     = sc_estGs.Y - sc_estGs.Y_bar;
        sc_estGs.H_childa = Spacecraft.delGdelX2(sc_estGs.X_bar1,scRec.xve(:,gsReceiveNum+1),scRec.xve(:,gsReceiveNum+1),scTrans.eReceive(:,gsReceiveNum+1),scTrans.gsReceive(:,gsReceiveNum+1),scTrans.t(gsReceiveNum+1) - scRec.t(gsReceiveNum+1), constant);
        % カルマンゲインの計算
        sc_estGs.K = sc_estGs.P1_bar * sc_estGs.H_childa.'/(sc_estGs.H_childa*sc_estGs.P1_bar*sc_estGs.H_childa.' + sc_estGs.R1);
        % 状態量の更新
        sc_estGs.X_hat1 = sc_estGs.X_bar1 + sc_estGs.K*sc_estGs.y;
        sc_estGs.P1 = (eye(6) - sc_estGs.K*sc_estGs.H_childa)*sc_estGs.P1_bar;

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
plot(time.list-time.list(1), error.initialClock - scEst.clockError)
xlabel('time [s]')
ylabel('clock error [s]')
hold off

% 軌道決定精度をplotする
% 探査機位置誤差の時間履歴
figure(2)
tiledlayout(2,1)
nexttile
title('position error of estimated value')
hold on
plot(time.list-time.list(1), scEst.state(1,:) - scTrue.state(1,:))
plot(time.list-time.list(1), scEst.state(2,:) - scTrue.state(2,:))
plot(time.list-time.list(1), scEst.state(3,:) - scTrue.state(3,:))
plot(time.list-time.list(1),...
    ( (scEst.state(1,:) - scTrue.state(1,:)).^2 + (scEst.state(2,:) - scTrue.state(2,:)).^2 + (scEst.state(3,:) - scTrue.state(3,:)).^2).^0.5)
plot(time.list-time.list(1), 450 * ones(length(time.list),1),'g')
plot(time.list-time.list(1), -450 * ones(length(time.list),1),'g')
xlabel('time [s]')
ylabel('position error [km]')
legend('x', 'y', 'z','length','450km')
% ylim([-1000 1000])
hold off
nexttile
title('downlink direction error')
hold on
plot(time.list(1:length(time.list)) - time.list(1), scEst.azmDown - scTrue.azmDown)
plot(time.list(1:length(time.list))- time.list(1), scEst.elvDown - scTrue.elvDown)
plot(time.list(1:length(time.list)) - time.list(1), ((scEst.azmDown - scTrue.azmDown).^2 + (scEst.elvDown - scTrue.elvDown).^2).^0.5);
plot(time.list(1:length(time.list))-time.list(1),   3e-7* ones(length(time.list),1),'g');
plot(time.list(1:length(time.list))-time.list(1),   -3e-7* ones(length(time.list),1)),'g';
hold off
xlabel('time [s]')
ylabel('angle [rad]')
legend('azimuth', 'elevation', 'angle','requirement')

% 地上局の推定値
figure(3)
title('position error of estimated value')
hold on
plot(time.list-time.list(1), scEstGs.state(1,:) - scTrue.state(1,:))
plot(time.list-time.list(1), scEstGs.state(2,:) - scTrue.state(2,:))
plot(time.list-time.list(1), scEstGs.state(3,:) - scTrue.state(3,:))
plot(time.list-time.list(1),...
    ( (scEstGs.state(1,:) - scTrue.state(1,:)).^2 + (scEstGs.state(2,:) - scTrue.state(2,:)).^2 + (scEstGs.state(3,:) - scTrue.state(3,:)).^2).^0.5)
xlabel('time [s]')
ylabel('position error [km]')
legend('x', 'y', 'z','length')


figure(4)
title('estemated/true position of spacecraft')
hold on
plot3(scTrue.state(1,:),scTrue.state(2,:),scTrue.state(3,:))
plot3(scEst.state(1,:),scEst.state(2,:),scEst.state(3,:))
plot3(scEstGs.state(1,:),scEstGs.state(2,:),scEstGs.state(3,:))
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
legend('true', 'estimated(sc)', 'estimated(gs)')
hold off




