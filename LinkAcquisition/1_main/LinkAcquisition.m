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
rng(1)

%% setting parameter
[constant,time,error,gs,sc,X_hat,P,P_list,R] = setparam(SSD);

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
scEst.clockError = X_hat(1);
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

%% ここからは，time stepごとの計算をしていく
transNum = 0;
receiveNum = 0;
for i = 1:length(time.list)-1
    % 探索一回につき観測一回
    if mod(i,time.obsStep)==1 || time.obsStep ==1
        transNum = transNum + 1;
        %% 4.ground station Send a signal around the estimated orbit to the spacecraft according to the search pattern.
        % 地上局が推定している探査機の軌道から，目標方向(stateEstOpn)と到達時刻(tEstOpn)を計算
        opnEstTemp = GroundStation.calcTarget(time.list(i),gsTrue.state(:,i),eTrue.state(:,i),scEstGs,time,constant);
        %[new]実際に時刻time.list(i)で向けるべきだった方向を計算する
        opnTrueTemp = GroundStation.calcTarget(time.list(i),gsTrue.state(:,i),eTrue.state(:,i),scTrue,time,constant);
        % 推定値周りに探索した時に，探索にかかる時間を求め，tTrans, stateTrans,
        % azmTrans, elvTransを求める
        [gsTrue,eTrue] = GroundStation.search(i,gsTrue,eTrue,gs,time,constant,opnEstTemp,opnTrueTemp);
        % レーザーを照射した時刻gsTrue.tTrans, とその時刻での
        %% 5. calculate observed value using the send signal
        % tTransにstateTransから照射された光が宇宙機に到達する時刻を求め, その時刻の宇宙機の状態量と観測量を計算する
        [scTrue,gsTrue] = Spacecraft.calcObservedValue(scTrue,gsTrue,eTrue,i,constant,time,gs,sc,error);
        obs.t(transNum) = scTrue.tReceive(i);
        obs.ltd(transNum) = scTrue.lengthObserved(i) /constant.lightSpeed;
        obs.azm(transNum) = scTrue.azmObserved(i);
        obs.elv(transNum) = scTrue.elvObserved(i);
        obs.xve(:,transNum) = gsTrue.stateTrans(:,i);
        obs.xvg(:,transNum) = eTrue.stateTrans(:,i); 
    end

    %% 6.estimate spacecraft orbit using observed value by EKF
    %%  次の時刻までに観測がなかった時 → リファレンス，誤差共分散を更新
    if time.list(i+1) < obs.t(receiveNum +1)
       dt = time.simDt;
       % リファレンスと誤差共分散行列を更新
       [X_hat, P] = Spacecraft.updateState(X_hat,P, dt,scEst.mu);
    else
    % 次の時刻までに観測があった時→dt1とdt2に分けて更新
        dt1 = obs.t(receiveNum +1) - time.list(i);
        % dt1までリファレンスと誤差共分散行列を更新
        [X_bar, P_bar] = Spacecraft.updateState(X_hat,P, dt1,scEst.mu);
%         % 観測ベクトル
        Y = [obs.ltd(receiveNum +1);  obs.azm(receiveNum +1); obs.elv(receiveNum +1)];
        % リファレンスの状態量の時の観測量
        Y_bar = Spacecraft.calcG(X_bar,obs.xve(:,receiveNum +1),obs.xvg(:,receiveNum +1),constant);
        y     = Y - Y_bar;
        H_childa = Spacecraft.delGdelX(X_bar,obs.xve(:,receiveNum +1),obs.xvg(:,receiveNum +1),constant);
        % カルマンゲインの計算と推定値の更新
        K = P_bar * H_childa.'/(H_childa*P_bar*H_childa.' + R);
        X_hat = X_bar + K*y;
        P = (eye(7) - K*H_childa)*P_bar;
        % dt2の間は普通にリファレンスと誤差共分散行列を更新
        dt2 = time.list(i+1) - obs.t(receiveNum +1);
        [X_hat, P] = Spacecraft.updateState(X_hat,P, dt2,scEst.mu);
        receiveNum = receiveNum + 1;  
    end
    
    
    
    % 記録する
   if i == length(time.list)
   else
   scEst.clockError(i+1)= X_hat(1);
   scEst.state(:,i+1) = X_hat(2:7);
   P_list(:,:,i+1) = P;
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
plot(time.list-time.list(1), 450 * ones(length(time.list),1))
xlabel('time [s]')
ylabel('position error [km]')
legend('x', 'y', 'z','length','450km')
% ylim([-1000 1000])
hold off
nexttile
title('position error of initail guess')
hold on
plot(time.list-time.list(1), scEstGs.state(1,:) - scTrue.state(1,:))
plot(time.list-time.list(1), scEstGs.state(2,:) - scTrue.state(2,:))
plot(time.list-time.list(1), scEstGs.state(3,:) - scTrue.state(3,:))
xlabel('time [s]')
ylabel('position error [km]')
legend('x', 'y', 'z')
% ylim([-1000 1000])
hold off


figure(3)
title('estemated/true position of spacecraft')
hold on
plot3(scTrue.state(1,:),scTrue.state(2,:),scTrue.state(3,:))
plot3(scEst.state(1,:),scEst.state(2,:),scEst.state(3,:))
plot3(scEstGs.state(1,:),scEstGs.state(2,:),scEstGs.state(3,:))
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
legend('true', 'estimated', 'initial guess')
hold off


%% 地球との相対位置を求める
figure(4)
title('estemated/true position of spacecraft')
hold on
plot3(scTrue.state(1,:)-eTrue.state(1,:)-gsTrue.state(1,:),scTrue.state(2,:)-eTrue.state(2,:)-gsTrue.state(2,:),scTrue.state(3,:)-eTrue.state(3,:)-gsTrue.state(3,:))
plot3(scEst.state(1,:)-eTrue.state(1,:)-gsTrue.state(1,:),scEst.state(2,:)-eTrue.state(2,:)-gsTrue.state(2,:),scEst.state(3,:)-eTrue.state(3,:)-gsTrue.state(3,:))
plot3(scEstGs.state(1,:)-eTrue.state(1,:)-gsTrue.state(1,:),scEstGs.state(2,:)-eTrue.state(2,:)-gsTrue.state(2,:),scEstGs.state(3,:)-eTrue.state(3,:)-gsTrue.state(3,:))
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
legend('true', 'estimated', 'initial guess')
hold off





% 可視化したいもの
% 地上局送信方向の履歴と真値，推定値の履歴
% 軌道決定精度
% 通信可否
% 観測量と観測精度







%% 線型近侍した目標値と実際の観測地点の差
% search中に出てくるやつ


