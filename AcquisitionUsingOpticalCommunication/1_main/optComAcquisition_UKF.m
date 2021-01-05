%{
誤差共分散が小さい時に，EKFではうまく収束しなかったので，EKFの代わりにUKFを使ってみる
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
rng(2)

%% 1.setting parameter and initial state
[constant,time,error,gs,sc,gsTrue,earth,scTrue,scEstByScUkf,scEstByGsUkf,~,ukf] = setparam(SSD);

%% 2. calculate true orbit of spacecraft
% 真値の計算
scTrue.calcOrbitTwoBody(constant.sunMu, error.dynamics)
%% 3. calculate initial guess of spacecraft orbit (without observation)
scEstByScUkf.calcOrbitTwoBody(constant.sunMu,0)
scEstByGsUkf.calcOrbitTwoBody(constant.sunMu,0)

%% ここからは，time stepごとの計算をしていく
for i = 1:length(time.list)-1
    %% 4 地上局からの送信
    % 探索を開始するtime.stepの計算
    if  i == 1  || i == time.lastSearch + time.obsStep
         % 探索範囲の設定. 今回の探索にかかる時間=次の探索が始まる時間を求める．
        [time,scEstByScUkf.R1wSc,scEstByScUkf.R2wSc] = GroundStation.setSearchArea(time,gs,SSD,scEstByGsUkf.P,scEstByScUkf.R1wSc,scEstByScUkf.R2wSc,error);
        % 地上局が推定している探査機の軌道から目標方向と到達時刻を計算する
        [gsTrue.opnEstTempT_ut,gsTrue.opnEstTempState_ut] ...
            = GroundStation.calcTarget(time.list(i),gsTrue.state(:,i),earth.state(:,i),scEstByGsUkf.state(:,i),scEstByGsUkf,time,constant,"estimated value");
        % 真の軌道から，本来の向けるべき方向を計算する
        [gsTrue.opnTrueTempT_ut,gsTrue.opnTrueTempState_ut] ...
            = GroundStation.calcTarget(time.list(i),gsTrue.state(:,i),earth.state(:,i),scTrue.state(:,i),scTrue,time,constant,"true value");
        % 推定値周りに探索して，送信時刻,送信方向,送信時の地上局の位置・速度を求める
        [gsTrue,earth] = gsTrue.search(i,earth,gs,time,constant,error);
        % 宇宙機に届く時刻と，宇宙機が受信する内容を求める
        [scTrue,gsTrue] = scTrue.receiveUplink(gsTrue,earth,constant,time);
        % 初めて，探査機が2wayを観測できる時間を求める
        if gsTrue.dr_counter == 1
            time.sc2wayget = scTrue.t_ur(gsTrue.ut_counter);
        end
        time.lastSearch = i;      
    end
    %% 5 宇宙機での状態量の推定 & Downlink
    % シグマ点列の計算
    scEstByScUkf.x_sp = Spacecraft.calcSigmaPoint(scEstByScUkf.X, scEstByScUkf.P,ukf);
    %% 観測がある場合, 1wayの観測がある場合，2wayの観測がある場合で場合分け      
    if time.list(i+1) < scTrue.t_ur(scTrue.ur_counter + 1) % 次の時刻まで観測がなかった場合
        % 状態量と誤差共分散行列を伝搬
        [scEstByScUkf.X, scEstByScUkf.P, scEstByScUkf.x_sp] = Spacecraft.timeUpdateUkf(scEstByScUkf.x_sp,constant, ukf, time.simDt, time.simDt, error);
    else % 観測があった時
        % 何度目の観測か
        scTrue.ur_counter = scTrue.ur_counter + 1;
        % 観測時刻までの時間
        time.scDt1 = scTrue.t_ur(scTrue.ur_counter) - time.list(i);
        % 観測からの時間
        time.scDt2 = time.list(i+1) - scTrue.t_ur(scTrue.ur_counter);
        if time.list(i+1) < time.sc2wayget
            % 観測は1way
            time.obsScType = 1;
        % 2wayの観測も得られる時，観測
        else
            % 観測は2way
            time.obsScType = 2;
%             time.obsScType = 1;
        end
        % 観測までは推定値と共分散を時間伝搬
        [scEstByScUkf.X, scEstByScUkf.P,~] = Spacecraft.timeUpdateUkf(scEstByScUkf.x_sp,constant, ukf, time.scDt1, time.simDt, error);
        % 観測量の計算
        scTrue.calcObservation_sc(scEstByScUkf,gsTrue,constant,error,sc,gs,time.obsScType); % time.obsScType=1: 1wayの観測, =2:2wayの観測
        % EKFで観測値を用いて推定値のupdate 
        scEstByScUkf.observationUpdateByScUkf(scTrue,earth, gsTrue,constant,time.obsScType,ukf,time)
        % シグマポイントを再計算 
        scEstByScUkf.x_sp = Spacecraft.calcSigmaPoint(scEstByScUkf.X, scEstByScUkf.P,ukf);
        % 観測から次ステップまでの推定値と共分散を時間伝搬. 
        [scEstByScUkf.X, scEstByScUkf.P,~] = Spacecraft.timeUpdateUkf(scEstByScUkf.x_sp,constant, ukf, time.scDt2, time.simDt, error);
        % downlink方向の決定. ダウンリンクを受ける時刻も計算
        [scTrue,gsTrue,earth] = scTrue.calcDownDirection(time.list(i+1),scTrue.state(:,i+1),scEstByScUkf.X(2:7),gsTrue,earth,time,constant);  
    end
    % 推定値を記録していく
    scEstByScUkf.clockError(i+1)= error.clock0- scEstByScUkf.X(1); %残りの誤差
    scEstByScUkf.state(:,i+1)= scEstByScUkf.X(2:7);
    scEstByScUkf.P_list(:,:,i+1) = scEstByScUkf.P;
    
    %% 6 地上局での状態量の推定 → 4でのuplinkにつながる
    % シグマポイントの計算
    scEstByGsUkf.x_sp = Spacecraft.calcSigmaPoint(scEstByGsUkf.X, scEstByGsUkf.P,ukf);
    % gsTrue.t_dr(gsTrue.dr_counter)が定義できないので，定義できるようになる前は場合を分けている
    if time.list(i+1) < scTrue.t_ur(1)
        % 状態量と誤差共分散行列を伝搬
        [scEstByGsUkf.X, scEstByGsUkf.P, scEstByGsUkf.x_sp] = Spacecraft.timeUpdateUkf(scEstByGsUkf.x_sp,constant, ukf, time.simDt, time.simDt, error);
    % 次の時刻まで観測がなかった時
    elseif time.list(i+1) < gsTrue.t_dr(gsTrue.dr_counter + 1)
        % 状態量と誤差共分散行列を伝搬
        [scEstByGsUkf.X, scEstByGsUkf.P, scEstByGsUkf.x_sp] = Spacecraft.timeUpdateUkf(scEstByGsUkf.x_sp,constant, ukf, time.simDt, time.simDt, error);
    % 次の時刻までに観測があった時．
    else
        gsTrue.dr_counter = gsTrue.dr_counter + 1;
        % 送信時刻の状態量を推定する．まず，送信時刻のアプリオリな状態量を得る
        if gsTrue.dr_counter == 1 % 初めての観測の時. 実際に送信した時刻は推定値からクロックのオフセットの分だけずれている
            % 探査機からのdownlinkがどの時刻のものと推定されているか
            time.scDtEstByGs = scTrue.t_dt(1) + scEstByGsUkf.clockError(i);
            % time.scDtEstByGsでの地上局が推定している宇宙機の状態量と誤差共分散行列を求める
            time.gsDt1 = time.scDtEstByGs - time.list(i);
            [scEstByGsUkf.X_dt, scEstByGsUkf.P_dt, scEstByGsUkf.x_sp_dt] = Spacecraft.timeUpdateUkf(scEstByGsUkf.x_sp,constant, ukf, time.gsDt1, time.simDt, error);      
        else %前の観測があった時刻から推定値を更新する %更新の間隔には，クロックのオフセットは乗らない         
            time.scDtEstByGsNew = scTrue.t_dt(gsTrue.dr_counter) + scEstByGsUkf.clockError(i);
            time.gsDt1 = time.scDtEstByGsNew - time.scDtEstByGs;
            time.scDtEstByGs = time.scDtEstByGsNew;
            [scEstByGsUkf.X_dt, scEstByGsUkf.P_dt, scEstByGsUkf.x_sp_dt] = Spacecraft.timeUpdateUkf(scEstByGsUkf.x_sp_dt,constant, ukf, time.gsDt1, time.simDt, error);                 
        end
        % 観測量の計算
        gsTrue.calcObservation_gs(scTrue,earth,constant,gs,sc,error);
        % EKFで観測量から推定量を計算する
        scEstByGsUkf.observationUpdateByGsUkf(gsTrue,earth,constant,ukf)
        % シグマポイントの再計算 % To Do
        scEstByGsUkf.x_sp_dt = Spacecraft.calcSigmaPoint(scEstByGsUkf.X_dt, scEstByGsUkf.P_dt,ukf);
        % 時刻time.list(i+1)での推定値と誤差共分散を求める(結果に使う)
        time.gsDt2 = time.list(i+1) - time.scDtEstByGs;
        [scEstByGsUkf.X, scEstByGsUkf.P, scEstByGsUkf.x_sp] = Spacecraft.timeUpdateUkf(scEstByGsUkf.x_sp_dt,constant, ukf, time.gsDt2, time.simDt, error);             
    end
    % 推定値を記録していく
    scEstByGsUkf.clockError(i+1)= error.clock0- scEstByGsUkf.X(1); %残りの誤差
    scEstByGsUkf.state(:,i+1)= scEstByGsUkf.X(2:7);
    scEstByGsUkf.P_list(:,:,i+1) = scEstByGsUkf.P;
    disp(i)  
end

showResult(scTrue,scEstByScUkf,scEstByGsUkf,error);
