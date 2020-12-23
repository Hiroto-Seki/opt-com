classdef Spacecraft < handle
%     このクラスのオブジェクト
%     scTrue: 真値
%     scInitiallGuess: 初期推定値
%     scEstByScEkf: 宇宙機による推定値(EKF)の状態量
%     scEstByGsEkf: 地上局による推定値(EKF)の状態量
%     scEstByScBatch: 宇宙機による推定値(Batch)の状態量
%     scEstByGsBatch: 地上局による推定値(Batch)の状態量
    properties
          % ------------全てのオブジェクトに共通であるもの-----------
          mu                     %Gravitational Constant of central body
          t                      % 基準時刻
          state                  % tに対応する位置・速度
          clockError             % 残りの時計誤差 
          % 時刻ごとに，ランダムに姿勢を与える．→QDセンサー上での送信方向の誤差を計算する．姿勢は推定しない．(ダイナミクスが早いので，観測が追い付かなそうなので)
          % ------------scTrueにのみあるもの (観測に関するもの)--------
          ur_counter             %uplinkを受信した回数を数える
          t_ur                   %uplinkを受信する時刻
          state_ur               %観測時の宇宙機の位置・速度
          attStateTrue_ur        %t_urに対応する姿勢(ロール・ピッチ・ヨー) 
          attStateObserved_ur    %t_urに対応する姿勢 
          lengthTrue_ur          %距離の真値
          lengthObserved_ur      %観測される距離
          ur2w_counter         %2wayが観測された回数
          length2wObserved_ur%観測される距離(2way)
          durationAtGs           %地上局が受信してから送信するまでの時間
          directionTrue_ur       %誤差がない時の慣性空間上での見かけの方向(azimuth,elevation)
          directionObserved_ur   %観測結果から計算される慣性空間上での見かけの方向(sttとqdの誤差が含まれる)
          directionAccuracy_ur   %観測誤差共分散行列の計算にしようする．測角の精度
          receivedPower_ur       %受信される電力
          accelTrue_ur           %加速度の真値
          accelObseved_ur        %加速度の観測値
          gsT_ur                 %uplinkに載っている，送信時刻
          eState_ur              %uplinkに載っている，送信時刻の地球の状態量   
          gsState_ur             %uplinkに載っている, 送信時刻の地上局の状態量
          transDirection_ur      %uplinkに載っている，送信方向(誤差を持つ)
          % ----------scTrueにあるもの(送信に関するもの)--------
          dt_counter             % downlinkを送信した回数を数える 
          t_dt                   % downlinkを送信した時刻  
          state_dt               % downlinkを送信した時刻の状態量(観測量の計算に用いる)
          accel_dt               % downlinkに載っている，送信時の加速度
          pointingError_dt       % downlinkの送信方向誤差
          tRec_dt                % downlinkが受信されるはずの時刻
          % --------- scEstByScEkf, scEstByGsEkf, scEstByScBatch,
          % scEstByGsBatch にあるもの (推定に使う) ------
          X                      % 状態量の推定値をまとめたベクトル 
          Y                      % 観測値
          y                      % Y - Y*
          P                      % 誤差共分散行列
          R2wGs                      % 観測誤差共分散行列(gsの観測)
          R1wSc                    % 観測誤差共分散行列(scの観測)
          R2wSc                    % 観測誤差共分散行列(scの観測)
          H                      % batchの時使う
          P_list                 % 観測誤差共分散行列のリスト
          X_dt                   % scEstByGsEkfが使用するダウンリンクを送信した時刻の推定状態量
          P_dt                   % scEstByGsEkfが使用するダウンリンクを送信した時刻の誤差共分散行列
          % ----------UKFに使うもの---------------------
          x_sp                   % Xのsigma point
          y_sp                   % シグマポイントに対応するY
          x_sp_dt
          
          
    end
 
    methods 
        function obj = Spacecraft(time,mu) %コンストラクタ
            obj.t          = time.list;
            obj.mu         = mu;
            obj.state      = zeros(6,length(obj.t));
            obj.clockError = zeros(1,length(obj.t));
        end 
        obj = calcOrbitTwoBody(obj,error)
        xvAtT = calcStateAtT_sc(obj,t,time)  
        [obj,gsTrue] = receiveUplink(obj,gsTrue,earth,constant,time) %uplinkを受信する．その時の観測量を求める
        obj = calcObservation_sc(obj,scEst,gsTrue,constant,error,sc,gs,type) %観測量の計算 obj = scTrue, type=1:1way, type=2:2way
        observationUpdateBySc(obj,scTrue,earth,gsTrue,constant,type) % (宇宙機による)EKFで観測量を用いて推定値と誤差共分散を更新. 1wayと2wayでtype分け
        observationUpdateByScUkf(obj,scTrue,earth, gsTrue,constant,type,ukf) % (宇宙機による)UKFで観測量を用いて推定値と誤差共分散を更新. 1wayと2wayでtype分け
        [obj,gsTrue,eTrue] = calcDownDirection(obj,t,scTrueAtT,scEstAtT,gsTrue,eTrue,scAtT,time,constant) % obj = scTrue
        observationUpdateByGs(obj,gsTrue,earth,constant) % (地上局による)EKFでの観測を用いて推定値と誤差共分散を更新
        observationUpdateByGsUkf(obj,gsTrue,earth,constant,ukf) %% (地上局による)UKFでの観測を用いて推定値と誤差共分散を更新
        
        
    end
    methods(Static)
        rotationMatrix = rotation(roll,pitch,yaw,order) %order=1:x軸→y軸→z軸の順番に回転する, order=2:z軸→y軸→x軸の順番に回転する
        A  = delFdelX(xv,mu) %運動方程式の微分
        [roll, pitch, yaw,P] = attDetermination(stt,sttError,qd,qdError,directionEstI,scfl)
        Y_star = calcG1w_ur(X_star,xve,xvg,constant,mu) % 推定値の時の観測量を計算(1way, 宇宙機による推定)
        Y_star = calcG2w_ur(X_star,xvet,xvgt,xver,xvgr,dtAtGs, dt2w, constant,mu)
        Y_star = calcG_dr(X_star,xv_ut,xv_dr,dtAtSc,constant,mu) % 推定値の時の観測量を計算(2way, 地上局による推定)
        [X, P] = timeUpdate(X, P, mu, Dt, dt)
        [X_new, P_new, x_sp_new] = calcTimeUpdateUkf(x_sp,mu, ukf,Dt,dt, error) % UKFに使う．シグマ点列と誤差共分散を時間更新
        H = delGdelX1w_ur(X_star,xve,xvg,constant,mu)   % 観測方程式の微分(1way, 宇宙機による推定)
        H = delGdelX2w_ur(X_star,xvet,xvgt,xver,xvgr, dt2w, constant,mu)
        H = delGdelX_dr(X_star,xv_ut,xv_dr,dtAtSc,constant,mu);  % 観測方程式の微分(2way, 地上局による推定)
        [opn_t,opn_stateE,opn_stateGs] = calcTarget(t,gs,e,scAtT,time,constant) % ダウンリンクが届く時刻とその時刻の地球,地上局の位置を求める
        X_sp = calcSigmaPoint(X,P,ukf) % シグマ点列の計算
%           xvAtT = calcStateAtT_sc(spacecraft,t,time)
%         [scTrue,gsTrue] = calcObservedValue(scTrue,gsTrue,eTrue,i,constant,time,gs,sc,error)  %one wayでの誤差ありの観測量を計算する
%         [scTrue,gsTrue] = calcObservedValue2way(scTrue,gsTrue,eTrue,i,constant,time,gs,sc,error,scTrans,gsReceiveNum); % 2wayの観測    
%         [X_bar, P_bar] = updateState(X_hat,P, dt,mu) % 探査機が推定する時にSTMを用いて状態量と誤差共分散行列を更新する
%         [X_bar, P_bar] = updateState2(X_hat,P, dt,mu) % 地上局が推定する時にSTMを用いて状態量と誤差共分散行列を更新する
%         A  = delFdelX(xv,mu) % 探査機が推定するEKFに用いる運動方程式の微分
%         A  = delFdelX2(xv,mu) % 地上局が推定するEKFに用いる運動方程式の微分
%         Y_bar = calcG(X_bar,xve,xvg,constant) %探査機が推定する際，リファレンスの状態量の時の観測量を計算
%         Y_bar = calcG2(X_bar1,xve0,xvg0,xve2,xvg2,duration, constant) %地上局が推定する際，リファレンスの状態量の時の観測量を計算
%         Y_bar2w = calcG2w(X_bar2w,xve,xvg,xveDr,xvgDr,dt,constant);
%         H_childa = delGdelX(X_bar,xve,xvg,constant) % 探査機の1wayでの観測の時の，ノミナル状態量の観測量微分
%         H_childa = delGdelX2(X_bar1,xve0,xvg0,xve2,xvg2,dt, constant) % 地上局の2wayでの観測の時の，ノミナル状態量の観測量微分
%         H_childa2w = delGdelX2w(xvs,xve,xvg,xveDr,xvgDr,dt,constant) %探査機の2wayでの観測の時の，ノミナル状態量の観測量微分
%         opn = calcTarget(t,gs,e,scAtT,time,constant)  
    end
end