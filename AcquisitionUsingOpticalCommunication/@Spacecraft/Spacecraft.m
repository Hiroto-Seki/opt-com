classdef Spacecraft < handle
%     このクラスのオブジェクト
%     scTrue: 真値
%     scInitiallGuess: 初期推定値
%     scEstByScEkf: 宇宙機による推定値(EKF)の状態量
%     scEstByGsEkf: 地上局による推定値(EKF)の状態量
%     scEstByScUkf: 宇宙機による推定値(UKF)の状態量
%     scEstByGsUkf: 地上局による推定値(UKF)の状態量
    properties
          % ------------全てのオブジェクトに共通であるもの-----------
%           mu                     % Gravitational Constant of central body
          t                      % 基準時刻
          state                  % tに対応する位置・速度
          clockError             % 残りの時計誤差 
          % 時刻ごとに，ランダムに姿勢を与える．→QDセンサー上での送信方向の誤差を計算する．姿勢は推定しない．(ダイナミクスが早いので，観測が追い付かなそうなので)
          % ------------scTrueにのみあるもの (観測に関するもの)--------
          ur_counter             %uplinkを受信した回数を数える
          ur_observability         %uplinkを観測できたかどうか. 1=観測できない. 2=SN比が低い 3=観測できた
          t_ur                   %uplinkを受信する時刻
          state_ur               %観測時の宇宙機の位置・速度
          attStateTrue_ur        %t_urに対応する姿勢(ロール・ピッチ・ヨー) 
          attStateObserved_ur    %t_urに対応する姿勢 
          lengthTrue_ur          %距離の真値
          lengthObserved_ur      %観測される距離
          ur2w_counter           %2wayが観測された回数
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
          recDownAngle_ur        %uplinkに載っている，地上局の観測量(downlinkの測角)
          recDownAngleAccuracy_ur%上の精度
          % ----------scTrueにあるもの(送信に関するもの)--------
          dt_counter             % downlinkを送信した回数を数える 
          t_dt                   % downlinkを送信した時刻  
          state_dt               % downlinkを送信した時刻の状態量(観測量の計算に用いる)
          accel_dt               % downlinkに載っている，宇宙機の観測量(宇宙機の加速度)
          recUpAngle_dt          % downlinkに載っている，宇宙機の観測量(uplinkの受信角度)
          recUpAngleAccuracy_dt  % 上の精度
          transUpAngle_dt        % downlinkに載っている，宇宙機の観測量(uplinkの送信角度)
          transUpAngleAccuracy_dt% 上の精度
          pointingError_dt       % downlinkの送信方向誤差
          targetError_dt         % downlinkの目標方向の誤差
          tRec_dt                % downlinkが受信されるはずの時刻
          % --------- scEstByScEkf, scEstByGsEkf(推定に使う) ------
          X                      % 状態量の推定値をまとめたベクトル 
          Y                      % 観測値
          y                      % Y - Y*
          P                      % 誤差共分散行列
          R                  % 観測誤差共分散行列(gsの観測)
          H                      % batchの時使う
          P_list                 % 観測誤差共分散行列のリスト
          X_dt                   % scEstByGsEkfが使用するダウンリンクを送信した時刻の推定状態量
          P_dt                   % scEstByGsEkfが使用するダウンリンクを送信した時刻の誤差共分散行列
          useObs                 % それぞれの観測を使うかどうか記述している
          % ----------UKFに使うもの---------------------
          x_sp                   % Xのsigma point
          y_sp                   % シグマポイントに対応するY
          x_sp_dt
%           % ----------Information filterに使うもの
%           S
%           W
%           S_dt
%           W_dt
          % 記録用
          targetError_dtList
          pointError_dtList
          
    end
 
    methods 
        function obj = Spacecraft(time) %コンストラクタ
            obj.t          = time.list;
%             obj.mu         = mu;
            obj.state      = zeros(6,length(obj.t));
            obj.clockError = zeros(1,length(obj.t));
        end 
        obj = calcOrbitTwoBody(obj,sunMu,error)
        xvAtT = calcStateAtT_sc(obj,t,time,constant)  
        [obj,gsTrue] = receiveUplink(obj,gsTrue,earth,constant,time) %uplinkを受信する．その時の観測量を求める
        obj = calcObservation_sc(obj,scEst,gsTrue,constant,error,sc,gs,type) %観測量の計算 obj = scTrue, type=1:1way, type=2:2way
        observationUpdateByScEkf(obj,scTrue,earth,gsTrue,constant,type,time,ekf) % (宇宙機による)EKFで観測量を用いて推定値と誤差共分散を更新. 1wayと2wayでtype分け
        observationUpdateByScUkf(obj,scTrue,earth, gsTrue,constant,type,ukf,time) % (宇宙機による)UKFで観測量を用いて推定値と誤差共分散を更新. 1wayと2wayでtype分け
        [obj,gsTrue,eTrue] = calcDownDirection(obj,t,scTrueAtT,scEstAtT,gsTrue,eTrue,scAtT,time,constant,error) % obj = scTrue
        observationUpdateByGsEkf(obj,gsTrue,earth,constant,ekf,scTrue,time) % (地上局による)EKFでの観測を用いて推定値と誤差共分散を更新
        observationUpdateByGsUkf(obj,gsTrue,earth,constant,ukf,scTrue) %% (地上局による)UKFでの観測を用いて推定値と誤差共分散を更新
        observationUpdateByScESIRF(obj,scTrue,earth, gsTrue,constant,type,time)
        
    end
    methods(Static)
        rotationMatrix = rotation(roll,pitch,yaw,order) %order=1:x軸→y軸→z軸の順番に回転する, order=2:z軸→y軸→x軸の順番に回転する
        A  = delFdelX(xv,mu) %運動方程式の微分
        [roll, pitch, yaw,P] = attDetermination(stt,sttError,qd,qdError,directionEstI,scfl)
        Y_star = calcG_ur(X_star,xvet,xvgt,xver,xvgr,dtAtGs, dt2w, constant,time,type) % 推定値の時の観測量を計算(1way, 宇宙機による推定)
        Y_star = calcG_dr(X_star,xv_ut,xv_dr,dtAtSc,constant) % 推定値の時の観測量を計算(2way, 地上局による推定)
        xvsc = timeUpdate_sc(xvsc,mu, Dt, dt)
        [X, P] = timeUpdateEkf(X, P, constant, Dt, dt,error)
        [X_new, P_new, x_sp_new] = timeUpdateUkf(x_sp,constant,ukf,Dt,dt, error) % UKFに使う．シグマ点列と誤差共分散を時間更新
        H = delGdelX_ur(X_star,xvet,xvgt,xver,xvgr, dt2w, constant,type,time)   % 観測方程式の微分(1way, 宇宙機による推定)
        H = delGdelX_dr(X_star,xv_ut,xv_dr,dtAtSc,constant);  % 観測方程式の微分(2way, 地上局による推定)
        [Yv,YStarv,Hm,Rm,obsNum,sigmaN] = alignReqInfo4Est(Y,YStar,H,R,obsType,estType,reqList); %y, Y, H, Rを過不足なく並べる. type="1u:1wayのuplink", "2u:2wayのuplink","2d:2wayのdownlink"
        [opn_t,opn_stateE,opn_stateGs] = calcTarget(t,gs,e,scAtT,time,constant) % ダウンリンクが届く時刻とその時刻の地球,地上局の位置を求める
        X_sp = calcSigmaPoint(X,P,ukf) % シグマ点列の計算
        dxdt = twobody_stateAndSTM(t,x,mu)
        STM = getSTM(x0,mu,t0, tf)
        [X, S] = timeUpdateESIRF(X, S, constant, Dt, dt,error)
        dGamma = calcDGamma(delFdelX,dt,tau)
    end
end