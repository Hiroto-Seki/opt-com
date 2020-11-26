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
          % ------------scTrueにのみあるもの (観測に関するもの)--------
          attState               % tに対応する姿勢(ロール・ピッチ・ヨー)  
          % 時刻ごとに，ランダムに姿勢を与える．→QDセンサー上での送信方向の誤差を計算する．姿勢は推定しない．(ダイナミクスが早いので，観測が追い付かなそうなので)
          % ------------scTrueにのみあるもの (観測に関するもの)--------
          ur_counter             %uplinkを受信した回数を数える
          t_ur                   %uplinkを受信する時刻
          state_ur               %観測時の宇宙機の位置・速度
          lengthTrue_ur          %距離の真値
          lengthObserved_ur      %観測される距離
          directionTrue_ur       %誤差がない時の慣性空間上での見かけの方向(正規化するかは要検討)
          directionObserved_ur   %観測される見かけ上の方向(sttとqdの誤差が含まれる)
          receivedPower_ur       %受信される電力
          accelTrue_ur           %加速度の真値
          accelObseved_ur        %加速度の観測値
          gsT_ur                 %uplinkに載っている，送信時刻
          eState_ur              %uplinkに載っている，送信時刻の地球の状態量   
          gsState_ur             %uplinkに載っている, 送信時刻の地上局の状態量
          transDirection_ur      %uplinkに載っている，送信方向(誤差を持つ)
          % ----------scTrue, scEstByScEkfにあるもの (送信に関するもの)--------
          dt_counter             % downlinkを送信した回数を数える (scTrueに代表させる)
          t_dt                   % downlinkを送信した時刻        (scTrueに代表させても良い)
          tRec_dt                % downlinkが受信されるはずの時刻
          gsStateRec_dt          % downlinkが受信されるはずの時刻の地上局の状態量
          eStateRec_dt           % downlinkが受信されるはずの時刻の地球の状態量
          direction_dt           % qdセンサー上のダウンリンク方向 (真値は真の方向，推定値は，実際に送信した方向→これに真の姿勢の回転行列をかけて，送信誤差を求める)
          % --------- scEstByScEkf, scEstByGsEkf, scEstByScBatch,
          % scEstByGsBatch にあるもの (推定に使う) ------
          X                      % 状態量の推定値をまとめたベクトル
          Y                      % 観測値
          P                      % 誤差共分散行列
          R                      % 観測誤差共分散行列
          phi                    % STM   
          X_bar                  % 時間伝搬後のリファレンスの状態量
          H_childa               % 観測方程式のリファレンス状態量での微分
          H                      %batchの時使う
          P_bar                  % 時間伝搬後の誤差共分散行列
          Y_bar                  % リファレンスの状態量での観測値
          y                      % 観測値とリファレンスの差分
          K                      % カルマンゲイン
          P_list                 % 観測誤差共分散行列のリスト
    end
 
    methods 
        function obj = Spacecraft(time,mu) %コンストラクタ
            obj.t          = time.list;
            obj.mu         = mu;
            obj.state      = zeros(6,length(obj.t));
            obj.attState   = zeros(3,length(obj.t));
            obj.clockError = zeros(1,length(obj.t));
        end 
        obj = calcOrbitTwoBody(obj,error)
        xvAtT = calcStateAtT_sc(obj,t,time)  
    end
    methods(Static)
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