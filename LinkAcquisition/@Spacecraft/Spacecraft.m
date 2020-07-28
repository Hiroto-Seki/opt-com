classdef Spacecraft < handle
    properties
          mu                      %Gravitational Constant of central body
          % tに対応する状態量
          t 
          state                   % tに対応する状態量
          azmDown                % ダウンリンクの方位角
          elvDown                % ダウンリンクの仰角
          eDown                  % ダウンリンク先の地球の状態量
          gsDown                 % ダウンリンク先の地上局の状態量
          tDown                  % ダウンリンク先の時刻
          
          % tReceiveに対応する状態量
          tReceive               %受信時刻
          stateReceive           %受信時刻の状態量
          receivedPower          %受信電力(tReceivedに対応)
          azmObserved            %観測された方位角
          elvObserved            %観測された仰角
          lengthObserved         %観測された距離
          clockError             %(推定値)時計誤差
          azmTrue                %観測された方位角(誤差なし)
          elvTrue                %観測された仰角(誤差なし)
          lengthTrue             %観測された距離(誤差なし)

          
          
    end
    methods 
        function obj = Spacecraft(time,mu,refTime) %コンストラクタ 
            if strcmp(refTime,"spacecraft")
                obj.t = time.list;
            else 
                %obj.tOpn = time.list;
            end
            obj.mu = mu;
            obj.tReceive = zeros(1,length(obj.t));
            obj.stateReceive = zeros(6,length(obj.t));
        end 
        obj = calcOrbitTwoBody(obj,timestep,error)
        obj = calcDownDirection(obj,t,gs,e,scAtT,time,constant,i)
        
    end
    methods(Static)
        xv = twobody(xv, mu,a_pertubation)
        xvAtT = calcStateSc(sc,t,time)    % ある時刻の状態量を計算する
        [scTrue,gsTrue] = calcObservedValue(scTrue,gsTrue,eTrue,i,constant,time,gs,sc,error)  %誤差ありの観測量を計算する
        [X_bar, P_bar] = updateState(X_hat,P, dt,mu) % 探査機が推定する時にSTMを用いて状態量と誤差共分散行列を更新する
        [X_bar, P_bar] = updateState2(X_hat,P, dt,mu) % 地上局が推定する時にSTMを用いて状態量と誤差共分散行列を更新する
        A  = delFdelX(xv,mu) % 探査機が推定するEKFに用いる運動方程式の微分
        A  = delFdelX2(xv,mu) % 地上局が推定するEKFに用いる運動方程式の微分
        Y_bar = calcG(X_bar,xve,xvg,constant) %探査機が推定する際，リファレンスの状態量の時の観測量を計算
        Y_bar = calcG2(X_bar1,xve0,xvg0,xve2,xvg2,duration, constant) %地上局が推定する際，リファレンスの状態量の時の観測量を計算
        H_childa = delGdelX(X_bar,xve,xvg,constant) % 探査機が推定するEKFのための状態量の時の観測量を計算
        H_childa = delGdelX2(X_bar1,xve0,xvg0,xve2,xvg2,dt, constant) % 地上局が推定するEKFのための状態量の時の観測量を計算
        opn = calcTarget(t,gs,e,scAtT,time,constant)  
    end
end