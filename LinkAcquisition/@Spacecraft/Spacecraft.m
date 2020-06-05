classdef Spacecraft < handle
    properties
          t 
          mu                      %Gravitational Constant of central body
          state
          tOpn
          stateAtTEstOpn         %gsTrue.tEstOpnでの状態量を求める
          tReceive               %受信時刻
          stateReceive           %受信時刻の状態量
          receivedPower          %受信電力(tReceivedに対応)
          azmObserved            %観測された方位角
          elvObserved            %観測された仰角
          lengthObserved         %観測された距離
          clockError             %(推定値)時計誤差
          
          
    end
    methods 
        function obj = Spacecraft(time,mu,refTime) %コンストラクタ 
            if strcmp(refTime,"spacecraft")
                obj.t = time.list;
            else 
                obj.tOpn = time.list;
            end
            obj.mu = mu;
            obj.stateAtTEstOpn = zeros(6,length(obj.t));
            obj.tReceive = zeros(1,length(obj.t));
            obj.stateReceive = zeros(6,length(obj.t));
        end 
        obj = calcOrbitTwoBody(obj,timestep,error)
    end
    methods(Static)
        xv = twobody(xv, mu,a_pertubation)
        xvAtT = calcStateSc(sc,t,time)    % ある時刻の状態量を計算する
        [scTrue,gsTrue] = calcObservedValue(scTrue,gsTrue,eTrue,i,constant,time,gs,sc,error)  %誤差なしの観測量を計算する
        [X_bar, P_bar] = updateState(X_hat,P, dt,mu) % STMを用いて状態量と誤差共分散行列を更新する
        A  = delFdelX(xv,mu) % EKFに用いる運動方程式の微分
        Y_bar = calcG(X_bar,xve,xvg,constant) %リファレンスの状態量の時の観測量を計算
        H_childa = delGdelX(X_bar,xve,xvg,constant) % EKFのための状態量の時の観測量を計算
    end
end