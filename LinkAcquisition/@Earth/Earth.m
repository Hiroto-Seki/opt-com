classdef Earth < handle
    properties
          t                       %時刻
          mu                      %中心天体
          state
          stateTrans              % 光を照射した時刻の状態量
    end
    methods 
        function obj = Earth(time,mu) %コンストラクタ 
            obj.t     = time.list;
            obj.mu = mu;
        end  
        getEphem(obj,time)  % 真値の場合はエフェメリスから状態量を取得する
        
    end
    methods(Static)
        xv = twobody(xv, mu,a_pertubation) %運動方程式
        xvAtT = calcStateE(sc,t,time)    % ある時刻の状態量を計算する
    end
end