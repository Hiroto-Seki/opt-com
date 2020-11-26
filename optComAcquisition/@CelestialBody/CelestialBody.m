classdef CelestialBody < handle
    % このクラスのオブジェクト
    % eTrue
    % saturn
    properties
        % 全てのオブジェクトに共通のプロパティ
          name               % 名前
          t                  % 基準時刻
          mu                 % 中心天体
          state              % 基準時刻に対応する状態量
        % earthのみにあるプロパティ(送受信に関するパラメーター)
          t_ut       
          state_ut
          t_dr
          state_dr
    end
    methods 
        function obj = CelestialBody(time,mu,Name) %コンストラクタ 
            obj.t     = time.list;
            obj.mu = mu;
            obj.name = Name;
        end  
        getEphem(obj,time)  % 真値の場合はエフェメリスから状態量を取得する
        xvAtT = calcStateAtT_cb(obj,t,time) 
    end
    methods(Static)
        xv = twobody(xv, mu, a_pertubation) % 運動方程式. 宇宙機も，同じダイナミクスに従う
%         xvAtT = calcStateE(sc,t,time)    % ある時刻の状態量を計算する
    end
end