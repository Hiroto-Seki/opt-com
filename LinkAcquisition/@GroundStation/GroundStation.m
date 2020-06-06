classdef GroundStation < handle
    properties
          lat                     %緯度
          lon
          alt
          t                       %時刻
          state
          % レーザー照射探索に関わるパラメータ
          tEstOpn               %光が相手の宇宙機に届くと推定される時刻 time_estimate_oponent
          stateEstOpn           %光が相手の宇宙機に届くと推定される時刻の推定位置
          tTrans                %真値に一番近い方向を照射する時の時刻
          stateTrans            %trans時刻の地上局位置速度
          azmTrans              %光照射方位角
          elvTrans              %光照射仰角
          pointingErrorTrans    %光照射角度誤差→受信電力の計算に用いる
          
    end
    methods 
        function obj = GroundStation(gs,constant,time) 
            %コンストラクタ，　緯度，経度，高度，時刻から慣性空間上の位置速度の状態量を得る
            obj.lat = gs.lat;
            obj.lon = gs.lon;
            obj.alt = gs.alt;
            obj.t   = time.list;
            bodyFixed = cspice_georec(gs.lat,gs.lon,gs.alt,constant.earthRadius,constant.earthF);
            bodyFixed = [bodyFixed;0;0;0]; 
            xform_ec2IAU = cspice_sxform('IAU_EARTH','ECLIPJ2000',time.list);
            % get ehemeris data of the ground station in ECLIPJ2000, earth center
            obj.state = zeros(6,time.stepNum+1);
            for i_1 = 1: time.stepNum+1
                obj.state(:,i_1) = xform_ec2IAU(:,:,i_1)*bodyFixed;
            end
            obj.tEstOpn = zeros(1,length(obj.t));
            obj.stateEstOpn = zeros(6,length(obj.t));
            obj.tTrans = zeros(1,length(obj.t));
            obj.stateTrans = zeros(6,length(obj.t));
            obj.azmTrans = zeros(1,length(obj.t));
            obj.elvTrans = zeros(1,length(obj.t));
            obj.pointingErrorTrans = zeros(1,length(obj.t));
        end
        
    end
    methods(Static)
        opn = calcTarget(t,gsAtT,eAtT,sc,time,constant)
        [gsTrue,eTrue] = search(i,gsTrue,eTrue,gs,time,constant,opnEstTemp,opnTrueTemp) 
        xv = earthRotation(pos0, t, constant) % t秒後の状態量の計算
    end
end