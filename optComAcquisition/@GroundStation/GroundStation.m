classdef GroundStation < handle
% このクラスのオブジェクト
%   gs1True: 1つ目の地上局
    properties
          lat      %緯度
          lon      %経度 
          alt      %高度
          t        %基準時刻
          state    %基準時刻の状態量
          % 送信に関するパラメーター
          ut_counter       %uplinkを送信した回数
          t_ut             %uplinkを送信した時刻
          state_ut         %uplinkを送信した時刻の状態量
          direction_ut     %uplinkの方向(1:azimuth,2:elevation)
          pointingError_ut %uplinkの送信指向誤差
          % 探索に関するパラメーター
          opnEstTempT_ut
          opnEstTempState_ut
          opnTrueTempT_ut
          opnTrueTempState_ut
          % 受信に関するパラメーター
          dr_counter           %downlinkを受信した回数
          t_dr                 %downlinkを受信した時刻
          state_dr             %downlinkを受信した時刻の状態量
          lengthTrue_dr        %観測誤差を含まない測距情報
          lenghtObserved_dr    %downlinkで観測される測距情報(宇宙機のクロック誤差が乗っている)
          directionTrue_dr     %downlink時の観測誤差を含まない見かけの宇宙機の方向
          powerObserved_dr     %downlinkで受信した電力強度
          directionObserved_dr %downlinkで観測される測角情報
          scAccel_dr           %downlinkされる宇宙機の加速度情報
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
        end
        xvAtT = calcStateAtT_gs(obj,t,time,constant) % 時刻tでの状態量を得る
%         obj = calcObservedValue(obj,scTrue,scEst,i,constant,time,gs,sc,error) % 観測量の計算
    end
    methods(Static)
        [opn_t,opn_state] = calcTarget(t,gsAtT,eAtT,sc,time,constant)
        [gsTrue,earth] = search(i,gsTrue,earth,gs,time,constant) 
        xv = earthRotation(pos0, t, constant) % t秒後の状態量の計算
    end 
end