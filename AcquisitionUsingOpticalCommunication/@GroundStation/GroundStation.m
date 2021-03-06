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
          directionAccuracy_ut %uplinkの方向の精度(Rに使う)
          pointingError_ut %uplinkの送信指向誤差(回線計算に使う)
          ut2w_counter     %宇宙機が2wayの観測をできる何番目のuplinkか
          ut2w_counterList %宇宙機が2wayの観測をできる何番目のuplinkか(何度目のダウンリンク直後の探索の送信か)を格納する
                           %ut_counter番目の要素が何番目の宇宙機が2wayの観測をできる何番目のuplinkかに対応させる
          % 探索に関するパラメーター
          opnEstTempT_ut
          opnEstTempState_ut
          opnTrueTempT_ut
          opnTrueTempState_ut
          % 受信に関するパラメーター
          dr_counter           %downlinkを受信した回数
          dr_observability     %downlinkを観測できたかどうか. 1=観測できない. 2=SN比が低い 3=観測できた
          t_dr                 %downlinkを受信した時刻
          state_dr             %downlinkを受信した時刻の状態量
          lengthTrue_dr        %観測誤差を含まない測距情報
          lengthObserved_dr    %downlinkで観測される測距情報(宇宙機のクロック誤差が乗っている)
          length2wObserved_dr  %2wayの測距(距離換算)
          durationAtSc         %宇宙機での受信から送信までの時間(2wayの測距に使う)
          directionTrue_dr     %downlink時の観測誤差を含まない見かけの宇宙機の方向
          receivedPower_dr     %downlinkで受信した電力強度
          snr_dr
          directionAccuracy_dr %測角の精度
          directionObserved_dr %downlinkで観測される測角情報
          scAccel_dr           %downlinkされる宇宙機の加速度情報
          scRecAngle_dr        %downlinkされる宇宙機の測角
          scRecAngleAccuracy_dr%上の精度
          transUpAngle_dr      %downlinkされる地上局のuplinkでの送信方向
          transUpAngleAccuracy_dr %上の精度
          % 記録用
          t_drList
          snr_drList
          
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
            xform_ec2IAU = cspice_sxform('IAU_EARTH','ECLIPJ2000', time.t0Ephemeris +time.list );
            % get ehemeris data of the ground station in ECLIPJ2000, earth center
            obj.state = zeros(6,time.stepNum+1);
            for i_1 = 1: time.stepNum+1
                obj.state(:,i_1) = xform_ec2IAU(:,:,i_1)*bodyFixed;
            end
        end
        xvAtT = calcStateAtT_gs(obj,t,time,constant) % 時刻tでの状態量を得る
        [obj,earth] = search(obj,i,earth,gs,time,constant,error) 
        obj = calcObservation_gs(obj,scTrue,earth,constant,gs,sc,error,scEstbyGs,time,scEstbySc) % 観測量の計算
        
    end
    methods(Static)
        [opn_t,opn_state,dtlt] = calcTarget(t,gsAtT,eAtT,scAtT,spacecraft,time,constant,valueType)
        xv = earthRotation(pos0, t, constant) % t秒後の状態量の計算
        [gs,time,R] = setSearchArea(time,gs,SSD,scEstByGsSeqP,error) %1回の観測にかかる時間を求める
    end 
end