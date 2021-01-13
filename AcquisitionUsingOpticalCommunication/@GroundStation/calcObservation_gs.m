    % 地上局の観測値の計算
    %           dr_counter           %downlinkを受信した回数
    %           t_dr                 %downlinkを受信した時刻
    %           state_dr             %downlinkを受信した時刻の状態量
    %           lengthTrue_dr        %観測誤差を含まない測距情報
    %           lenghtObserved_dr    %downlinkで観測される測距情報(宇宙機のクロック誤差が乗っている)
    %           length2wTrue_dr      %2wayの測距(距離換算)
    %           length2wObserved_dr  %2wayの測距(距離換算)
    %           directionTrue_dr     %downlink時の観測誤差を含まない見かけの宇宙機の方向
    %           powerObserved_dr     %downlinkで受信した電力強度
    %           directionObserved_dr %downlinkで観測される測角情報
    %           scAccel_dr           %downlinkされる宇宙機の加速度情報



    function calcObservation_gs(obj,scTrue,earth,constant,gs,sc,error,scEstByGs,time,scEstBySc)
    % 何度目の受信か
    dr_counter = obj.dr_counter;

    % 地上局からの送信時の状態量
    t_ut    = obj.t_ut(dr_counter);
    xvgs_tu = obj.state_ut(:,dr_counter);
    xve_ut  = earth.state_ut(:,dr_counter);

    % 宇宙機の受信時の状態量
    t_ur    = scTrue.t_ur(dr_counter);
    xvsc_ur = scTrue.state_ur(:,dr_counter);
    
    % 宇宙機の受信時の推定状態量
    xvscEst_ur = scEstBySc.calcStateAtT_sc(t_ur,time,constant);

    % 宇宙機の送信時の状態量
    t_dt    = scTrue.t_dt(dr_counter);
    xvsc_dt = scTrue.state_dt(:,dr_counter);

    % 地上局の受信時の状態量
    t_dr    = obj.t_dr(dr_counter);
    xvgs_dr = obj.state_dr(:,dr_counter);
    xve_dr  = earth.state_dr(:,dr_counter);
    
    
    
    % 地上局が推定する．宇宙機の送信時の状態量
    xvsc_dt_estByGs = scEstByGs.calcStateAtT_sc(t_dt,time,constant);
    % 地上局が推定する宇宙機の方向
    gs2scD_estByGs  = xvsc_dt_estByGs - xve_dr - xvgs_dr; 

    
    

    %% 測角
    gs2scD = xvsc_dt - xve_dr - xvgs_dr;  % vector(position and velocity) of GroundStation To SpaceCraft when Downlink
    directionTrue = gs2scD(1:3) + (xve_dr(4:6) + xvgs_dr(4:6)) * ( norm(gs2scD(1:3))/constant.lightSpeed ) ;
    obj.directionTrue_dr(:,dr_counter) = [atan2(directionTrue(2),directionTrue(1));
                                        atan(directionTrue(3)/(directionTrue(1)^2 +directionTrue(2)^2 )^0.5)];
    

                                    
    % 観測誤差の計算
    pointingError_dt = scTrue.pointingError_dt(dr_counter);
    % 指向誤差損失
    Lp = calcPointingLoss(pointingError_dt,sc.gamma,sc.alpha,sc.aperture,sc.wavelength_down);
    % 自由空間損失
    Ls = (sc.wavelength_down/(4 * pi * (norm(gs2scD(1:3))  * 1e3)))^2;
    % 受信電力の計算
    obj.receivedPower_dr(dr_counter) = sc.slotPower * sc.tAntGain * sc.tEff * Lp * Ls * sc.atmosphereEff *  gs.rAntGain * gs.rEff; 
    qdIl = obj.receivedPower_dr(dr_counter) * gs.qdS; %入射光による電流
    Snr = (gs.qdGain * qdIl)^2 /...
          (gs.qdGain^2 * 2 * constant.elementaryCharge * (qdIl + gs.qdId) * gs.qdBw * gs.qdF + gs.qdIj^2);
    
    % 地上局の推定する方向の誤差
    gsRecDirectionError = acos(gs2scD_estByGs.' * gs2scD/ norm(gs2scD_estByGs)/norm(gs2scD));
      
    if gsRecDirectionError > gs.qdFov || 1 > Snr
        disp("Downlink is out of the FOV of ground station")
        obj.dr_observability(dr_counter) = 1;
    elseif gs.reqSnr_down > Snr
        disp("Downlink of signal to noise ratio is too low")
        obj.dr_observability(dr_counter) = 2;
    else
        disp("Downlink of signal is observed")
        obj.dr_observability(dr_counter) = 3;
    end
    
    % 観測誤差(地上局はQDの精度=慣性空間での測角精度とする)
    obj.directionAccuracy_dr(dr_counter) = gs.qdFov /Snr;
    obj.directionObserved_dr(:,dr_counter) = obj.directionTrue_dr(:,dr_counter) + randn(2,1) * obj.directionAccuracy_dr(dr_counter);
    
    % downlinkに載っている情報
    % uplinkを受けてから受信するまでの時間の分の移動を補正した値にする
    dtAtSc = t_dt - t_ur;
    % direction_utの補正
    posGs2Sc_Est = xvscEst_ur(1:3) - (xvgs_dr(1:3)+xve_dr(1:3));
    velGs2Sc_Est = xvscEst_ur(4:6) - (xvgs_dr(4:6)+xve_dr(4:6));
    dAzm_ut = dtAtSc * ( - velGs2Sc_Est(1) * sin(scTrue.transUpAngle_dt(1,dr_counter)) ...
        +  velGs2Sc_Est(2) * cos(scTrue.transUpAngle_dt(1,dr_counter)))...
        /((t_ur - t_ut)* constant.lightSpeed * cos(scTrue.transUpAngle_dt(2,dr_counter)));
    dElv_ut = dtAtSc * ( -   (velGs2Sc_Est(1) * posGs2Sc_Est(1) + velGs2Sc_Est(2) * posGs2Sc_Est(2))/(posGs2Sc_Est(1)^2+posGs2Sc_Est(2)^2)^0.5...
        * sin(scTrue.transUpAngle_dt(2,dr_counter)) ...
        +  velGs2Sc_Est(3) * cos(scTrue.transUpAngle_dt(2,dr_counter)))...
        /((t_ur - t_ut)* constant.lightSpeed);     
    obj.scAccel_dr(:,dr_counter) = scTrue.accel_dt(:,dr_counter);
    obj.scRecAngle_dr(:,dr_counter) = scTrue.recUpAngle_dt(:,dr_counter) - [dAzm_ut;dElv_ut]; %uplinkの受信角度
    obj.transUpAngle_dr(:,dr_counter) = scTrue.transUpAngle_dt(:,dr_counter) + [dAzm_ut;dElv_ut];
    obj.transUpAngleAccuracy_dr(:,dr_counter) = scTrue.transUpAngleAccuracy_dt(:,dr_counter);
    obj.scRecAngleAccuracy_dr(dr_counter) = scTrue.recUpAngleAccuracy_dt(dr_counter);    
    

    % 測距(1way) % クロックのオフセットは乗っている
    obj.lengthTrue_dr(dr_counter) = (t_dr - t_dt) * constant.lightSpeed;
    obj.lengthObserved_dr(dr_counter) = obj.lengthTrue_dr(dr_counter) + ( - error.clock0  + randn * error.randomClock) * constant.lightSpeed;

    % 測距(2way). 宇宙機の受信してから送信するのにかかる時間を含めている
    obj.length2wObserved_dr(dr_counter) = (t_dr - t_ut) * constant.lightSpeed + randn * error.randomClock * constant.lightSpeed;
    
    obj.durationAtSc(dr_counter) = t_dt - t_ur;

end