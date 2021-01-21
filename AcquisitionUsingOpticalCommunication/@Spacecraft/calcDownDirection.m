% 宇宙機から地上局のダウンリンク方向を計算する
% 3種類の値を計算する
% 宇宙機から地上局の実際の方向: 地上局での観測量の計算に用いる．
% 宇宙機から実際に送信した方向: 実際に送信した方向を計算する．pointing lossの計算に用いる
% 宇宙機が推定した地上局方向: 宇宙機から実際に送信した方向．これに姿勢の誤差が乗って，実際の方向が決まる

% 入力
% obj = scTrue
% t   = 送信時刻
% scTrueAtT: 送信時刻の宇宙機の位置の真値
% scEstAtT: 送信時刻の宇宙機の位置の推定値
% gsTrue: 地上局のクラス
% eTrue: 地球のクラス
% time
% constant

% 出力(obj = scTrue)
% obj.t_dt(obj.dt_counter)
% obj.state_dt(:,obj.dt_counter)
% obj.tRec_dt(obj.dt_counter)          : ダウンリンクが地上局に届く時刻の計算
% obj.pointingError_dt(obj.dt_counter) :pointing Errorの計算
% gsTrue.state_dr(:,obj.dt_counter)
% eTrue.state_dr(:,obj.dt_counter)

function [obj,gsTrue,earth] = calcDownDirection(obj,t,scTrueAtT,scEstAtT,gsTrue,earth,time,constant,error)
    % 何度目の送信か
    dt_counter = obj.dt_counter + 1;
    obj.dt_counter = dt_counter;
    
    % 宇宙機の受信回数と送信回数が一致しているのか確認(念のため)
    if obj.dt_counter == obj.ur_counter
    else
        disp("dt_counter and ur_counter is not same")
    end
    
    % 送信した時刻での宇宙機の状態量(真値)を求める．(地上局での観測量の計算に用いる)
    obj.t_dt(dt_counter) = t;
    obj.state_dt(:,dt_counter) = scTrueAtT;
    
    % 宇宙機の観測量をdownlinkして地上局と共有する(To Do: uplink受信時刻とDownlink送信時刻の差分を補正した方がいいが，interpを使うと逆にブレるので，時間差が小さいので，今回は補正しない)
    obj.accel_dt(:,dt_counter)   = obj.accelObseved_ur(:,dt_counter); %加速度
    obj.recUpAngle_dt(:,dt_counter) = obj.directionObserved_ur(:,dt_counter); %uplinkの受信角度
    obj.transUpAngle_dt(:,dt_counter) = obj.transDirection_ur(:,dt_counter);
    % 観測の分散も送信する
    obj.recUpAngleAccuracy_dt(dt_counter) = obj.directionAccuracy_ur(dt_counter);
    obj.transUpAngleAccuracy_dt(dt_counter) = gsTrue.directionAccuracy_ut(dt_counter);
    
    % 時刻tで宇宙機(真値)から送信された光が地上局に届く時刻とその時刻での地球・地上局の位置・速度を求める
    opnDownTrue = Spacecraft.calcTarget(t,gsTrue,earth,scTrueAtT,time,constant);
    obj.tRec_dt(dt_counter) = opnDownTrue.t;
    gsTrue.t_dr(dt_counter) = obj.tRec_dt(dt_counter); 
    gsTrue_dr = opnDownTrue.stateGs;
    gsTrue.state_dr(:,dt_counter) = gsTrue_dr;
    eTrue_dr  = opnDownTrue.stateE;
    earth.state_dr(:,dt_counter) = eTrue_dr;
    % 目標方向の真値を求める(光行差を含む)
    scTrue2gsI = eTrue_dr(1:3) + gsTrue_dr(1:3) - scTrueAtT(1:3); %光行差含まない
    directionTrueI = scTrue2gsI + scTrueAtT(4:6) * norm(scTrue2gsI) / constant.lightSpeed ; %光行差を含む
   
     % 時刻tで宇宙機(推定値)から送信された光が地上局に届く時刻とその時刻での地球・地上局の位置・速度を求める
    opnDownEst  = Spacecraft.calcTarget(t,gsTrue,earth,scEstAtT,time,constant);
    gsEst_dr    = opnDownEst.stateGs;
    eEst_dr     = opnDownEst.stateE; 
    
    %% 実際に送信した方向を求める
    % 推定する慣性系での方向を求める
    % 慣性座標系での宇宙機から地上局方向の推定値
    scEst2gsI = eEst_dr(1:3) + gsEst_dr(1:3) - scEstAtT(1:3); %光行差を含まない
    directionEstI = scEst2gsI + scEstAtT(4:6) * norm(scEst2gsI)/constant.lightSpeed;
    % 推定する姿勢で，FLT座標系に変換した方向 (姿勢の推定値はuplinkを受信した時刻のものと同じにしている)
    rollEst  = obj.attStateObserved_ur(1,dt_counter);
    pitchEst = obj.attStateObserved_ur(2,dt_counter);
    yawEst   = obj.attStateObserved_ur(3,dt_counter);
    directionEstFLT =  (Spacecraft.rotation(rollEst,pitchEst,yawEst,1))\directionEstI;      
    
    % 実際の姿勢で慣性系での方向を求める(姿勢の真値はuplinkを受信した時刻のものと同じにしている) 制御誤差も載せる
    rollTrue  = obj.attStateTrue_ur(1,dt_counter) + error.scPoint * randn ;
    pitchTrue = obj.attStateTrue_ur(2,dt_counter) + error.scPoint * randn;
    yawTrue   = obj.attStateTrue_ur(3,dt_counter) + error.scPoint * randn;    
    directionTransI = (Spacecraft.rotation(rollTrue,pitchTrue,yawTrue,1))* directionEstFLT ;

    % 送信方向の誤差を求める
    pointingError = abs( acos(directionTrueI.' * directionTransI/( norm(directionTrueI)* norm(directionTransI)) ));
    obj.pointingError_dt(dt_counter) = pointingError;
    
    
    
    
    
    end