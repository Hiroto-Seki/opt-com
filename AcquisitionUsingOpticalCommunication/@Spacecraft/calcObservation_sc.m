% 宇宙機による観測値を計算する
% 入力
% 出力
% 中間パラメーター

function obj = calcObservation_sc(obj,scEst,gsTrue,constant,error,sc,gs,type)      
      %% 1way, 2wayに共通
      %% 1wayの測距
      % 観測値
      obj.lengthObserved_ur(obj.ur_counter)   =  obj.lengthTrue_ur(obj.ur_counter) + constant.lightSpeed * (error.clock0 + error.randomClock * randn);
      
      %% 慣性空間上での見かけのuplink受信方向
      % 慣性空間上での見かけのuplink受信方向の真値
      directionTrueI =  ( obj.eState_ur(1:3,obj.ur_counter) +  obj.gsState_ur(1:3,obj.ur_counter) - obj.state_ur(1:3,obj.ur_counter))...
                   + 1/constant.lightSpeed * obj.state_ur(4:6,obj.ur_counter) * norm(obj.eState_ur(1:3,obj.ur_counter) +  obj.gsState_ur(1:3,obj.ur_counter) - obj.state_ur(1:3,obj.ur_counter)) ; 
      % 正規化する
      directionTrueI_normalized = directionTrueI/norm(directionTrueI);
      %宇宙機の姿勢の真値(推定はしない) だいたい地上局方向に光学系の座標系が向くように設定する.(厳密には一致させなくても良い)
      % 姿勢は，宇宙機の軌道の推定値をもとに決定している. 軌道はここでは推定しない
      directionEstI = ( obj.eState_ur(1:3,obj.ur_counter) +  obj.gsState_ur(1:3,obj.ur_counter) -scEst.X(2:4))...
                   + 1/constant.lightSpeed * scEst.X(5:7) * norm(obj.eState_ur(1:3,obj.ur_counter) +  obj.gsState_ur(1:3,obj.ur_counter) - scEst.X(2:4)) ;      
      directionEstI_normalized = directionEstI/norm(directionEstI);
      %テキトーな値にするため，ランダムな成分を加えている
      randomAtt =  50 * 10^-6;
      rollTrue  =  asin(-directionEstI_normalized(2)) + randn * randomAtt;
      pitchTrue = atan2(directionEstI_normalized(1),directionEstI_normalized(3)) + randn *  randomAtt;
      yawTrue   = 0 + randn * randomAtt;
      % 光学系(Fight Laser-communication Terminal)座標系での真の見かけのuplink受信方向(見かけの相対位置)
      directionTrueFLT =  (Spacecraft.rotation(rollTrue,pitchTrue,yawTrue,1))\directionTrueI; 
%       directionTrueFLT_normalized =  (Spacecraft.rotation(rollTrue,pitchTrue,yawTrue,1))\directionTrueI_normalized; 
      % pointingLossの計算
      Lp = calcPointingLoss(gsTrue.pointingError_ut(obj.ur_counter),gs.gamma,gs.alpha,gs.tAperture,gs.wavelength_up);
      % 自由空間損失
      Ls = (gs.wavelength_up/(4 * pi * (obj.lengthTrue_ur(obj.ur_counter)  * 1e3)))^2;
      % 受信電力強度
      obj.receivedPower_ur(obj.ur_counter) = gs.slotPower * gs.tAntGain * gs.tEff * Lp * Ls * gs.atmosphereEff *  sc.rAntGain * sc.rEff; 
      % QDセンサーの精度 
      qdIl = obj.receivedPower_ur(obj.ur_counter) * sc.qdS; %入射光による電流
      Snr = (sc.qdGain * qdIl)^2 /...
          (sc.qdGain^2 * 2 * constant.elementaryCharge * (qdIl + sc.qdId) * sc.qdBw * sc.qdF + sc.qdIj^2);
      % QDセンサーの値(x,y)を姿勢の真値から求める. 誤差を含む
      qdXY_noError = - sc.fL /(directionTrueFLT(3)*1e3) * (directionTrueFLT(1:2)*1e3); %単位がmであることに注意
      % もし，QDセンサー内に入らなかったはいらなかったり,SN比が低すぎたら
      if norm(qdXY_noError) > sc.qdPhi || 1 > Snr
          disp("Uplink doesn't enter the FOV of QD")
          obj.ur_observability(obj.ur_counter) = 1;
      elseif sc.reqSnr_up > Snr
          disp("Uplink signal to noise ratio is too low")
          obj.ur_observability(obj.ur_counter) = 2;
      else
          disp("Uplink signal is observed")
          obj.ur_observability(obj.ur_counter) = 3;
      end 
      
      %% 宇宙機の姿勢を推定する
      % 　※宇宙機の軌道の推定値が悪いと，qdセンサーの値を姿勢に変換する精度も落ちてしまうので，宇宙機の軌道の精度が悪い場合は，sttのみの観測を使い，宇宙機の軌道の精度がいい場合は，qdも姿勢決定にしようする
      sttObserved = [rollTrue;pitchTrue;yawTrue] + randn(3,1)*error.stt;
      %   軌道の推定値のズレがどの程度観測に影響しそうか
      rollEst  = sttObserved(1);
      pitchEst = sttObserved(2);
      yawEst   = sttObserved(3);
      obj.directionAccuracy_ur(obj.ur_counter) = sqrt(error.stt^2 * 3/2  + ((1/Snr)*sc.qdFov)^2);

      
      
      % 見かけのuplink方向を記録する
%       obj.directionTrue_ur(:,obj.ur_counter) = [ atan2( directionTrueI_normalized(2),directionTrueI_normalized(1) ); asin(directionTrueI_normalized(3))];
      obj.directionTrue_ur(:,obj.ur_counter) = directionTrueI_normalized;
      
      % 計算の簡略化
%       obj.directionObserved_ur(:,obj.ur_counter) = obj.directionTrue_ur(:,obj.ur_counter) + randn(2,1)* obj.directionAccuracy_ur(obj.ur_counter);
      obj.directionObserved_ur(:,obj.ur_counter) = obj.directionTrue_ur(:,obj.ur_counter) + randn(3,1)* obj.directionAccuracy_ur(obj.ur_counter);
      
%       obj.directionObserved_ur(:,obj.ur_counter) = [ atan2( directionEstI_normalizedNew(2),directionEstI_normalizedNew(1) ); asin(directionEstI_normalizedNew(3))];  
      % 姿勢を記録する
      obj.attStateTrue_ur(:,obj.ur_counter) = [rollTrue;pitchTrue;yawTrue];       
      obj.attStateObserved_ur(:,obj.ur_counter) = [rollEst;pitchEst;yawEst]; 

         
    %% 加速度
    velAccel = CelestialBody.twobody(obj.state_ur(:,obj.ur_counter),constant.sunMu,0);
    obj.accelTrue_ur(:,obj.ur_counter) = velAccel(4:6);
    obj.accelObseved_ur(:,obj.ur_counter) = (Spacecraft.rotation(rollEst,pitchEst,yawEst,1)) * ( (Spacecraft.rotation(rollTrue,pitchTrue,yawTrue,1))\obj.accelTrue_ur(:,obj.ur_counter)) + randn(3,1) *  error.accel;
    
    
    %% 2way のみ
    if type == 2
%     % 2wayの測距(往復時間,滞在時間)
        obj.length2wObserved_ur(obj.ur_counter) = (obj.t_ur(obj.ur_counter) - obj.t_dt(obj.ur2w_counter) + error.randomClock * randn )*constant.lightSpeed ;
        obj.durationAtGs(obj.ur_counter) = gsTrue.t_ut(obj.ur_counter) - gsTrue.t_dr(obj.ur2w_counter);
      % 地上局の観測量も送信されている
      obj.recDownAngle_ur(:,obj.ur_counter) = gsTrue.directionObserved_dr(:,obj.ur2w_counter);
      obj.recDownAngleAccuracy_ur(obj.ur_counter) = gsTrue.directionAccuracy_dr(obj.ur2w_counter);
    end

end
