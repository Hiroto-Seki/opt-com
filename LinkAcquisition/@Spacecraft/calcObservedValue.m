% 観測時間と観測量の計算をする
% 入力
% i: 何番目の計算か time.list(i)に対応する
% gsTrue.tTrans: 地上局が光を送信する時刻
% gsTrue.stateTrans: 地上局が光を送信する時刻の地上局の状態量
% eTrue.stateTrans: 地上局が光を送信する時刻の地球の状態量
% scTue.state: 探査機の軌道



function [scTrue,gsTrue] = calcObservedValue(scTrue,gsTrue,eTrue,i,constant,time,gs,sc,error)
% 伝搬時間許容誤差
 relTimeDelayError = 1e-6;
 timeDelayErrorTemp = 100;
% 伝搬時間の計算
 timeDelayTemp = 1/constant.lightSpeed * ...
         ((eTrue.stateTrans(1,i) + gsTrue.stateTrans(1,i) - scTrue.state(1,i))^2 + ...
          (eTrue.stateTrans(2,i) + gsTrue.stateTrans(2,i) - scTrue.state(2,i))^2 +....
          (eTrue.stateTrans(3,i) + gsTrue.stateTrans(3,i) - scTrue.state(3,i))^2)^0.5;
      
  while timeDelayErrorTemp > relTimeDelayError
     % 伝搬時間後の状態量を求める       
     xvsc = Spacecraft.calcStateSc(scTrue,gsTrue.tTrans(i) + timeDelayTemp,time);
     % RLT項を含まない伝搬遅延
     timeDelayNew = 1/constant.lightSpeed * ...
          ((eTrue.stateTrans(1,i) + gsTrue.stateTrans(1,i) - xvsc(1))^2 + ...
          (eTrue.stateTrans(2,i) + gsTrue.stateTrans(2,i) - xvsc(2))^2 +....
          (eTrue.stateTrans(3,i) + gsTrue.stateTrans(3,i) - xvsc(3))^2)^0.5;
     timeDelayErrorTemp = abs(timeDelayNew - timeDelayTemp);
     timeDelayTemp = timeDelayNew;
  end 
  %観測時の状態量
  scTrue.tReceive(i) = gsTrue.tTrans(i)  + timeDelayTemp;
  scTrue.stateReceive(:,i) = xvsc;  
  % 地上局→探査機の送信光の指向誤差の計算
  % 実際の方向の計算
  relPos = scTrue.stateReceive(1:3,i) - eTrue.stateTrans(1:3,i) - gsTrue.stateTrans(1:3,i); % 相対位置 地上局→探査機
  azm = atan2(relPos(2),relPos(1));
  elv = atan(relPos(3) /(relPos(1)^2 + relPos(2)^2)^0.5 );
  gsTrue.pointingErrorTrans(i) = ((gsTrue.azmTrans(i) - azm)^2 + (gsTrue.elvTrans(i) - elv)^2)^0.5;
  
  % QDセンサの観測誤差の計算
  % 指向誤差損失の計算
  Lp = calcPointingLoss(gsTrue.pointingErrorTrans(i),gs);
  % 自由空間損失
  Ls = (gs.wavelength/(4 * pi * (timeDelayTemp * constant.lightSpeed * 1e3)))^2;
  % 受信電力強度
  scTrue.receivedPower(i) = gs.peakPower * gs.tAntGain * gs.tEff * Lp * Ls * gs.atmosphereEff *  sc.rAntGain * sc.rEff; 
  % QDセンサの精度(1σ)
  QdAccuracy = (sc.qdIj^2 ...
      + 2 * constant.elementaryCharge * sc.qdId * sc.qdBw ...
      + 2 * constant.elementaryCharge * scTrue.receivedPower(i) * sc.qdBw )^0.5 ...
      / (scTrue.receivedPower(i) * sc.qdS) * sc.qdFov;
  % 観測量の計算
  scTrue.lengthObserved(i) = (timeDelayTemp + error.initialClock + error.randomClock*randn) * constant.lightSpeed;
  xve = eTrue.stateTrans(:,i);
  xvg = gsTrue.stateTrans(:,i);
  xvs = scTrue.stateReceive(:,i);
  scTrue.azmObserved(i)  =  atan2(xve(2) + xvg(2) - xvs(2) + xvs(5)*timeDelayTemp...
                                , xve(1) + xvg(1) - xvs(1) + xvs(4)*timeDelayTemp) ...
                                +QdAccuracy * randn;
  scTrue.elvObserved(i)  = atan( (xve(3) + xvg(3) - xvs(3) + xvs(6)*timeDelayTemp)/...
                                   ((xve(2) + xvg(2) - xvs(2) + xvs(5) *timeDelayTemp)^2 ...
                                  + (xve(1) + xvg(1) - xvs(1) + xvs(4) *timeDelayTemp)^2)^0.5 )...
                                  + QdAccuracy * randn;
  
    

 
end