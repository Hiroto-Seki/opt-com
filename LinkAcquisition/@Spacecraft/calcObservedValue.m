% 観測時間と観測量の計算をする
% 入力
% i: 何番目の計算か time.list(i)に対応する
% gsTrue.tTrans: 地上局が光を送信する時刻
% gsTrue.stateTrans: 地上局が光を送信する時刻の地上局の状態量
% eTrue.stateTrans: 地上局が光を送信する時刻の地球の状態量
% scTue.state: 探査機の軌道



function [scTrue,gsTrue] = calcObservedValue(scTrue,gsTrue,eTrue,i,constant,time,gs,sc,error)
  % 光を照射した時刻，位置と探査機の軌道から光が届く時の探査機の時刻と位置を計算する
  receive = GroundStation.calcTarget(gsTrue.tTrans(i),gsTrue.stateTrans(:,i),eTrue.stateTrans(:,i),scTrue,time,constant);
  ltd = receive.t - gsTrue.tTrans(i);
  %観測時の状態量
  scTrue.tReceive(i) = receive.t;
  scTrue.stateReceive(:,i) = receive.state;  
  % 地上局→探査機の送信光の指向誤差の計算
  % 実際の方向の計算
  relPos = scTrue.stateReceive(1:3,i) - eTrue.stateTrans(1:3,i) - gsTrue.stateTrans(1:3,i); % 相対位置 地上局→探査機
  azm = atan2(relPos(2),relPos(1));
  elv = atan(relPos(3) /(relPos(1)^2 + relPos(2)^2)^0.5 );
%   scatter(azm,elv)
  gsTrue.pointingErrorTrans(i) = ((gsTrue.azmTrans(i) - azm)^2 + (gsTrue.elvTrans(i) - elv)^2)^0.5;
  
  % QDセンサの観測誤差の計算
  % 指向誤差損失の計算
  Lp = calcPointingLoss(gsTrue.pointingErrorTrans(i),gs.gamma,gs.alpha,gs.tAperture,gs.wavelength);
  % 自由空間損失
  Ls = (gs.wavelength/(4 * pi * (ltd * constant.lightSpeed * 1e3)))^2;
  % 受信電力強度
  scTrue.receivedPower(i) = gs.peakPower * gs.tAntGain * gs.tEff * Lp * Ls * gs.atmosphereEff *  sc.rAntGain * sc.rEff; 
  % QDセンサの精度(1σ)
  QdAccuracy = (sc.qdIj^2 ...
      + 2 * constant.elementaryCharge * sc.qdId * sc.qdBw ...
      + 2 * constant.elementaryCharge * scTrue.receivedPower(i) *sc.qdS * sc.qdBw )^0.5 ...
      / (scTrue.receivedPower(i) * sc.qdS) * sc.qdFov;
  % 観測量の角度精度
  angleError = (QdAccuracy ^2 + error.stt^2)^0.5;
%   disp(angleError)
  % 観測量の計算
  % 誤差なし
  scTrue.lengthTrue(i) = ltd * constant.lightSpeed;
  xve = eTrue.stateTrans(:,i);
  xvg = gsTrue.stateTrans(:,i);
  xvs = scTrue.stateReceive(:,i);
  scTrue.azmTrue(i)  =  atan2(xve(2) + xvg(2) - xvs(2) + xvs(5)*ltd...
                                , xve(1) + xvg(1) - xvs(1) + xvs(4)*ltd);
  scTrue.elvTrue(i)  = atan( (xve(3) + xvg(3) - xvs(3) + xvs(6)*ltd)/...
                                   ((xve(2) + xvg(2) - xvs(2) + xvs(5) *ltd)^2 ...
                                  + (xve(1) + xvg(1) - xvs(1) + xvs(4) *ltd)^2)^0.5 );
  % 観測量誤差あり
  scTrue.lengthObserved(i) = scTrue.lengthTrue(i) + (error.initialClock + error.randomClock*randn) * constant.lightSpeed;
  scTrue.azmObserved(i)  =  scTrue.azmTrue(i)+angleError * randn;
  scTrue.elvObserved(i)  = scTrue.elvTrue(i)+ angleError * randn;
  
  % 加速度センサの値
  velAccel = Spacecraft.twobody(scTrue.stateReceive(:,i), constant.sunMu, 0);
  scTrue.accelObserved(:,i) = velAccel(4:6) + randn(3,1) * norm(velAccel(4:6))* error.stt + randn(3,1)* error.accel;
   
 
end