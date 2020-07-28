%% 時刻t(i)の状態量をもとにダウンリンクの観測量を計算する

function calcObservedValue(obj, scTrue,scEst,i,constant,time,gs,sc,error)
% 探査機のダウンリンク方向の誤差の計算 
pointingError = ((scEst.azmDown(i) - scTrue.azmDown(i)).^2 + (scEst.elvDown(i) - scTrue.elvDown(i)).^2).^0.5;
Lp = calcPointingLoss(pointingError,sc.gamma,sc.alpha,sc.aperture,sc.wavelength);
ltd = scTrue.tDown(i) - time.list(i);
Ls = (sc.wavelength/(4 * pi * (ltd * constant.lightSpeed * 1e3)))^2;
obj.receivedPower(i) = sc.peakPower * sc.tAntGain * sc.tEff * Lp * Ls * sc.atmosphereEff *  gs.rAntGain * gs.rEff; 

% QDセンサの精度の計算
QdAccuracy = (gs.qdIj^2 ...
      + 2 * constant.elementaryCharge * gs.qdId * gs.qdBw ...
      + 2 * constant.elementaryCharge * obj.receivedPower(i) *gs.qdS * gs.qdBw )^0.5 ...
      / (obj.receivedPower(i) * gs.qdS) * gs.qdFov;


% 距離情報(時計誤差を含む)
obj.lengthObserved(i) =  (ltd - error.initialClock + error.randomClock*randn) * constant.lightSpeed;
% 角度情報(誤差なし)
direction = (scTrue.state(1:3,i) - scTrue.eDown(1:3,i) - scTrue.gsDown(1:3,i) ) + ltd * (scTrue.eDown(4:6,i) + scTrue.gsDown(4:6,i) );

obj.azmObserved(i) = atan2(direction(2),direction(1)) + QdAccuracy * randn;
obj.elvObserved(i) = atan(direction(3)/ (direction(1)^2 +direction(2)^2)^0.5 ) + QdAccuracy * randn;

 
end