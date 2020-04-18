% 地上局→探査機のone-wayの通信で，時計誤差を見積もる

function calcClockOffset(obj,time,earthState,gsState,scState,constant,i,clockErrorCorrection)
relDt = 1;
relTolDt = 0.0001;
iterNum = 1;
while (relDt > relTolDt) && (iterNum < 100)
% for i = 1:timeNumber
    t3 = time.list(i) ;
    t2 = t3 - obj.ltd(i) - clockErrorCorrection;
    %%  時刻t2での推定値の位置, 速度を求める 
    if t2 > time.list(1)
            % 構造体の状態量の内一番近い時刻のものを探す
            closeTimeIndex = i - round(t2/time.simDt);                    % 一番近い時刻がt(closeTimeIndex)
            closeTimeOffset = t2 - time.list(closeTimeIndex);      % 一番近い時刻から伝搬しなければいけない時間
             % earthについて伝播遅延(temp)時間前の状態量を得る．
            xve = [earthState.pos(:,closeTimeIndex);earthState.vel(:,closeTimeIndex)];
            k1e = orbitalState.twobody(xve,constant.sunMu);
            k2e = orbitalState.twobody(xve+0.5*closeTimeOffset*k1e,constant.sunMu);
            k3e = orbitalState.twobody(xve+0.5*closeTimeOffset*k2e,constant.sunMu);
            k4e = orbitalState.twobody(xve+closeTimeOffset*k3e,constant.sunMu);
            xve = xve + closeTimeOffset/6*(k1e+2*k2e+2*k3e+k4e); 
         % state構造体の状態量を時間帯範囲内にtimeDelayTemp時間遡った時刻が収まっていない時
     else
            % state構造体の初期時刻との差分.初期時刻から伝搬しなければならない時間(マイナス→逆方向伝搬になる)
            goBackTime = t2 -  time.list(1); 
            % 初期時刻から遡る
            xve =  [earthState.pos(:,1);earthState.vel(:,1)];
            % time.simDtずつ逆方向伝搬する．
            while goBackTime < 0
                if abs(goBackTime) > time.simDt
                    backTimeStep = -time.simDt;
                else
                    backTimeStep = goBackTime;
                end
                k1e = orbitalState.twobody(xve,constant.sunMu);
                k2e = orbitalState.twobody(xve+0.5*backTimeStep*k1e,constant.sunMu);
                k3e = orbitalState.twobody(xve+0.5*backTimeStep*k2e,constant.sunMu);
                k4e = orbitalState.twobody(xve+backTimeStep*k3e,constant.sunMu);
                xve = xve + backTimeStep/6*(k1e+2*k2e+2*k3e+k4e);
                goBackTime = goBackTime + time.simDt;
            end           
     end
     % ground sationの伝搬時間前の位置速度を求める
     [xg,rg] = groundState.earthRotation(gsState.pos(:,i), -obj.ltd(i), constant);
     
     % 時刻t3の宇宙機の時刻t2の地上局に対する相対位置を求める
     r2 = xve(1:3) + xg;
     v2 = xve(4:6) + rg;
     r23 = scState.pos(:,i) - r2;
%      r2r2Dot = r2.'*v2;
%      r3r3Dot = scState.pos(:,i).' * scState.vel(:,i);
     pDot23 = r23.' * v2/norm(r23);    
    % 時計誤差が求まる
     obj.clockError(i) = - (t3 - t2 - norm(r23)/constant.lightSpeed)/(1-pDot23/constant.lightSpeed) + clockErrorCorrection;
     relDt = abs(obj.clockError(i) - clockErrorCorrection);
     clockErrorCorrection = obj.clockError(i);
     obj.clockErrorLog(i,iterNum) = clockErrorCorrection; 
     iterNum = iterNum +1;
end


end