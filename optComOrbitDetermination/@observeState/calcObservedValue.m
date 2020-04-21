% calc length and direction from sc to gs observed value
function calcObservedValue(obj,time,ephemData,scState,constant,error,i)
     % 伝搬遅延誤差許容値
     relTimeDelayError = 1e-4;
     timeDelayErrorTemp = 100;
     % 伝搬遅延の計算 
     timeDelayTemp = 1/constant.lightSpeed * ...
         ((ephemData.earth(1,i) + ephemData.gs(1,i) - scState.pos(1,i))^2 + ...
          (ephemData.earth(2,i) + ephemData.gs(2,i) - scState.pos(2,i))^2 +....
          (ephemData.earth(3,i) + ephemData.gs(3,i) - scState.pos(3,i))^2)^0.5;
     while timeDelayErrorTemp > relTimeDelayError
         % 伝搬時間だけ遡った時間の地球の位置速度を求める         
         % state構造体の状態量を時間帯範囲内にtimeDelayTemp時間遡った時刻が収まっている時
         if (time.list(i) - timeDelayTemp) > time.list(1)
             % 構造体の状態量の内一番近い時刻のものを探す
            closeTimeIndex = i - round(timeDelayTemp/time.simDt);                    % 一番近い時刻がt(closeTimeIndex)
            closeTimeOffset = (time.list(i) -timeDelayTemp) - time.list(closeTimeIndex);      % 一番近い時刻から伝搬しなければいけない時間
             % earthについて伝播遅延(temp)時間前の状態量を得る．
            xve = ephemData.earth(:,closeTimeIndex);
            k1e = orbitalState.twobody(xve,constant.sunMu,0);
            k2e = orbitalState.twobody(xve+0.5*closeTimeOffset*k1e,constant.sunMu,0);
            k3e = orbitalState.twobody(xve+0.5*closeTimeOffset*k2e,constant.sunMu,0);
            k4e = orbitalState.twobody(xve+closeTimeOffset*k3e,constant.sunMu,0);
            xve = xve + closeTimeOffset/6*(k1e+2*k2e+2*k3e+k4e); 
            % ground stationについて伝播遅延(temp)時間前の状態量を得る．
            xvg = ephemData.gs(:,closeTimeIndex);
            k1g = groundState.calcEarthRotation(xvg, constant);
            k2g = groundState.calcEarthRotation(xvg+0.5*closeTimeOffset*k1g, constant);
            k3g = groundState.calcEarthRotation(xvg+0.5*closeTimeOffset*k2g,constant);
            k4g = groundState.calcEarthRotation(xvg+closeTimeOffset*k3g,constant);
            xvg = xvg + closeTimeOffset/6*(k1g+2*k2g+2*k3g+k4g); 
         % state構造体の状態量を時間帯範囲内にtimeDelayTemp時間遡った時刻が収まっていない時
         else
            % state構造体の初期時刻との差分.初期時刻から伝搬しなければならない時間(マイナス→逆方向伝搬になる)
            goBackTime = (time.list(i) - timeDelayTemp) -  time.list(1); 
            % 初期時刻から遡る
            xve =  ephemData.earth(:,1);
            xvg =  ephemData.gs(:,1);
            % time.simDtずつ逆方向伝搬する．
            while goBackTime < 0
                if abs(goBackTime) > time.simDt
                    backTimeStep = -time.simDt;
                else
                    backTimeStep = goBackTime;
                end
                k1e = orbitalState.twobody(xve,constant.sunMu,0);
                k2e = orbitalState.twobody(xve+0.5*backTimeStep*k1e,constant.sunMu,0);
                k3e = orbitalState.twobody(xve+0.5*backTimeStep*k2e,constant.sunMu,0);
                k4e = orbitalState.twobody(xve+backTimeStep*k3e,constant.sunMu,0);
                xve = xve + backTimeStep/6*(k1e+2*k2e+2*k3e+k4e);
                k1g = groundState.calcEarthRotation(xvg, constant);
                k2g = groundState.calcEarthRotation(xvg+0.5*backTimeStep*k1g, constant);
                k3g = groundState.calcEarthRotation(xvg+0.5*backTimeStep*k2g,constant);
                k4g = groundState.calcEarthRotation(xvg+backTimeStep*k3g,constant);
                xvg = xvg + backTimeStep/6*(k1g+2*k2g+2*k3g+k4g); 
                goBackTime = goBackTime + time.simDt;
            end           
         end
         % RLT項を含まない伝搬遅延
         timeDelayNew = 1/constant.lightSpeed * ...
         ((xve(1) + xvg(1) - scState.pos(1,i))^2 + ...
          (xve(2) + xvg(2) - scState.pos(2,i))^2 +....
          (xve(3) + xvg(3)- scState.pos(3,i))^2)^0.5;
         timeDelayErrorTemp = abs(timeDelayNew - timeDelayTemp);
         timeDelayTemp = timeDelayNew;
     end
     obj.xve(:,i) = xve;
     obj.xvg(:,i) = xvg;
     obj.ltd(i) = timeDelayTemp + error.clock; %時計誤差の分だけ観測量はずれる
     obj.length(i) = obj.ltd(i) * constant.lightSpeed;
    % 角度も計算
    obj.azimuth(i) = atan2(xve(2) + xvg(2) - scState.pos(2,i) + scState.vel(2,i)*timeDelayTemp...
                                , xve(1) + xvg(1) - scState.pos(1,i) + scState.vel(1,i)*timeDelayTemp) + randn * error.angle;
    obj.elevation(i) = atan( (xve(3) + xvg(3) - scState.pos(3,i) + scState.vel(3,i)*timeDelayTemp)/...
                                   ((xve(2) + xvg(2) - scState.pos(2,i) + scState.vel(2,i) *timeDelayTemp)^2 ...
                                  + (xve(1) + xvg(1) - scState.pos(1,i) + scState.vel(1,i) *timeDelayTemp)^2)^0.5 ) + randn * error.angle;
    
%  end
end