% 相対論効果の項は入れていないので，あとで入れる予定
% 真値の計算

function fromSc2Gs(obj,time,earthState,gsState,scState,constant,clockError,i)


     % 伝搬遅延誤差許容値
     relTimeDelayError = 1e-4;
     timeDelayErrorTemp = 100;
     % 伝搬遅延の計算 
     timeDelayTemp = 1/constant.lightSpeed * ...
         ((earthState.pos(1,i) + gsState.pos(1,i) - scState.pos(1,i))^2 + ...
          (earthState.pos(2,i) + gsState.pos(2,i) - scState.pos(2,i))^2 +....
          (earthState.pos(3,i) + gsState.pos(3,i) - scState.pos(3,i))^2)^0.5;
     while timeDelayErrorTemp > relTimeDelayError
         % 伝搬時間だけ遡った時間の地球の位置速度を求める         
         % state構造体の状態量を時間帯範囲内にtimeDelayTemp時間遡った時刻が収まっている時
         if (time.list(i) - timeDelayTemp) > time.list(1)
             % 構造体の状態量の内一番近い時刻のものを探す
            closeTimeIndex = i - round((time.list(i) - timeDelayTemp)/time.simDt);                    % 一番近い時刻がt(closeTimeIndex)
            closeTimeOffset = (time.list(i) -timeDelayTemp) - time.list(closeTimeIndex);      % 一番近い時刻から伝搬しなければいけない時間
             % earthについて伝播遅延(temp)時間前の状態量を得る．
            xve = [earthState.pos(:,closeTimeIndex);earthState.vel(:,closeTimeIndex)];
            k1e = orbitalState.twobody(xve,constant.sunMu,0);
            k2e = orbitalState.twobody(xve+0.5*closeTimeOffset*k1e,constant.sunMu,0);
            k3e = orbitalState.twobody(xve+0.5*closeTimeOffset*k2e,constant.sunMu,0);
            k4e = orbitalState.twobody(xve+closeTimeOffset*k3e,constant.sunMu,0);
            xve = xve + closeTimeOffset/6*(k1e+2*k2e+2*k3e+k4e); 
            
%             % ground stationについて伝播遅延(temp)時間前の状態量を得る．
%             [xg,vg] = groundState.earthRotation(gsState.pos(:,closeTimeIndex), closeTimeOffset, constant);
%             xvg = [xg;vg];
         % state構造体の状態量を時間帯範囲内にtimeDelayTemp時間遡った時刻が収まっていない時
         else
            % state構造体の初期時刻との差分.初期時刻から伝搬しなければならない時間(マイナス→逆方向伝搬になる)
            goBackTime = (time.list(i) - timeDelayTemp) -  time.list(1); 
            % 初期時刻から遡る
            xve =  [earthState.pos(:,1);earthState.vel(:,1)];
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
                goBackTime = goBackTime + time.simDt;
            end           
         end
         % ground sationの伝搬時間前の位置速度を求める
         [xg,vg] = groundState.earthRotation(gsState.pos(:,i), -timeDelayTemp, constant);
         xvg = [xg;vg];
         % RLT項を含まない伝搬遅延
         timeDelayNew = 1/constant.lightSpeed * ...
         ((xve(1) + xvg(1) - scState.pos(1,i))^2 + ...
          (xve(2) + xvg(2) - scState.pos(2,i))^2 +....
          (xve(3) + xvg(3)- scState.pos(3,i))^2)^0.5;
         timeDelayErrorTemp = abs(timeDelayNew - timeDelayTemp);
         timeDelayTemp = timeDelayNew;
     end
     obj.ltd(i) = timeDelayTemp + clockError; %時計誤差の分だけ観測量はずれる
     obj.length(i) = obj.ltd(i) * constant.lightSpeed;
    % 角度も計算
    obj.azimuth(i) = atan2(xve(2) + xvg(2) - scState.pos(2,i) + scState.vel(2,i)/constant.lightSpeed *timeDelayTemp*constant.lightSpeed...
                                , xve(1) + xvg(1) - scState.pos(1,i) + scState.vel(1,i)/constant.lightSpeed *timeDelayTemp*constant.lightSpeed);
    obj.elevation(i) = atan( (xve(3) + xvg(3) - scState.pos(3,i) + scState.vel(3,i)/constant.lightSpeed *timeDelayTemp*constant.lightSpeed)/...
                                   ((xve(2) + xvg(2) - scState.pos(2,i) + scState.vel(2,i)/constant.lightSpeed *timeDelayTemp*constant.lightSpeed)^2 ...
                                  + (xve(1) + xvg(1) - scState.pos(1,i) + scState.vel(1,i)/constant.lightSpeed *timeDelayTemp*constant.lightSpeed)^2)^0.5 );
    
%  end
end