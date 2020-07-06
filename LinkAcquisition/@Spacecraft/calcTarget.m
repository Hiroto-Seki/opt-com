% 時刻tでの地上局位置，地球位置と，探査機の軌道から，探査機に光が届く時刻とその時の探査機の状態量を求める
% 入力:
% t: 時刻
% gs : 地上局の軌道(class)
% e: 地球の軌道(class)
% sc: 時刻tでの探査機
% time:
% constant
% i                : 探査機の時刻obj.t(i)に対応する値を計算する
% 出力
% opn.t : 地上局に光が届く時間
% opn.stateGs      : 地上局に光が届いた時の地上局の位置
% opn.stateE       : 地上局に光が届いた時の地球の位置

function opn = calcTarget(t,gs,e,scAtT,time,constant)
     % 伝搬遅延誤差許容値
     relTimeDelayError = 1e-4;
     timeDelayErrorTemp = 100;
     % 伝搬遅延の計算 
     timeDelayTemp = 1000;
     while timeDelayErrorTemp > relTimeDelayError
         % 伝搬時間後の状態量を求める       
         % 伝搬時刻後までのリストを持っている時
         if (t + timeDelayTemp) < e.t(length(e.t))
            % 構造体の状態量の内一番近い時刻のものを探す
            closeTimeIndex = round((t + timeDelayTemp - time.list(1))/time.simDt)+1;     % 一番近い時刻がt(closeTimeIndex)
            closeTimeOffset = (t + timeDelayTemp) - time.list(closeTimeIndex);      % 一番近い時刻から伝搬しなければいけない時間
            % 地球について伝搬時間(temp)後の状態量を得る
            xve = e.state(:,closeTimeIndex);
            k1e = e.twobody(xve,e.mu,0);
            k2e = e.twobody(xve+0.5*closeTimeOffset*k1e,e.mu,0);
            k3e = e.twobody(xve+0.5*closeTimeOffset*k2e,e.mu,0);
            k4e = e.twobody(xve+closeTimeOffset*k3e,e.mu,0);
            xve = xve + closeTimeOffset/6*(k1e+2*k2e+2*k3e+k4e); 
            % 地上局について伝搬時間(temp)後の状態量を得る
            xvgs = gs.state(:,closeTimeIndex);
            xvgs = gs.earthRotation(xvgs(1:3), closeTimeOffset, constant);
   
         % 伝搬時刻後の状態量を持っていない時
         else
            restTime = t + timeDelayTemp -  e.t(length(e.t)); 
            % 最終時刻から伝搬していく
            xve =  e.state(:,length(e.t));
            xvgs = gs.state(:,length(e.t));
            xvgs = gs.earthRotation(xvgs(1:3), restTime, constant);
            % time.simDtずつ伝搬する．
            while restTime > 0
                if restTime > time.simDt
                    timeStep = time.simDt;
                else
                    timeStep = restTime;
                end
                k1e = e.twobody(xve,e.mu,0);
                k2e = e.twobody(xve+0.5*timeStep*k1e,e.mu,0);
                k3e = e.twobody(xve+0.5*timeStep*k2e,e.mu,0);
                k4e = e.twobody(xve+timeStep*k3e,e.mu,0);
                xve = xve + timeStep/6*(k1e+2*k2e+2*k3e+k4e);
                restTime = restTime - time.simDt;
            end           
         end
         % RLT項を含まない伝搬遅延
         timeDelayNew = 1/constant.lightSpeed * ...
         ((xve(1) + xvgs(1) - scAtT(1))^2 + ...
          (xve(2) + xvgs(2) - scAtT(2))^2 +....
          (xve(3) + xvgs(3) - scAtT(3))^2)^0.5;
         timeDelayErrorTemp = abs(timeDelayNew - timeDelayTemp);
         timeDelayTemp = timeDelayNew;
     end
     opn.t    = t + timeDelayTemp;
     opn.stateGs = xvgs;
     opn.stateE = xve;
end