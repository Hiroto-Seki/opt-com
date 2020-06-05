% calculate laser transmit direction using gs estimates value
% 入力:
% i                : 地上局の時刻obj.t(i)に対応する値を計算する
% scEstGs          : 地上局が推定している宇宙機
% eTrue            : 地球
% 出力
% obj.stateEstOpn  : obj.t(i)に照射された光が宇宙機に届いたと推定される時の宇宙機の位置速度
% obj.tEstOpn      : obj.t(i)に照射された光が宇宙機に届くと推定される時刻

function calcTarget(obj,i,scEstGs,eTrue,time,constant)
     % 伝搬遅延誤差許容値
     relTimeDelayError = 1e-4;
     timeDelayErrorTemp = 100;
     % 伝搬遅延の計算 
     timeDelayTemp = 1/constant.lightSpeed * ...
         ((eTrue.state(1,i) + obj.state(1,i) - scEstGs.state(1,i))^2 + ...
          (eTrue.state(2,i) + obj.state(2,i) - scEstGs.state(2,i))^2 +....
          (eTrue.state(3,i) + obj.state(3,i) - scEstGs.state(3,i))^2)^0.5;
     while timeDelayErrorTemp > relTimeDelayError
         % 伝搬時間後の状態量を求める       
         % 伝搬時刻後までのリストを持っている時
         if (obj.t(i) + timeDelayTemp) < scEstGs.t(length(scEstGs.t))
            % 構造体の状態量の内一番近い時刻のものを探す
            closeTimeIndex = i + round(timeDelayTemp/time.simDt);                    % 一番近い時刻がt(closeTimeIndex)
            closeTimeOffset = (time.list(i) + timeDelayTemp) - time.list(closeTimeIndex);      % 一番近い時刻から伝搬しなければいけない時間
            % spacecraftについて伝搬時刻後(temp)時間後の状態量を得る
            xvsc = scEstGs.state(:,closeTimeIndex);
            k1sc = scEstGs.twobody(xvsc,scEstGs.mu,0);
            k2sc = scEstGs.twobody(xvsc+0.5*closeTimeOffset*k1sc,scEstGs.mu,0);
            k3sc = scEstGs.twobody(xvsc+0.5*closeTimeOffset*k2sc,scEstGs.mu,0);
            k4sc = scEstGs.twobody(xvsc+closeTimeOffset*k3sc,scEstGs.mu,0);
            xvsc = xvsc + closeTimeOffset/6*(k1sc+2*k2sc+2*k3sc+k4sc); 
         % 伝搬時刻後の状態量を持っていない時
         else
            restTime = (obj.t(i) + timeDelayTemp) -  scEstGs.t(length(scEstGs.t)); 
            % 最終時刻から伝搬していく
            xvsc =  scEstGs.state(:,length(scEstGs.t));
            % time.simDtずつ伝搬する．
            while restTime > 0
                if restTime > time.simDt
                    timeStep = time.simDt;
                else
                    timeStep = restTime;
                end
                k1sc = scEstGs.twobody(xvsc,scEstGs.mu,0);
                k2sc = scEstGs.twobody(xvsc+0.5*timeStep*k1sc,scEstGs.mu,0);
                k3sc = scEstGs.twobody(xvsc+0.5*timeStep*k2sc,scEstGs.mu,0);
                k4sc = scEstGs.twobody(xvsc+timeStep*k3sc,scEstGs.mu,0);
                xvsc = xvsc + timeStep/6*(k1sc+2*k2sc+2*k3sc+k4sc);
                restTime = restTime - time.simDt;
            end           
         end
         % RLT項を含まない伝搬遅延
         timeDelayNew = 1/constant.lightSpeed * ...
         ((eTrue.state(1,i) + obj.state(1,i) - xvsc(1))^2 + ...
          (eTrue.state(2,i) + obj.state(2,i) - xvsc(2))^2 +....
          (eTrue.state(3,i) + obj.state(3,i) - xvsc(3))^2)^0.5;
         timeDelayErrorTemp = abs(timeDelayNew - timeDelayTemp);
         timeDelayTemp = timeDelayNew;
     end
     obj.tEstOpn(:,i)     = obj.t(i) + timeDelayTemp;
     obj.stateEstOpn(:,i) = xvsc;
end