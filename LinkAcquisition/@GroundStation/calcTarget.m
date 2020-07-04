% 時刻tでの地上局位置，地球位置と，探査機の軌道から，探査機に光が届く時刻とその時の探査機の状態量を求める
% 入力:
% t: 時刻
% gsAtT : 時刻tでの地上局
% eAtT: 時刻tでの地球
% sc: 探査機の軌道と時刻情報を含む
% time:
% constant
% i                : 地上局の時刻obj.t(i)に対応する値を計算する
% scEstGs          : 地上局が推定している宇宙機
% eTrue            : 地球
% 出力
% opn.t : 
% opn.state      : 

function opn = calcTarget(t,gsAtT,eAtT,sc,time,constant)
     % 伝搬遅延誤差許容値
     relTimeDelayError = 1e-4;
     timeDelayErrorTemp = 100;
     % 伝搬遅延の計算 
     timeDelayTemp = 1000;
     while timeDelayErrorTemp > relTimeDelayError
         % 伝搬時間後の状態量を求める       
         % 伝搬時刻後までのリストを持っている時
         if (t + timeDelayTemp) < sc.t(length(sc.t))
            % 構造体の状態量の内一番近い時刻のものを探す
            closeTimeIndex = round((t + timeDelayTemp - time.list(1))/time.simDt)+1;                    % 一番近い時刻がt(closeTimeIndex)
            closeTimeOffset = (t + timeDelayTemp) - time.list(closeTimeIndex);      % 一番近い時刻から伝搬しなければいけない時間
            % spacecraftについて伝搬時刻後(temp)時間後の状態量を得る
            xvsc = sc.state(:,closeTimeIndex);
            k1sc = sc.twobody(xvsc,sc.mu,0);
            k2sc = sc.twobody(xvsc+0.5*closeTimeOffset*k1sc,sc.mu,0);
            k3sc = sc.twobody(xvsc+0.5*closeTimeOffset*k2sc,sc.mu,0);
            k4sc = sc.twobody(xvsc+closeTimeOffset*k3sc,sc.mu,0);
            xvsc = xvsc + closeTimeOffset/6*(k1sc+2*k2sc+2*k3sc+k4sc); 
         % 伝搬時刻後の状態量を持っていない時
         else
            restTime = t + timeDelayTemp -  sc.t(length(sc.t)); 
            % 最終時刻から伝搬していく
            xvsc =  sc.state(:,length(sc.t));
            % time.simDtずつ伝搬する．
            while restTime > 0
                if restTime > time.simDt
                    timeStep = time.simDt;
                else
                    timeStep = restTime;
                end
                k1sc = sc.twobody(xvsc,sc.mu,0);
                k2sc = sc.twobody(xvsc+0.5*timeStep*k1sc,sc.mu,0);
                k3sc = sc.twobody(xvsc+0.5*timeStep*k2sc,sc.mu,0);
                k4sc = sc.twobody(xvsc+timeStep*k3sc,sc.mu,0);
                xvsc = xvsc + timeStep/6*(k1sc+2*k2sc+2*k3sc+k4sc);
                restTime = restTime - time.simDt;
            end           
         end
         % RLT項を含まない伝搬遅延
         timeDelayNew = 1/constant.lightSpeed * ...
         ((eAtT(1) + gsAtT(1) - xvsc(1))^2 + ...
          (eAtT(2) + gsAtT(2) - xvsc(2))^2 +....
          (eAtT(3) + gsAtT(3) - xvsc(3))^2)^0.5;
         timeDelayErrorTemp = abs(timeDelayNew - timeDelayTemp);
         timeDelayTemp = timeDelayNew;
     end
     opn.t    = t + timeDelayTemp;
     opn.state = xvsc;
end