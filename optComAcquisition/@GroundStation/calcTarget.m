% 時刻tでの地上局位置，地球位置と，探査機の軌道から，探査機に光が届く時刻とその時の探査機の状態量を求める
% 入力:
% t: 時刻
% gsAtT : 時刻tでの地上局
% eAtT: 時刻tでの地球
% scAtT: 時刻tでの宇宙機
% sc: 探査機の軌道と時刻情報を含む
% time:
% constant
% 出力
% opn_t : 
% opn_state      : 

function [opn_t,opn_state] = calcTarget(t,gsAtT,eAtT,spacecraft,time,constant)
     % 伝搬遅延誤差許容値
     relTimeDelayError = 1e-15;
     timeDelayErrorTemp = 10;
     % 伝搬遅延の計算 
     timeDelayTemp = 1000; %初期化.
     while timeDelayErrorTemp > relTimeDelayError
         % 伝搬時間後の状態量を求める 
         xvsc = spacecraft.calcStateAtT_sc(t + timeDelayTemp,time);
         % RLT項を含まない伝搬遅延
         timeDelayNew = 1/constant.lightSpeed * ...
         ((eAtT(1) + gsAtT(1) - xvsc(1))^2 + ...
          (eAtT(2) + gsAtT(2) - xvsc(2))^2 +....
          (eAtT(3) + gsAtT(3) - xvsc(3))^2)^0.5;
         timeDelayErrorTemp = abs(timeDelayNew - timeDelayTemp);
         timeDelayTemp = timeDelayNew;
     end
     opn_t    = t + timeDelayTemp;
     opn_state = xvsc;

end