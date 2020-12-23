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

function opn = calcTarget(t,gsTrue,earth,scAtT,time,constant)
%      % 伝搬遅延誤差許容値
%      relTimeDelayError = 1e-15;
%      timeDelayErrorTemp = 10;
%      % 伝搬遅延の計算 
%      timeDelayTemp = 1000;
%      while timeDelayErrorTemp > relTimeDelayError
%          % 伝搬時間後の状態量を求める
%         xve  = earth.calcStateAtT_cb(t + timeDelayTemp,time);
%         xvgs = gsTrue.calcStateAtT_gs(t + timeDelayTemp,time,constant);  
%         timeDelayNew = 1/constant.lightSpeed * ...
%          ((xve(1) + xvgs(1) - scAtT(1))^2 + ...
%           (xve(2) + xvgs(2) - scAtT(2))^2 +....
%           (xve(3) + xvgs(3) - scAtT(3))^2)^0.5;
%         timeDelayErrorTemp = abs(timeDelayNew - timeDelayTemp);
%         timeDelayTemp = timeDelayNew;
%      end
%      opn.t    = t + timeDelayTemp;
%      opn.stateGs = xvgs;
%      opn.stateE = xve;
     
     % 伝搬遅延誤差許容値
     tol = 1e-13;
     corr = 1;
     % イタレーション回数
     iter = 0;
     iterMax = 10;
     % 伝搬遅延の計算 
     dtlt_length = 0; %初期化.
     while (corr > tol) && (iter < iterMax)
         % 伝搬時間後の状態量を補間して計算する
         xve = interp1(earth.t, earth.state.', t + dtlt_length/constant.lightSpeed, 'spline','extrap').';
         xvgs = interp1(gsTrue.t, gsTrue.state.', t + dtlt_length/constant.lightSpeed, 'spline','extrap').';
         % 伝搬遅延
         dtlt_length_New = norm(xve(1:3) + xvgs(1:3) - scAtT(1:3) );
         corr = abs(dtlt_length_New - dtlt_length);
         dtlt_length = dtlt_length_New; 
         iter = iter + 1;
     end
     opn.t = t + dtlt_length/constant.lightSpeed;
     opn.stateGs = interp1(gsTrue.t, gsTrue.state.', opn.t, 'spline','extrap').';
     opn.stateE  = interp1(earth.t, earth.state.', opn.t, 'spline','extrap').';
end