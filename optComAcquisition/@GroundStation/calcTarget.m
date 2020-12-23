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

function [opn_t,opn_state,dtlt] = calcTarget(t,gsAtT,eAtT,spacecraft,time,constant)
     % 伝搬遅延誤差許容値
     tol = 1e-13;
     corr = 1;
     % イタレーション回数
     iter = 0;
     iterMax = 10;
     % 伝搬遅延の計算 
     dtlt_length = 0; %初期化.
     while (corr > tol) && (iter < iterMax)
         % 伝搬時間後の状態量を求める 
%          xvsc = spacecraft.calcStateAtT_sc(t + dtlt_length/constant.lightSpeed,time);
         % 試しにinterpで計算してみる
         xvsc = interp1(spacecraft.t, spacecraft.state.', t + dtlt_length/constant.lightSpeed, 'spline','extrap').';
         % RLT項を含まない伝搬遅延
         dtlt_length_New = norm(eAtT(1:3) + gsAtT(1:3) - xvsc(1:3) );
         corr = abs(dtlt_length_New - dtlt_length);
         dtlt_length = dtlt_length_New; 
         iter = iter + 1;
     end
     opn_t    = t + dtlt_length/constant.lightSpeed;
     opn_state = interp1(spacecraft.t, spacecraft.state.', opn_t, 'spline','extrap').';
%      opn_state = spacecraft.calcStateAtT_sc(opn_t,time);
     dtlt = dtlt_length/constant.lightSpeed;
end