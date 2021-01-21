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

function [opn_t,opn_state,dtlt] = calcTarget(t,gsAtT,eAtT,scAtT,spacecraft,time,constant,valueType)
     % 伝搬遅延誤差許容値
     tol = 1e-13;
     corr = 1;
     % イタレーション回数
     iter = 0;
     iterMax = 10;
     % 伝搬遅延の計算 
     dtlt_length = 0; %初期化.
     while (corr > tol) && (iter < iterMax)
         if strcmp(valueType,"true value") %真値の場合は，配列を補間して求める→時間かかるので，伝搬して求める
            xvsc = interp1(spacecraft.t, spacecraft.state.', t + dtlt_length/constant.lightSpeed, 'spline','extrap').';
         else %推定値の場合は現在時刻から伝搬する
            xvsc = Spacecraft.timeUpdate_sc(scAtT,constant.sunMu, dtlt_length/constant.lightSpeed, min(time.simDt*10,200)); %そこまで精度はいらないのでタイムステップを大きく取っている
         end         % RLT項を含まない伝搬遅延
         dtlt_length_New = norm(eAtT(1:3) + gsAtT(1:3) - xvsc(1:3) );
         corr = abs(dtlt_length_New - dtlt_length);
         dtlt_length = dtlt_length_New; 
         iter = iter + 1;
     end
     opn_t    = t + dtlt_length/constant.lightSpeed;
     if strcmp(valueType,"true value") %真値の場合は，配列を補間して求める→時間かかるので伝搬して求める
            opn_state = interp1(spacecraft.t, spacecraft.state.', opn_t, 'spline','extrap').';
     else %推定値の場合は現在時刻から伝搬する
            opn_state = Spacecraft.timeUpdate_sc(scAtT,constant.sunMu, dtlt_length/constant.lightSpeed, min(time.simDt*10,200)); %そこまで精度はいらないのでタイムステップを大きく取っている
     end         % RLT項を含まない伝搬遅延
     
     
     dtlt = dtlt_length/constant.lightSpeed;
end