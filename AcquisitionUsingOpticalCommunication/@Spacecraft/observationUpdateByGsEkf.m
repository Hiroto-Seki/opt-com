%---- 中間パラメーター -----
% X_star : アプリオリな推定値
% P_bar  : 誤差共分散行列(観測前)
% Y      : 観測値
% Y_star : アプリオリな推定値に対応する観測値
% y      : Y - Y_star
% H      : 観測方程式を状態量(=アプリオリな推定値)で微分したもの
% R      : 観測誤差共分散行列
% K      : カルマンゲイン
% x      : X - X_star
% X      : 更新した推定値
% P      : 更新した誤差共分散行列

function observationUpdateByGsEkf(obj,gsTrue,earth,constant,ekf)
    % 何度目の観測か
    dr_counter = gsTrue.dr_counter;
    %% 必要な変数を取得
    % uplinkを送信した時の地球+地上局の位置・速度
    xv_ut = earth.state_ut(:,dr_counter) +  gsTrue.state_ut(:,dr_counter);
    % downlinkを受信した時刻の地球+地上局の位置・速度
    xv_dr = earth.state_dr(:,dr_counter) +  gsTrue.state_dr(:,dr_counter);
    % 宇宙機が受信してから送信するまでの時間
    dtAtSc = gsTrue.durationAtSc(dr_counter);

    % X_starとP_barを取得
    X_star = obj.X_dt;
    P_bar  = obj.P_dt;
    % Rを取得．QDセンサーの精度を反映
    R = obj.R2wGs;
    R(1,1) = gsTrue.directionAccuracy_dr(dr_counter)^2;
    R(2,2) = R(1,1);
    obj.R2wGs = R;
    %% Yを取得. 
    Y = [gsTrue.directionObserved_dr(:,dr_counter);... % 測角
         gsTrue.scAccel_dr(:,dr_counter);...       % 加速度計
         gsTrue.lengthObserved_dr(dr_counter);...
         gsTrue.length2wObserved_dr(dr_counter)];    % 測距(1way)
    %% Y_starを計算
    Y_star = Spacecraft.calcG_dr(X_star,xv_ut,xv_dr,dtAtSc,constant);
    y = Y - Y_star; 
    % 角度の不連続性を解消
    for y_i = 1 %1番目は受信側の方位角の測角
        y(y_i) = mod(y(y_i) + pi, 2*pi) - pi;
    end
    obj.y    = y;

    
    % Hを計算
    H = Spacecraft.delGdelX_dr(X_star,xv_ut,xv_dr,dtAtSc,constant);
    
    % 観測残差及び残差検定
    for k = length(y):-1:1
        S = (H*P_bar*H.' + R);
        if ekf.sigmaN < abs(y(k))/sqrt(S(k,k)) %3シグマに設定している
            y(k) = [];
            H(k,:) = [];
            R(k,:) = [];
            R(:,k) = [];
        end  
    end
    
   %% 有効な観測がない時は例外処理をする
   if isempty(y)
       X = X_star;
       P = P_bar;
   else
       % Kを計算
       K = P_bar * H.'/(H*P_bar*H.' + R);
       % XとPを計算
       X = X_star + K * y;
       P = (eye(7) - K * H)*P_bar;
   end
    obj.X_dt = X;
    obj.P_dt = P;


end