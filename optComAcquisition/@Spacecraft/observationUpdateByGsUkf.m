function observationUpdateByGsUkf(obj,gsTrue,earth,constant,ukf)
    % 何度目の観測か
    dr_counter = gsTrue.dr_counter;
    %% 必要な変数を取得
    % uplinkを送信した時の地球+地上局の位置・速度
    xv_ut = earth.state_ut(:,dr_counter) +  gsTrue.state_ut(:,dr_counter);
    % downlinkを受信した時刻の地球+地上局の位置・速度
    xv_dr = earth.state_dr(:,dr_counter) +  gsTrue.state_dr(:,dr_counter);
    % 宇宙機が受信してから送信するまでの時間
    dtAtSc = gsTrue.durationAtSc(dr_counter);

    % シグマ点列を取得
    x_sp   = obj.x_sp_dt;
    % シグマ点列の重み平均を取得
    x_mean = obj.X_dt;
    %% 観測値を取得 
    Y = [gsTrue.directionObserved_dr(:,dr_counter);... % 測角
         gsTrue.scAccel_dr(:,dr_counter);...       % 加速度計
         gsTrue.lengthObserved_dr(dr_counter);...
         gsTrue.length2wObserved_dr(dr_counter)];    % 測距(1way)
    %% シグマ点列での観測ベクトルを取得
    y_sp = zeros(length(Y),size(x_sp,2));
    for i = 1:size(x_sp,2)
        x_spi = x_sp(:,i);
        y_sp(:,i) = Spacecraft.calcG_dr(x_spi,xv_ut,xv_dr,dtAtSc,constant);
        if i == 1
            y_mean = ukf.w0_m * y_sp(:,i);
        else
            y_mean = y_mean + ukf.wi_m * y_sp(:,i);
        end
    end
    %% 共分散行列及び相互共分散行列
    for j = 1:size(x_sp,2)
        if j == 1
            Pvv = ukf.w0_c * (y_sp(:,j) - y_mean) * (y_sp(:,j) - y_mean).';
            Pxy = ukf.w0_c * (x_sp(:,j) - x_mean) * (y_sp(:,j) - y_mean).';
        else
            Pvv = Pvv + ukf.wi_c * (y_sp(:,j) - y_mean) * (y_sp(:,j) - y_mean).';
            Pxy = Pxy + ukf.wi_c * (x_sp(:,j) - x_mean) * (y_sp(:,j) - y_mean).';           
        end
    end
       
    % Rを取得．QDセンサーの精度を反映
    R = obj.R2wGs;
    R(1,1) = gsTrue.directionAccuracy_dr(dr_counter)^2;
    R(2,2) = R(1,1);
    obj.R2wGs = R;
    
    Pvv = Pvv + R;
    
    %% 観測残差及び残差検定
    v = Y - y_mean;
    
    % ここで，方位角の測角の部分で2piの不連続を回避するためにの計算をする．
    for v_i =  1 %1番目は受信側の方位角の測角，6番目は送信側の方位角の測角
        v(v_i) = mod(v(v_i) + pi, 2*pi) - pi;
    end    
   for k = length(v):-1:1
        if ukf.sigmaN < abs(v(k))/sqrt(Pvv(k,k))
            % 観測を棄却する
%             v(k) = 0;
            v(k) = [];
            Pvv(k,:) = [];
            Pvv(:,k) = [];
            Pxy(:,k) = [];
        end
   end
    
   % 有効な観測がない時は例外処理
   if isempty(v)
       obj.X_dt = x_mean;
       obj.P_dt = obj.P_dt;
   else
       %% カルマンゲインの計算
       K = Pxy /Pvv;
       %% 状態ベクトルと共分散行列の観測更新
       obj.X_dt = x_mean + K * v;
       obj.P_dt = obj.P_dt - K *Pxy.';
   end   
       
    % デバッグ用に記録
    obj.y = v;  



end