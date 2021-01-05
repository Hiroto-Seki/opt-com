
function observationUpdateByScUkf(obj,scTrue,earth, gsTrue,constant,type,ukf,time)
    % 何度目の観測か
    ur_counter = scTrue.ur_counter;
    ur2w_counter = scTrue.ur2w_counter;
    % シグマ点列を取得
    x_sp   = obj.x_sp;
    % シグマ点列の重み平均を取得
    x_mean = obj.X;
    
    % 観測方程式に出てくるもの
    xve_ut = earth.state_ut(:,ur_counter);
    xvg_ut = gsTrue.state_ut(:,ur_counter);
    if type ==2
        xve_dr = earth.state_dr(:,ur2w_counter);
        xvg_dr = gsTrue.state_dr(:,ur2w_counter);
        dtAtGs = gsTrue.t_ut(ur_counter) - gsTrue.t_dr(ur2w_counter);
        dt2w  = scTrue.t_ur(ur_counter) - scTrue.t_dt(ur2w_counter);        
    end
    %% 観測値を取得
    if type == 1 %1wayの観測
        Y = [scTrue.directionObserved_ur(:,ur_counter);... % 測角
            scTrue.accelObseved_ur(:,ur_counter);...       % 加速度計
            scTrue.transDirection_ur(:,ur_counter);...     % uplink方向
            scTrue.lengthObserved_ur(ur_counter)];         % 測距 
    else %2wayの観測
        Y = [scTrue.directionObserved_ur(:,ur_counter);... % 測角
            scTrue.accelObseved_ur(:,ur_counter);...       % 加速度計
            scTrue.transDirection_ur(:,ur_counter);...     % uplink方向
            scTrue.lengthObserved_ur(ur_counter);          % 測距1way
            scTrue.length2wObserved_ur(ur_counter);];      % 測距2way
    end
    
    %% シグマ点列による観測ベクトルの計算
    y_sp = zeros(length(Y),size(x_sp,2));
    for i = 1:size(x_sp,2)
        x_spi = x_sp(:,i);
        % 観測によって場合分け
        if type == 1
            y_sp(:,i) = Spacecraft.calcG1w_ur(x_spi,xve_ut,xvg_ut,constant);
        else
            y_sp(:,i) = Spacecraft.calcG2w_ur(x_spi,xve_ut,xvg_ut,xve_dr,xvg_dr,...
                                        dtAtGs, dt2w,constant,time);
        end
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
    
    % 誤差共分散に観測誤差共分散の項を追加
    if type == 1 %1wayの観測
        R = obj.R1wSc;
    else %2wayの観測
        R = obj.R2wSc;
    end
    R(1,1) = scTrue.directionAccuracy_ur(ur_counter)^2;
    R(2,2) = R(1,1); 
    Pvv = Pvv + R;
    

    %% 観測残差及び残差検定
    v = Y - y_mean;
    % ここで，方位角の測角の部分で2piの不連続を回避するためにの計算をする．
    for v_i = [1,6] %1番目は受信側の方位角の測角，6番目は送信側の方位角の測角
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
    
    %有効な観測がない時は例外処理
   if isempty(v)
       obj.X = x_mean;
       obj.P = obj.P;
   else
        %% カルマンゲインの計算
        K = Pxy /Pvv;

        %% 状態ベクトルと共分散行列の観測更新
        obj.X = x_mean + K * v;
        obj.P = obj.P - K *Pxy.';
   end      
    % デバッグ用に記録
    obj.y = v;

end