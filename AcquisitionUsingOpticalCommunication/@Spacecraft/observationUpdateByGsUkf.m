function observationUpdateByGsUkf(obj,gsTrue,earth,constant,ukf,scTrue)
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
    
    %% Rを取得
    obj.R.direction_ur = gsTrue.scRecAngleAccuracy_dr(dr_counter)^2;
    obj.R.direction_dr = gsTrue.directionAccuracy_dr(dr_counter)^2;
    obj.R.direction_ut = gsTrue.transUpAngleAccuracy_dr(dr_counter)^2;
    
    %% 観測値を取得 
    Y.direction_ur = gsTrue.scRecAngle_dr(:,dr_counter); % 測角
    Y.direction_ut = gsTrue.transUpAngle_dr(:,dr_counter);     % uplink方向
    Y.accel_ur     = gsTrue.scAccel_dr(:,dr_counter);      % 加速度計
    Y.direction_dr = gsTrue.directionObserved_dr(:,dr_counter);         % uplinkされる，地上局での観測
    Y.length1w_dr  = gsTrue.lengthObserved_dr(dr_counter);
    Y.length2w_dr  = gsTrue.length2wObserved_dr(dr_counter);
    %% シグマ点列での観測ベクトルを取得
    % 初期化
    y_mean.azm_ur = 0;
    y_mean.elv_ur = 0;
    y_mean.azm_ut = 0;
    y_mean.elv_ut = 0;
    y_mean.accel_ur = zeros(3,1);
    y_mean.azm_dr = 0;
    y_mean.elv_dr = 0;
    y_mean.length1w_dr = 0;
    y_mean.length2w_dr = 0;
    
    %% 観測できているかの場合わけ
   switch scTrue.ur_observability(dr_counter)
        case 1 %uplinkが観測されていない
            switch gsTrue.dr_observability(dr_counter)
                case 1
                    obsType = "2d_noObsU_noObsD";
                case 2
                    obsType = "2d_noObsU_lowD";
                case 3
                    obsType = "2d_noObsU_obsD";
                otherwise
                    disp("obsType is not set correctly")
            end  
        case 2 %uplinkのSNRが低い
            switch gsTrue.dr_observability(dr_counter)
                case 1
                    obsType = "2d_lowU_noObsD";
                case 2
                    obsType = "2d_lowU_lowD";
                case 3
                    obsType = "2d_lowU_obsD";
                otherwise
                    disp("obsType is not set correctly")
            end  
        case 3 %uplinkが観測されている
            switch gsTrue.dr_observability(dr_counter)
                case 1
                    obsType = "2d_obsU_noObsD";
                case 2
                    obsType = "2d_obsU_lowD";
                case 3
                    obsType = "2d_obsU_obsD";
                otherwise
                    disp("obsType is not set correctly")
            end
        otherwise
            disp("obsType is not set correctly")
   end 
    
   % 有効な観測があるかどうかで場合わけ
   [~,~,~,~,obsNum] = Spacecraft.alignReqInfo4Est([],[],[],[],obsType,"obsCheck",obj.useObs); %有効な観測があるかだけを確認している
   if obsNum == 0
       obj.X_dt = x_mean;
       obj.P_dt = obj.P_dt;
       
   else
        for i = 1:size(x_sp,2)
            x_spi = x_sp(:,i);
            y_sp(i) = Spacecraft.calcG_dr(x_spi,xv_ut,xv_dr,dtAtSc,constant);
            [~,y_spM(:,i),~,~] = Spacecraft.alignReqInfo4Est(Y,y_sp(i),[],obj.R,obsType,"ukf",obj.useObs);
            if i == 1
                wm = ukf.w0_m;
            else
                wm = ukf.wi_m;
            end
            y_mean.azm_ur = y_mean.azm_ur + wm * y_sp(i).azm_ur;
            y_mean.elv_ur = y_mean.elv_ur + wm * y_sp(i).elv_ur;
            y_mean.azm_ut = y_mean.azm_ut + wm * y_sp(i).azm_ut;
            y_mean.elv_ut = y_mean.elv_ut + wm * y_sp(i).elv_ut;
            y_mean.accel_ur = y_mean.accel_ur  + wm * y_sp(i).accel_ur;
            y_mean.azm_dr = y_mean.azm_dr + wm * y_sp(i).azm_dr;
            y_mean.elv_dr = y_mean.elv_dr + wm * y_sp(i).elv_dr; 
            y_mean.length1w_dr = y_mean.length1w_dr + wm * y_sp(i).length1w_dr;
            y_mean.length2w_dr = y_mean.length2w_dr + wm * y_sp(i).length2w_dr;
        end

        % Y, Y_meanをベクトル化，Rを行列に
        [Yv,y_meanV,~,Rm] = Spacecraft.alignReqInfo4Est(Y,y_mean,[],obj.R,obsType,"ukf",obj.useObs);

        %% 共分散行列及び相互共分散行列
        for j = 1:size(x_sp,2)
            if j == 1
                Pvv = ukf.w0_c * (y_spM(:,j) - y_meanV) * (y_spM(:,j) - y_meanV).';
                Pxy = ukf.w0_c * (x_sp(:,j) - x_mean) * (y_spM(:,j) - y_meanV).';
            else
                Pvv = Pvv + ukf.wi_c * (y_spM(:,j) - y_meanV) * (y_spM(:,j) - y_meanV).';
                Pxy = Pxy + ukf.wi_c * (x_sp(:,j) - x_mean) * (y_spM(:,j) - y_meanV).';           
            end
        end


        Pvv = Pvv + Rm;

        %% 観測残差及び残差検定
        v = Yv - y_meanV;

        % デバッグ用に記録
        obj.y = v; 

       for k = length(v):-1:1
            if ukf.sigmaN < abs(v(k))/sqrt(Pvv(k,k))
                % 観測を棄却する
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
 

   end

end