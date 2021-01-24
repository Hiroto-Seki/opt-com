function [Yv,YStarv,Hm,Rm,obsNum,sigmaN] = alignReqInfo4Est(Y,Y_star,H,R,obsType,estType,reqList)

% ここで，そもそも選択によらず観測できないものは除く
switch obsType
    %% 1wayの宇宙機側の観測
    case "1u_Obs" 
        reqList.length2w_ur  = 0; %1wayuplinkの観測にはない
        reqList.direction_dr = 0; %1wayuplinkの観測にはない
        reqList.length1w_dr  = 0; %1wayuplinkの観測にはない
        reqList.length2w_dr  = 0; %1wayuplinkの観測にはない
    case  "1u_lowSnr"
        reqList.direction_ut = 0; %SN比が低くて観測できない
        reqList.length1w_ur  = 0; %SN比が低くて観測できない
        reqList.length2w_ur  = 0; %1wayuplinkの観測にはない
        reqList.direction_dr = 0; %1wayuplinkの観測にはない
        reqList.length1w_dr  = 0; %1wayuplinkの観測にはない
        reqList.length2w_dr  = 0; %1wayuplinkの観測にはない
    case  "1u_noObs"
        reqList.direction_ur  = 0; %視野に入らないので観測できない
        reqList.direction_ut  = 0; %視野に入らないので観測できない
        reqList.length1w_ur  = 0;  %視野に入らないので観測できない
        reqList.length2w_ur  = 0; %1wayuplinkの観測にはない
        reqList.direction_dr = 0;%1wayuplinkの観測にはない
        reqList.length1w_dr  = 0;%1wayuplinkの観測にはない
        reqList.length2w_dr  = 0;%1wayuplinkの観測にはない
    %% 2wayの宇宙機側の観測
    case  "2u_obsD_obsU"   %downlinkも観測できてuplinkも観測できた場合
        reqList.length1w_dr  = 0;
        reqList.length2w_dr  = 0;
%         reqList.length1w_ur  = 0; %1wayを使わないという選択
    case  "2u_obsD_lowU"   %downlinkは観測できてuplinkはSN比が低い場合
        reqList.direction_ut = 0; %SN比が低くて観測できない
        reqList.length1w_ur  = 0; %SN比が低くて観測できない
        reqList.length2w_ur  = 0; %SN比が低くて観測できない
        reqList.direction_dr = 0; %SN比が低くて観測できない
        reqList.length1w_dr  = 0; %宇宙機側の観測にはない
        reqList.length2w_dr  = 0; %宇宙機側の観測にはない
    case  "2u_obsD_noObsU" %downlinkは観測できてuplinkは観測できない場合
        reqList.direction_ur  = 0; %視野に入らないので観測できない
        reqList.direction_ut  = 0; %視野に入らないので観測できない
        reqList.length1w_ur  = 0;  %視野に入らないので観測できない
        reqList.length2w_ur  = 0;  %視野に入らないので観測できない
        reqList.direction_dr = 0;  %視野に入らないので観測できない
        reqList.length1w_dr  = 0;  %uplinkの観測にはない
        reqList.length2w_dr  = 0;  %uplinkの観測にはない             
    case  "2u_lowD_obsU"  %downlinkはSN比が低い(信号が受信できない)で，uplinkは観測できる  
        reqList.length2w_ur  = 0; %地上局が信号を受信できないので，downlinkとuplinkの対応が取れない
        reqList.length1w_dr  = 0; %uplinkの観測にはない
        reqList.length2w_dr  = 0; %uplinkの観測にはない
    case  "2u_lowD_lowU"   %downlinkはSN比が低い(信号が受信できない)で，uplinkもSN比が低い
        reqList.direction_ut = 0; %SN比が低くて観測できない
        reqList.length1w_ur  = 0; %SN比が低くて観測できない
        reqList.length2w_ur  = 0; %SN比が低くて観測できない
        reqList.direction_dr = 0; %SN比が低くて観測できない
        reqList.length1w_dr  = 0; %宇宙機側の観測にはない
        reqList.length2w_dr  = 0; %宇宙機側の観測にはない
    case  "2u_lowD_noObsU"  %downlinkはSN比が低い(信号が受信できない)で，uplinkは観測できない
        reqList.direction_ur  = 0; %視野に入らないので観測できない
        reqList.direction_ut  = 0; %視野に入らないので観測できない
        reqList.length1w_ur  = 0;  %視野に入らないので観測できない
        reqList.length2w_ur  = 0;  %視野に入らないので観測できない
        reqList.direction_dr = 0;  %視野に入らないので観測できない
        reqList.length1w_dr  = 0;  %uplinkの観測にはない
        reqList.length2w_dr  = 0;  %uplinkの観測にはない   
    case  "2u_noObsD_obsU"  %downlinkは観測できなくて，uplinkは観測できる 
        reqList.length2w_ur  = 0; %実質1way
        reqList.direction_dr = 0; %Downlinkが観測できていないので
        reqList.length1w_dr  = 0; %1wayuplinkの観測にはない
        reqList.length2w_dr  = 0; %1wayuplinkの観測にはない
    case  "2u_noObsD_lowU"  %downlinkは観測できなくて，uplinkはSN比が低い
        reqList.direction_ut = 0; %SN比が低くて観測できない
        reqList.length1w_ur  = 0; %SN比が低くて観測できない
        reqList.length2w_ur  = 0; %実質1way
        reqList.direction_dr = 0; %実質1way
        reqList.length1w_dr  = 0; %uplinkの観測にはない
        reqList.length2w_dr  = 0; %uplinkの観測にはない
    case  "2u_noObsD_noObsU" %downlinkは観測できなくて，uplinkも観測できない
        reqList.direction_ur  = 0; %視野に入らないので観測できない
        reqList.direction_ut  = 0; %視野に入らないので観測できない
        reqList.length1w_ur  = 0;  %視野に入らないので観測できない
        reqList.length2w_ur  = 0;  %実質1way
        reqList.direction_dr = 0;  %実質1way
        reqList.length1w_dr  = 0;  %uplinkの観測にはない
        reqList.length2w_dr  = 0;  %uplinkの観測にはない
    %% 2wayの地上局側の観測
    case "2d_obsU_obsD"   %uplinkを観測できて，downkinkも観測できる
        reqList.length1w_ur  = 0; %地上局の観測にはない
        reqList.length2w_ur  = 0; %地上局の観測にはない
    case "2d_obsU_lowD"   %uplinkを観測できて，downkinkはSN比が低い
        reqList.direction_ur  = 0; %downlinkのSN比が低いので受信できない
        reqList.direction_ut  = 0; %downlinkのSN比が低いので受信できない
        reqList.accel_ur      = 0; %downlinkのSN比が低いので受信できない
        reqList.length1w_ur  = 0;  %地上局の観測にはない
        reqList.length2w_ur  = 0;  %地上局の観測にはない
        reqList.length1w_dr  = 0;  %downlinkのSN比が低いので受信できない
        reqList.length2w_dr  = 0;  %downlinkのSN比が低いので受信できない
     case "2d_obsU_noObsD"   %uplinkを観測できて，downkinkは観測できない
        reqList.direction_ur  = 0; %downlinkが観測できないので，何も観測できない
        reqList.direction_ut  = 0; %downlinkが観測できないので，何も観測できない
        reqList.accel_ur      = 0; %downlinkが観測できないので，何も観測できない
        reqList.length1w_ur  = 0;  %地上局の観測にはない
        reqList.length2w_ur  = 0;  %地上局の観測にはない
        reqList.direction_dr = 0;  %downlinkが観測できないので，何も観測できない
        reqList.length1w_dr  = 0;  %downlinkが観測できないので，何も観測できない
        reqList.length2w_dr  = 0;  %downlinkが観測できないので，何も観測できない
     case "2d_lowU_obsD"   %uplinkのSN比が低く，downkinkは観測できる
        reqList.direction_ut  = 0; %uplinkで信号を受信できないので
        reqList.length1w_ur  = 0;  %地上局の観測にはない
        reqList.length2w_ur  = 0;  %地上局の観測にはない
        reqList.length2w_dr  = 0;  %uplinkで信号を受信できないので
     case "2d_lowU_lowD"   %uplinkのSN比が低く,downkinkはSN比が低い  
        reqList.direction_ur  = 0; %downlinkのSN比が低いので受信できない
        reqList.direction_ut  = 0; %downlinkのSN比が低いので受信できない
        reqList.accel_ur      = 0; %downlinkのSN比が低いので受信できない
        reqList.length1w_ur  = 0;  %地上局の観測にはない
        reqList.length2w_ur  = 0;  %地上局の観測にはない
        reqList.length1w_dr  = 0;  %downlinkのSN比が低いので受信できない
        reqList.length2w_dr  = 0;  %downlinkのSN比が低いので受信できない
     case "2d_lowU_noObsD"   %uplinkのSN比が低く,downkinkは観測できない 
        reqList.direction_ur  = 0; %downlinkが観測できないので，何も観測できない
        reqList.direction_ut  = 0; %downlinkが観測できないので，何も観測できない
        reqList.accel_ur      = 0; %downlinkが観測できないので，何も観測できない
        reqList.length1w_ur  = 0;  %地上局の観測にはない
        reqList.length2w_ur  = 0;  %地上局の観測にはない
        reqList.direction_dr = 0;  %downlinkが観測できないので，何も観測できない
        reqList.length1w_dr  = 0;  %downlinkが観測できないので，何も観測できない
        reqList.length2w_dr  = 0;  %downlinkが観測できないので，何も観測できない
     case "2d_noObsU_obsD"   %uplinkの観測ができず，downkinkは観測できる
        reqList.direction_ur  = 0; %uplinkで観測できていないので
        reqList.direction_ut  = 0; %uplinkで観測できていないので
        reqList.length1w_ur  = 0;  %地上局の観測にはない
        reqList.length2w_ur  = 0;  %地上局の観測にはない
        reqList.length2w_dr  = 0;  %uplinkで観測できていないので        
     case "2d_noObsU_lowD"   %uplinkの観測できず，downkinkはSN比が低い  
        reqList.direction_ur  = 0; %downlinkのSN比が低いので受信できない
        reqList.direction_ut  = 0; %downlinkのSN比が低いので受信できない
        reqList.accel_ur      = 0; %downlinkのSN比が低いので受信できない
        reqList.length1w_ur  = 0;  %地上局の観測にはない
        reqList.length2w_ur  = 0;  %地上局の観測にはない
        reqList.length1w_dr  = 0;  %downlinkのSN比が低いので受信できない
        reqList.length2w_dr  = 0;  %downlinkのSN比が低いので受信できない
     case "2d_noObsU_noObsD"   %uplinkの観測できず，downkinkは観測できない
        reqList.direction_ur  = 0; %downlinkが観測できないので，何も観測できない
        reqList.direction_ut  = 0; %downlinkが観測できないので，何も観測できない
        reqList.accel_ur      = 0; %downlinkが観測できないので，何も観測できない
        reqList.length1w_ur  = 0;  %地上局の観測にはない
        reqList.length2w_ur  = 0;  %地上局の観測にはない
        reqList.direction_dr = 0;  %downlinkが観測できないので，何も観測できない
        reqList.length1w_dr  = 0;  %downlinkが観測できないので，何も観測できない
        reqList.length2w_dr  = 0;  %downlinkが観測できないので，何も観測できない
    otherwise
        disp("observation type is not set correctly")
end


% 初期化
Yv     = [];
YStarv = [];
Hm     = [];
Rv     = [];
Rm     = []; 
sigmaN = [];

% 観測の数
obsNum = reqList.direction_ur + ...
         reqList.direction_ut + ... 
         reqList.accel_ur     + ...     
         reqList.length1w_ur  + ...
         reqList.length2w_ur  + ...
         reqList.direction_dr + ...
         reqList.length1w_dr  + ...
         reqList.length2w_dr;    

if strcmp(estType,"obsCheck")
    %観測の数を数えただけ 
else
    % 使用する観測のみを付け足していく
    if reqList.direction_ur == 1 % 要素は3つ分
        Yv     = [Yv;Y.direction_ur];
%         YStarv = [YStarv;Y_star.azm_ur;Y_star.elv_ur];
        YStarv = [YStarv;Y_star.direction_ur;];
        if strcmp(estType,"ekf")
%             Hm     = [Hm;H.azm_ur;H.elv_ur];
            Hm     = [Hm;H.direction_ur];
        end
        Rv     = [Rv;[1;1;1] * R.direction_ur * 1];
        sigmaN = [sigmaN;ones(3,1) * 3]; % [10,1]
    end

    if reqList.direction_ut == 1 % 要素は3つ分
        Yv     = [Yv;Y.direction_ut];
%         YStarv = [YStarv;Y_star.azm_ut;Y_star.elv_ut];
        YStarv = [YStarv;Y_star.direction_ut];
        if strcmp(estType,"ekf")
%             Hm     = [Hm;H.azm_ut;H.elv_ut];
            Hm     = [Hm;H.direction_ut];
        end
        Rv     = [Rv; [1;1;1] *R.direction_ut]; % [10000;1]
        sigmaN = [sigmaN;ones(3,1) * 3];
    end

    if reqList.accel_ur == 1 % 要素は3つ分
        Yv     = [Yv;Y.accel_ur];
        YStarv = [YStarv;Y_star.accel_ur];
        if strcmp(estType,"ekf")
            Hm     = [Hm;H.accel_ur];
        end
        Rv     = [Rv; ones(3,1) * R.accel_ur * 2];
        sigmaN = [sigmaN;ones(3,1) * 3];
    end

    if reqList.length1w_ur == 1 % 要素は1つ分
        Yv     = [Yv;Y.length1w_ur];
        YStarv = [YStarv;Y_star.length1w_ur];
        if strcmp(estType,"ekf")
            Hm     = [Hm;H.length1w_ur];
        end
        Rv     = [Rv; ones(1,1) * R.length1w_ur * 1];  % 1000
        sigmaN = [sigmaN;ones(1,1) * 3];
    end

    if reqList.length2w_ur == 1 % 要素は1つ分
        Yv     = [Yv;Y.length2w_ur];
        YStarv = [YStarv;Y_star.length2w_ur];
        if strcmp(estType,"ekf")
            Hm     = [Hm;H.length2w_ur];
        end
        Rv     = [Rv; ones(1,1) * R.length2w_ur * 1]; % 100
        sigmaN = [sigmaN;ones(1,1) * 3];
    end

    if reqList.direction_dr == 1 % 要素は3つ分
        Yv     = [Yv;Y.direction_dr];
%         YStarv = [YStarv;Y_star.azm_dr;Y_star.elv_dr];
        YStarv = [YStarv;Y_star.direction_dr];
        if strcmp(estType,"ekf")
%             Hm     = [Hm;H.azm_dr;H.elv_dr];
            Hm     = [Hm;H.direction_dr];
        end
        Rv     = [Rv; ones(3,1) * R.direction_dr * 2];
        sigmaN = [sigmaN;ones(3,1) * 3];
    end

    if reqList.length1w_dr == 1 % 要素は1つ分
        Yv     = [Yv;Y.length1w_dr];
        YStarv = [YStarv;Y_star.length1w_dr];
        if strcmp(estType,"ekf")
            Hm     = [Hm;H.length1w_dr];
        end
        Rv     = [Rv; ones(1,1) * R.length1w_dr * 2];
        sigmaN = [sigmaN;ones(1,1) * 3];
    end

    if reqList.length2w_dr == 1 % 要素は1つ分
        Yv     = [Yv;Y.length2w_dr];
        YStarv = [YStarv;Y_star.length2w_dr];
        if strcmp(estType,"ekf")
            Hm     = [Hm;H.length2w_dr];
        end
        Rv     = [Rv; ones(1,1) * R.length2w_dr * 2];
        sigmaN = [sigmaN;ones(1,1) * 3];
    end


    % Rv(ベクトル)を対角行列に変換する
    lenR = length(Rv);
    Rm   = zeros(lenR);
    for i = 1:lenR
        Rm(i,i) = Rv(i);
    end
end

end
