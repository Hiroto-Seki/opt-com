function [Yv,YStarv,Hm,Rm] = alignReqInfo4Est(Y,Y_star,H,R,type,reqList)

% ここで，そもそも選択によらず観測できないものは除く
switch type
    case "1u"
        reqList.length2w_ur  = 0; %この関数の外には影響ないので，書き換えてOK
        reqList.direction_dr = 0;
        reqList.length1w_dr  = 0;
        reqList.length2w_dr  = 0;
    case "2u"
        reqList.length1w_dr  = 0;
        reqList.length2w_dr  = 0;
    case "2d"
        reqList.length1w_ur  = 0;
        reqList.length2w_ur  = 0;
end

% 初期化
Yv     = [];
YStarv = [];
Hm     = [];
Rv     = [];

% 使用する観測のみを付け足していく

if reqList.direction_ur == 1 % 要素は2つ分
    Yv     = [Yv;Y.direction_ur];
    YStarv = [YStarv;Y_star.azm_ur;Y_star.elv_ur];
    Hm     = [Hm;H.azm_ur;H.elv_ur];
    Rv     = [Rv;ones(2,1) * R.direction_ur];
end

if reqList.direction_ut == 1 % 要素は2つ分
    Yv     = [Yv;Y.direction_ut];
    YStarv = [YStarv;Y_star.azm_ut;Y_star.elv_ut];
    Hm     = [Hm;H.azm_ut;H.elv_ut];
    Rv     = [Rv; ones(2,1) *R.direction_ut];
end

if reqList.accel_ur == 1 % 要素は3つ分
    Yv     = [Yv;Y.accel_ur];
    YStarv = [YStarv;Y_star.accel_ur];
    Hm     = [Hm;H.accel_ur];
    Rv     = [Rv; ones(3,1) * R.accel_ur];
end

if reqList.length1w_ur == 1 % 要素は1つ分
    Yv     = [Yv;Y.length1w_ur];
    YStarv = [YStarv;Y_star.length1w_ur];
    Hm     = [Hm;H.length1w_ur];
    Rv     = [Rv; ones(1,1) * R.length1w_ur];
end

if reqList.length2w_ur == 1 % 要素は1つ分
    Yv     = [Yv;Y.length2w_ur];
    YStarv = [YStarv;Y_star.length2w_ur];
    Hm     = [Hm;H.length2w_ur];
    Rv     = [Rv; ones(1,1) * R.length2w_ur];
end

if reqList.direction_dr == 1 % 要素は2つ分
    Yv     = [Yv;Y.direction_dr];
    YStarv = [YStarv;Y_star.azm_dr;Y_star.elv_dr];
    Hm     = [Hm;H.azm_dr;H.elv_dr];
    Rv     = [Rv; ones(2,1) * R.direction_dr];
end

if reqList.length1w_dr == 1 % 要素は1つ分
    Yv     = [Yv;Y.length1w_dr];
    YStarv = [YStarv;Y_star.length1w_dr];
    Hm     = [Hm;H.length1w_dr];
    Rv     = [Rv; ones(1,1) * R.length1w_dr];
end

if reqList.length2w_dr == 1 % 要素は1つ分
    Yv     = [Yv;Y.length2w_dr];
    YStarv = [YStarv;Y_star.length2w_dr];
    Hm     = [Hm;H.length2w_dr];
    Rv     = [Rv; ones(1,1) * R.length2w_dr];
end


% Rv(ベクトル)を対角行列に変換する
lenR = length(Rv);
Rm   = zeros(lenR);
for i = 1:lenR
    Rm(i,i) = Rv(i);
end

end
