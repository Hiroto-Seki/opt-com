function getTrajectoryEarthRotation(obj,time,constant)
% èâä˙âªÇ∑ÇÈ
    obj.t = time.list;
    obj.pos = zeros(3,length(time.list));
    obj.vel = zeros(3,length(time.list));
    for i = 1:length(time.list)
        t = time.list(i);
        [gsPos,gsVel] = obj.earthRotation(obj.pos0, t, constant);
        obj.pos(:,i) = gsPos;
        obj.vel(:,i) = gsVel;
    end
end

