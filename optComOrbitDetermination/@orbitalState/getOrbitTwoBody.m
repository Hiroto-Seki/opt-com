% 概要: time.list分の探査機の位置を伝搬して求めていく
function getOrbitTwoBody(obj,time,constant)
    mu = constant.sunMu;
    %まずはobjを初期化
    obj.t = time.list;
    obj.pos = zeros(3,length(time.list));
    obj.vel = zeros(3,length(time.list));
    % 初期値の入力
    obj.pos(:,time.refId) = reshape(obj.pos0,[3,1]);
    obj.vel(:,time.refId) = reshape(obj.vel0,[3,1]);
    
    %順伝搬
    xv = [obj.pos(:,time.refId);obj.vel(:,time.refId)] ;
    timestep = time.simDt;
    for j = time.refId+1:length(time.list)
        k1 = obj.twobody(xv,mu);
        k2 = obj.twobody(xv+0.5*timestep*k1,mu);
        k3 = obj.twobody(xv+0.5*timestep*k2,mu);
        k4 = obj.twobody(xv+timestep*k3,mu);
        xv = xv + timestep/6*(k1+2*k2+2*k3+k4); 
        obj.pos(:,j) = xv(1:3);
        obj.vel(:,j) = xv(4:6);
    end
    %逆伝搬
    xv =[obj.pos(:,time.refId);obj.vel(:,time.refId)];
    timestep = -time.simDt;
    for j = time.refId-1 : -1:  1
        k1 = obj.twobody(xv,mu);
        k2 = obj.twobody(xv+0.5*timestep*k1,mu);
        k3 = obj.twobody(xv+0.5*timestep*k2,mu);
        k4 = obj.twobody(xv+timestep*k3,mu);
        xv = xv + timestep/6*(k1+2*k2+2*k3+k4); 
        obj.pos(:,j) = xv(1:3);
        obj.vel(:,j) = xv(4:6);
    end
end



