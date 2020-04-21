% 概要: time.list分の探査機の位置を伝搬して求めていく
function calcOrbitTwoBody_pertubation(obj,time,constant, error)
    mu = constant.sunMu;
    %まずはobjを初期化
    obj.t = time.list;
    obj.pos = zeros(3,length(time.list));
    obj.vel = zeros(3,length(time.list));
    % 初期値の入力
    obj.pos(:,1) = reshape(obj.pos0,[3,1]);
    obj.vel(:,1) = reshape(obj.vel0,[3,1]); 
    %順伝搬
    xv = [obj.pos(:,1);obj.vel(:,1)] ;
    timestep = time.simDt;
    for j = 2:length(time.list)
        a_pertubation = [error.dynamics * randn ; error.dynamics * randn; error.dynamics * randn];
        k1 = obj.twobody(xv,mu,a_pertubation);
        k2 = obj.twobody(xv+0.5*timestep*k1,mu,a_pertubation);
        k3 = obj.twobody(xv+0.5*timestep*k2,mu,a_pertubation);
        k4 = obj.twobody(xv+timestep*k3,mu,a_pertubation);
        xv = xv + timestep/6*(k1+2*k2+2*k3+k4); 
        obj.pos(:,j) = xv(1:3);
        obj.vel(:,j) = xv(4:6);
    end
end



