% 概要: time.list分の探査機の位置を伝搬して求めていく
function calcOrbitTwoBody(obj,dynamicsError)
    mu = obj.mu;
    % 初期値
    state0 = obj.state;
    %まずはobjを初期化
    obj.state = zeros(6,length(obj.t));
    % 初期値の入力
    obj.state(:,1) = reshape(state0,[6,1]);
    %順伝搬
    xv = obj.state(:,1);
    for j = 2:length(obj.t)
        timestep = obj.t(j) - obj.t(j-1);
        a_pertubation = [dynamicsError * randn ; dynamicsError * randn; dynamicsError * randn];
        k1 = CelestialBody.twobody(xv,mu,a_pertubation);
        k2 = CelestialBody.twobody(xv+0.5*timestep*k1,mu,a_pertubation);
        k3 = CelestialBody.twobody(xv+0.5*timestep*k2,mu,a_pertubation);
        k4 = CelestialBody.twobody(xv+timestep*k3,mu,a_pertubation);
        xv = xv + timestep/6*(k1+2*k2+2*k3+k4); 
        obj.state(:,j) = xv;
    end
end
