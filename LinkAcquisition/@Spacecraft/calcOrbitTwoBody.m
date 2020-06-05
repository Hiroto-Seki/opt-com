% 概要: time.list分の探査機の位置を伝搬して求めていく
function calcOrbitTwoBody(obj,timestep,dynamicsError)
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
        a_pertubation = [dynamicsError * randn ; dynamicsError * randn; dynamicsError * randn];
        k1 = obj.twobody(xv,mu,a_pertubation);
        k2 = obj.twobody(xv+0.5*timestep*k1,mu,a_pertubation);
        k3 = obj.twobody(xv+0.5*timestep*k2,mu,a_pertubation);
        k4 = obj.twobody(xv+timestep*k3,mu,a_pertubation);
        xv = xv + timestep/6*(k1+2*k2+2*k3+k4); 
        obj.state(:,j) = xv;
    end
end
