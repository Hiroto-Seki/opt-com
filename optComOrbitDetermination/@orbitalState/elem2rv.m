% basiliskのorbitalMotion.elem2rv_parabを参考に作成
% angleのところがあっているか不安なので，天体と軌道の力学で確認する
function elem2rv(obj)

%     Function: elem2rv
%     Purpose: Translates the orbit elements
%             a   - semi-major axis           (km)
%             e   - eccentricity
%             i   - inclination               (rad)
%             AN  - ascending node            (rad)
%             AP  - argument of periapses     (rad)
%             f   - true anomaly angle        (rad)
%     to the inertial Cartesian position and velocity vectors.
%     The attracting body is specified through the supplied
%     gravitational constant mu (units of km^3/s^2).
% 
%     Inputs:
%         obj = orbital elements
%     Outputs:
%         obj.pos = position vector
%         obj.vel = velocity vector

    obj.pos = zeros(3,1);
    obj.vel = zeros(3,1);
%     # Calculate the semilatus rectum and the radius #
    p = obj.a * (1.0 - obj.e * obj.e);
    r = p / (1.0 + obj.e * cos(obj.f));
    theta = obj.AP + obj.f; %true latitude angle
    obj.pos(1) = r * (...
        cos(theta) * cos(obj.AN) - cos(obj.i) * sin(theta) * sin(...
            obj.AN));
    obj.pos(2) = r * (...
        cos(theta) * sin(obj.AN) + cos(obj.i) * sin(theta) * cos(...
            obj.AN));
    obj.pos(3) = r * (sin(theta) * sin(obj.i));

    h = sqrt(obj.mu * p);
    obj.vel(1) = -obj.mu / h * (cos(obj.AN) * (obj.e * sin(obj.AP) + sin(theta)) + ...
                         cos(obj.i) * ( ...
                             obj.e * cos(obj.AP) + cos(theta)) * sin(obj.AN));
    obj.vel(2) = -obj.mu / h * (sin(obj.AN) * (obj.e * sin(obj.AP) + sin(theta)) - ...
                         cos(obj.i) * (obj.e * cos(obj.AP) + cos(theta)) * ...
                         cos(obj.AN)); ...
    obj.vel(3) = obj.mu / h * (obj.e * cos(obj.AP) + cos(theta)) * sin(obj.i);
    

end