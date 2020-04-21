% 地球の時点による加速度
function xv = calcEarthRotation(xv, constant)
    x = xv(1);
    y = xv(2);
    z = xv(3);
    u = xv(4);
    v = xv(5);
    w = xv(6);
    omega = constant.earthRotation;
    phi = constant.eathAxis;
    ax = -omega^2 * x;
    ay = -omega^2 * (y*cos(phi) - z*sin(phi)) * cos(phi);
    az = omega^2 *  (y*cos(phi) - z*sin(phi)) * sin(phi);
    xv = [u;v;w;ax;ay;az];
end