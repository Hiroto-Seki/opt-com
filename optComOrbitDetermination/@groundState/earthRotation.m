function [pos,vel] = earthRotation(pos0, t, constant)
    x0 = pos0(1);
    y0 = pos0(2);
    z0 = pos0(3);
    omega = constant.earthRotation;
    phi = constant.eathAxis;
    x = x0*cos(omega * t) - (y0*cos(phi) - z0*sin(phi))*sin(omega * t);
    y = ((y0*cos(phi) - z0*sin(phi))*cos(omega * t) + x0*sin(omega * t))*cos(phi) + (z0*cos(phi) + y0*sin(phi)*sin(phi));
    z =  (z0*cos(phi) + y0*sin(phi)*cos(phi)) - ((y0*cos(phi) - z0*sin(phi))*cos(omega * t) + x0*sin(omega * t))*sin(phi);
    
    u = -omega * x0*sin(omega * t) - (y0*cos(phi) - z0*sin(phi))* omega * cos(omega * t);
    v = ((y0*cos(phi) - z0*sin(phi))* (-omega) * sin(omega * t) + x0* omega * cos(omega * t))*cos(phi);
    w = - ((y0*cos(phi) - z0*sin(phi))* (-omega) *sin(omega * t) + x0* (omega) * cos(omega * t))*sin(phi);
    pos = [x;y;z];
    vel = [u;v;w];
end