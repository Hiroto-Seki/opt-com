% “ñ‘Ì–â‘è
function xv = twobody(xv, mu, a_pertubation)
    x = xv(1);
    y = xv(2);
    z = xv(3);
    u = xv(4);
    v = xv(5);
    w = xv(6);

    r    = [x; y; z];
    R    = norm(r);
    a_0    = -mu*r/R^3;         
    a    = a_0 +a_pertubation;
    xv = [u;v;w;a];
end