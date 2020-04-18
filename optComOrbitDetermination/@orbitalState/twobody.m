% “ñ‘Ì–â‘è
function xv = twobody(xv, mu)
    x = xv(1);
    y = xv(2);
    z = xv(3);
    u = xv(4);
    v = xv(5);
    w = xv(6);

    r    = [x; y; z];
    R    = norm(r);
    a    = -mu*r/R^3;         
   % a_pertubation = [error * randn ; error * randn; error * randn] ;
%     a    = a_0 +a_pertubation;
    xv = [u;v;w;a];
end