% EKFのための観測式の微分を求める
function H_childa = delGdelX(X_bar,xve,xvg,constant)
    xvs = X_bar(2:7);
    c = constant.lightSpeed;
    H_childa = zeros(3,7);
    % D:距離
    D = ((xve(1) + xvg(1) - xvs(1))^2 + (xve(2) + xvg(2) - xvs(2))^2 +(xve(3) + xvg(3) - xvs(3))^2)^0.5;
    % Dのxs微分
    delDdelXs = (xvs(1) - xve(1) - xvg(1))/D;
    delDdelYs = (xvs(2) - xve(2) - xvg(2))/D;
    delDdelZs = (xvs(3) - xve(3) - xvg(3))/D;
    
    %% Tの微分
    % TのdeltaT微分
    H_childa(1,1) = 1;
    % Tのxs微分
    H_childa(1,2) = delDdelXs/c;
    % Tのys微分
    H_childa(1,3) = delDdelYs/c;
    % Tのzs微分
    H_childa(1,4) = delDdelZs/c;
    % Tのus微分
    H_childa(1,5) = 0;
    % Tのus微分
    H_childa(1,6) = 0;
    % Tのws微分
    H_childa(1,7) = 0;
    
    %% thetaの微分
    %共通で出てくる項はまとめる
    % tan(theta)の分母->t_d (theta_denominator)
    t_d = c*(xve(1) + xvg(1) - xvs(1)) + D*xvs(4);
    % tan(theta)の分詞->t_n (theta_numerator)
    t_n = c*(xve(2) + xvg(2) - xvs(2)) + D*xvs(5);
    
    % thetaのdeltaT微分
    H_childa(2,1) = 0;
    % thetaのxs微分
    H_childa(2,2) = (t_n^2 + t_d^2 ) * (      (xvs(5) * delDdelXs) * t_d - t_n*(-c + xvs(4)*delDdelXs));
    % thetaのys微分
    H_childa(2,3) = (t_n^2 + t_d^2 ) * ( (-c + xvs(5) * delDdelYs) * t_d - t_n*    ( xvs(4)*delDdelYs));
    % thetaのzs微分
    H_childa(2,4) = (t_n^2 + t_d^2 ) * (      (xvs(5) * delDdelZs) * t_d - t_n*    ( xvs(4)*delDdelZs));
    % thetaのus微分
    H_childa(2,5) = (t_n^2 + t_d^2 ) * (                                 - t_n*    ( D               ));
    % thetaのus微分
    H_childa(2,6) = (t_n^2 + t_d^2 ) * (                       (D) * t_d                              );
    % thetaのws微分
    H_childa(2,7) = 0;    

    %% phiの微分
    % tan(phi)の微分
    % tan(phi)の分母->p_d (phi_denominator)
    % 分母のx側の項 ->p_d_x
    p_d_x = c*(xve(1) + xvg(1) - xvs(1)) + D*xvs(4);
    % 分母のy側の項 ->p_d_y
    p_d_y = c*(xve(2) + xvg(2) - xvs(2)) + D*xvs(5);
    p_d = (p_d_x^2 + p_d_y^2)^0.5;
    % tan(theta)の分詞->p_n (phi_numerator)
    p_n = c*(xve(3) + xvg(3) - xvs(3)) + D*xvs(6);
    
    % phiのdeltaT微分
    H_childa(3,1) = 0;
    % phiのxs微分
    H_childa(3,2) = (p_n^2 + p_d^2) * (     (xvs(6) * delDdelXs) * p_d - p_n/p_d *(p_d_x*(-c +xvs(4)*delDdelXs) + p_d_y*(     xvs(5)*delDdelXs)));
    % phiのys微分
    H_childa(3,3) = (p_n^2 + p_d^2) * (     (xvs(6) * delDdelYs) * p_d - p_n/p_d *(p_d_x*(   +xvs(4)*delDdelYs) + p_d_y*( -c +xvs(5)*delDdelYs)));
    % phiのzs微分
    H_childa(3,4) = (p_n^2 + p_d^2) * ( (-c +xvs(6) * delDdelZs) * p_d - p_n/p_d *(p_d_x*(   +xvs(4)*delDdelZs) + p_d_y*(     xvs(5)*delDdelZs)));
    % phiのus微分
    H_childa(3,5) = (p_n^2 + p_d^2) * (                                - p_n/p_d *(p_d_x*D)                                                     );
    % phiのus微分 (ここから)
    H_childa(3,6) = (p_n^2 + p_d^2) * (                                - p_n/p_d *(p_d_y*D)                                                     );
    % phiのws微分
    H_childa(3,7) = (p_n^2 + p_d^2) * (                     ( D ) * p_d                                                                         ); 
    
end