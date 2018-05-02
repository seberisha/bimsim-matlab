function v_c = standSph2Cart(v_s)
% convert the vector from cartesian to spherical coordinates
% r, theta, phi -> x, y, z (where theta = 0 to 2*pi)
v_c(1) = v_s(1) * cos(v_s(2)) * sin(v_s(3));
v_c(2) = v_s(1)* sin(v_s(2)) * sin(v_s(3));
v_c(3) = v_s(1) * cos(v_s(3));