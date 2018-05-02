function v_s = standCart2Sph(v_c)
% convert the vector from cartesian to spherical coordinates
% x, y, z -> r, theta, phi (where theta = 0 to 2*pi)
v_s(1) = norm(v_c); %r
v_s(2) = atan2(v_c(2),v_c(1)); %theta
v_s(3) = acos(v_c(3)/v_s(1));%phi

       