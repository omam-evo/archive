function e = e_vartheta_a_b(theta, a, b)
   e = (exp(-1/2*norminv(theta).^2)./theta/sqrt(2*pi)).^a .* (-norminv(theta)).^b;
end