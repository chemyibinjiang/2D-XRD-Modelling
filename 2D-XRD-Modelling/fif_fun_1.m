function sum_part_x = fif_fun_1( theta,front_constant,A,a,nx )

sum_part_x=(1-exp((nx).*sin(theta).*front_constant.*a.*A))./(1-exp(sin(theta).*front_constant.*a.*A)).*conj((1-exp((nx).*sin(theta).*front_constant.*a.*A))./(1-exp(sin(theta).*front_constant.*a.*A)));
sum_part_x(isnan(sum_part_x))=(nx)^2;
sum_part_x(isinf(sum_part_x))=(nx)^2;
end

