function sum_part_y = fif_fun_2( theta,front_constant,A,B,b,aerfa,ny)

sum_part_y=(1-exp((ny).*sin(theta).*(front_constant.*((-b).*sin(aerfa).*(sqrt(1-A.^2)).*cos(B)+b.*cos(aerfa).*A))))./(1-exp(sin(theta).*(front_constant.*((-b).*sin(aerfa).*(sqrt(1-A.^2)).*cos(B)+b.*cos(aerfa).*A)))).*conj((1-exp((ny).*sin(theta).*(front_constant.*((-b).*sin(aerfa).*(sqrt(1-A.^2)).*cos(B)+b.*cos(aerfa).*A))))./(1-exp(sin(theta).*(front_constant.*((-b).*sin(aerfa).*(sqrt(1-A.^2)).*cos(B)+b.*cos(aerfa).*A)))));
sum_part_y(isnan(sum_part_y))=(ny)^2;
sum_part_y(isinf(sum_part_y))=(ny)^2;
end



