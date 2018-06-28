function sum_part_z= fif_fun_4( theta,front_constant,A,B,a,b,nz,aerfa,c_2,k_a,k_b )
SIA=(sqrt(1-A.^2));
siAsiB=SIA.*sin(B);
sum_part_z=(1-exp((nz).*sin(theta).*front_constant.*(c_2.*siAsiB+k_a*a.*A+k_b*((-b).*sin(aerfa).*(sqrt(1-A.^2)).*cos(B)+b.*cos(aerfa).*A))))./(1-exp(sin(theta).*front_constant.*(c_2.*siAsiB+k_a*a.*A+k_b*((-b).*sin(aerfa).*(sqrt(1-A.^2)).*cos(B)+b.*cos(aerfa).*A)))).*conj((1-exp((nz).*sin(theta).*front_constant.*(c_2.*siAsiB+k_a*a.*A+k_b*((-b).*sin(aerfa).*(sqrt(1-A.^2)).*cos(B)+b.*cos(aerfa).*A))))./(1-exp(sin(theta).*front_constant.*(c_2.*siAsiB+k_a*a.*A+k_b*((-b).*sin(aerfa).*(sqrt(1-A.^2)).*cos(B)+b.*cos(aerfa).*A)))));
sum_part_z(isnan(sum_part_z))=(nz)^2;
sum_part_z(isinf(sum_part_z))=(nz)^2;
end

