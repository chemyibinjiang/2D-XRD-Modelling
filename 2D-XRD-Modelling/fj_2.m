function conju_sum = fj_2( theta,coef_2,point_num2,front_constant,A,B,a,b,x,y,z,aerfa,k_a,k_b,c_2)

conju_sum_1=0;
SIA=(sqrt(1-A.^2));
siAsiB=SIA.*sin(B);
for p=1:size(x,1)
    part1=A*x(p).*sin(theta).*front_constant.*a;
    part2=sin(theta).*(front_constant.*((-b).*sin(aerfa).*(sqrt(1-A.^2)).*cos(B)+b.*cos(aerfa).*A)*y(p));
    part3=sin(theta).*front_constant.*(c_2.*siAsiB+k_a*a.*A+k_b*((-b).*sin(aerfa).*(sqrt(1-A.^2)).*cos(B)+b.*cos(aerfa).*A))*z(p);
    conju_sum_1=conju_sum_1+coef_2(p,point_num2).*exp(part1+part2+part3);
end
conju_sum=conju_sum_1.*conj(conju_sum_1);
end

