function f= XAtomScattering( theta,lambda,atom,para)
%parameter    atom a1 b1 a2 b2 a3 b3 a4 b4 c
% global para;
s=sin(theta)/lambda;
%f=para(1)*exp(-para(2)*s.^2)+para(3)*exp(-para(4)*s.^2)+para(5)*exp(-para(6)*s.^2)+para(7)*exp(-para(8)*s.^2)+para(9);
f=bsxfun(@plus,bsxfun(@times,exp(-para(atom,3)*s.^2),para(atom,2))+bsxfun(@times,exp(-para(atom,5)*s.^2),para(atom,4))+bsxfun(@times,exp(-para(atom,7)*s.^2),para(atom,6))+bsxfun(@times,exp(-para(atom,9)*s.^2),para(atom,8)),para(atom,10));
%f = bsxfun(@plus,bsxfun(@times,exp(-1*bsxfun(@times,ones(31,1)*s.^2,reshape(para(atom,3),1,1,[]))),reshape(para(atom,2),1,1,[])) + bsxfun(@times,exp(-1*bsxfun(@times,ones(31,1)*s.^2,reshape(para(atom,5),1,1,[]))),reshape(para(atom,4),1,1,[]))  + bsxfun(@times,exp(-1*bsxfun(@times,ones(31,1)*s.^2,reshape(para(atom,7),1,1,[]))),reshape(para(atom,6),1,1,[])) + bsxfun(@times,exp(-1*bsxfun(@times,ones(31,1)*s.^2,reshape(para(atom,9),1,1,[]))),reshape(para(atom,8),1,1,[])),reshape(para(atom,10),1,1,[]));

end

