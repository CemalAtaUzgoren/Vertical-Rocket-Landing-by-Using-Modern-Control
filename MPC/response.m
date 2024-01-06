function zdot = response(K_T,Ac,Bc,z,t)

Z=[z(1);z(2);z(3);z(4);z(5);z(6);z(7);z(8);z(9)];

u1=-K_T(1,1)*Z(1,1)-K_T(1,2)*Z(2,1)-K_T(1,3)*Z(3,1)-K_T(1,4)*Z(4,1)-K_T(1,5)*Z(5,1)-K_T(1,6)*Z(6,1)-K_T(1,7)*Z(7,1)-K_T(1,8)*Z(8,1)-K_T(1,9)*Z(9,1);
u2=-K_T(2,1)*Z(1,1)-K_T(2,2)*Z(2,1)-K_T(2,3)*Z(3,1)-K_T(2,4)*Z(4,1)-K_T(2,5)*Z(5,1)-K_T(2,6)*Z(6,1)-K_T(2,7)*Z(7,1)-K_T(2,8)*Z(8,1)-K_T(2,9)*Z(9,1);

U=[u1;u2];

system_eq = Ac*Z + Bc*U;
zdot=zeros(9,1);

for i=1:1:size(system_eq,1)
    zdot(i,1)=system_eq(i,1);
end

end