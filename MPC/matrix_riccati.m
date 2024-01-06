function xdot = matrix_riccati(Ac,Bc,Q,Ro,x,t)

ro_inv=1/Ro;
R=[x(1),x(2), x(3), x(4), x(5) ,x(6), x(7), x(8), x(9);
   x(2),x(10),x(11),x(12),x(13),x(14),x(15),x(16),x(17);
   x(3),x(11),x(18),x(19),x(20),x(21),x(22),x(23),x(24);
   x(4),x(12),x(19),x(25),x(26),x(27),x(28),x(29),x(30);
   x(5),x(13),x(20),x(26),x(31),x(32),x(33),x(34),x(35);
   x(6),x(14),x(21),x(27),x(32),x(36),x(37),x(38),x(39); 
   x(7),x(15),x(22),x(28),x(33),x(37),x(40),x(41),x(42);
   x(8),x(16),x(23),x(29),x(34),x(38),x(41),x(43),x(44);
   x(9),x(17),x(24),x(30),x(35),x(39),x(42),x(44),x(45)];


riccati_eq = (-1)*( Ac.'*R  +  R*Ac - R*Bc* ro_inv* Bc.'*R + Q);

count=1;
    for l=1:1:size(Ac,1) %sutun
        for k=l:1:size(Ac,2)%satir
            xdot_pre(count,1)=riccati_eq(k,l);
            count=count+1;
        end
    end
    xdot=xdot_pre;
end
