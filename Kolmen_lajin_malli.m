%T‰ss‰ koodissa tarkastellaan kolmeen lajin mallin ominaisarvoja,
%ratkaisuk‰yri‰ ja kolmiulotteisia R,C,P-k‰yri‰ ja tutkitaan
%tasapainopisteiden luonnetta.

clear all

%M‰‰ritell‰‰n laskuissa tarvittavat parametrit ja muuttujat. Tutkittavien
%parametrien edest‰ poistetaan prosenttimerkit koodin alusta ja 
%lopussa olevasta funktiosta.

syms C P R nonegative;

%Parametrit, kun tutkitaan asymptoottisesti stabiilia tasapainotilaa.

R_0=0.8; C_0=0.8; x_p=0.01; x_c=0.4; y_c=3.25; y_p=5;

%Parametrit, kun tutkitaan jaksollista ratkaisua.

%R_0=0.8; C_0=0.8; x_p=0.5; x_c=0.7; y_c=8; y_p=15;

%Parametrit kaoottiselle ratkaisulle.

%R_0=0.161; C_0=0.5; x_c=0.4; y_c=2.01; y_p=5;
%x_p = 0.08;
%x_p = 0.22;

%M‰‰ritell‰‰n differentaaliyht‰lˆ.

F(R,C,P)=[R*(1-R)-(x_c*y_c*C*R)/(R+R_0); 
    x_c*C*(-1+(y_c*R)/(R+R_0))-(x_p*y_p*C*P)/(C+C_0); 
    x_p*P*(-1+(y_p*C)/(C+C_0))];

%Ratkaistaan differentaaliyht‰lˆn tasapainopisteet ja luodaan niille omat
%vektorit.

p = solve(F == [0;0;0],R,C,P);
TP1=double([p.R(1) p.C(1) p.P(1)])
TP2=double([p.R(2) p.C(2) p.P(2)])
TP3=double([p.R(3) p.C(3) p.P(3)])
TP4=double([p.R(4) p.C(4) p.P(4)])
TP5=double([p.R(5) p.C(5) p.P(5)])
TP6=double([p.R(6) p.C(6) p.P(6)])


%Lasketaan Jaakobin matriisit eri tasapainopisteiss‰.

DF(R,C,P)=jacobian(F);
A1 = DF(p.R(1),p.C(1),p.P(1));
A2 = DF(p.R(2),p.C(2),p.P(2));
A3 = DF(p.R(3),p.C(3),p.P(3));
A4 = DF(p.R(4),p.C(4),p.P(4));
A5 = DF(p.R(5),p.C(5),p.P(5));
A6 = DF(p.R(6),p.C(6),p.P(6));

%Lasketaan matriisien ominaisarvot.

omarvo1 = double(eig(A1))
omarvo2 = double(eig(A2))
omarvo3 = double(eig(A3))
omarvo4 = double(eig(A4))
omarvo5 = double(eig(A5))
omarvo6 = double(eig(A6))

%Tutkitaan, ovatko ominaisarvojen reaaliosat negatiivisia. 
%T‰h‰n k‰ytet‰‰n Matlabin isAlways-komentoa, joka kertoo, 
%onko annettu ehto tosi vai ep‰tosi.

TP1ominaisarvo1=isAlways(real(omarvo1(1)) <= 0);
TP1ominaisarvo2=isAlways(real(omarvo1(2)) <= 0);
TP1ominaisarvo3=isAlways(real(omarvo1(3)) <= 0);

TP2ominaisarvo1=isAlways(real(omarvo2(1)) <= 0);
TP2ominaisarvo2=isAlways(real(omarvo2(2)) <= 0);
TP2ominaisarvo3=isAlways(real(omarvo2(3)) <= 0);

TP3ominaisarvo1=isAlways(real(omarvo3(1)) <= 0);
TP3ominaisarvo2=isAlways(real(omarvo3(2)) <= 0);
TP3ominaisarvo3=isAlways(real(omarvo3(3)) <= 0);

TP4ominaisarvo1=isAlways(real(omarvo4(1)) <= 0);
TP4ominaisarvo2=isAlways(real(omarvo4(2)) <= 0);
TP4ominaisarvo3=isAlways(real(omarvo4(3)) <= 0);

TP5ominaisarvo1=isAlways(real(omarvo5(1)) <= 0);
TP5ominaisarvo2=isAlways(real(omarvo5(2)) <= 0);
TP5ominaisarvo3=isAlways(real(omarvo5(3)) <= 0);

TP6ominaisarvo1=isAlways(real(omarvo6(1)) <= 0);
TP6ominaisarvo2=isAlways(real(omarvo6(2)) <= 0);
TP6ominaisarvo3=isAlways(real(omarvo6(3)) <= 0);

%Tulostetaan jokaisen tasapainopisteen ominaisarvojen totuusarvot. Jos
%kaikki totuusarvot ovat ykkˆsi‰, on kyseess‰ asymptoottisesti stabiili
%tasapainopiste.

Tasapainopiste1=[TP1ominaisarvo1 TP1ominaisarvo2 TP1ominaisarvo3]
Tasapainopiste2=[TP2ominaisarvo1 TP2ominaisarvo2 TP2ominaisarvo3]
Tasapainopiste3=[TP3ominaisarvo1 TP3ominaisarvo2 TP3ominaisarvo3]
Tasapainopiste4=[TP4ominaisarvo1 TP4ominaisarvo2 TP4ominaisarvo3]
Tasapainopiste5=[TP5ominaisarvo1 TP5ominaisarvo2 TP5ominaisarvo3]
Tasapainopiste6=[TP6ominaisarvo1 TP6ominaisarvo2 TP6ominaisarvo3]

%Ratkaistaan differentiaaliyht‰lˆ k‰ytt‰m‰ll‰ Matlabin ode45-komentoa, joka
%k‰ytt‰‰ lopussa olevaa funktiota.

x1 = [0.4;0.1; 0.8];
[t,x] = ode45(@funktio,[0 1000],x1);

%Piirret‰‰n R, C ja P kuvaajat ajan funktiona.
figure(1)
plot(t,x(:,1),t,x(:,2),t,x(:,3))
legend('R','C','P')

%Piirret‰‰n kolmiulotteinen R,C,P-k‰yr‰ ratkaisusta.
figure(2)
plot3(x(:,1),x(:,2),x(:,3))
xlabel('R')
ylabel('C')
zlabel('P')

%Differentiaaliyht‰lˆn ratkaisemiseen k‰ytetty funktio

function dy = funktio(t,x)

%Parametrit, kun tutkitaan asymptoottisesti stabiilia tasapainopistett‰.
%R_0=0.8; C_0=0.8; x_p=0.01; x_c=0.4; y_c=3.25; y_p=5;

%Parametrit, kun tutkitaan jaksollista ratkaisua.
%R_0=0.8; C_0=0.8; x_p=0.5; x_c=0.7; y_c=8; y_p=15;

%Parametrit kaoottiselle ratkaisulle.
%R_0=0.161; C_0=0.5; x_c=0.4; y_c=2.01; y_p=5;
%x_p = 0.08;
%x_p = 0.22;


dy = [x(1)*(1-x(1))-(x_c*y_c*x(2)*x(1))/(x(1)+R_0); 
    x_c*x(2)*(-1+(y_c*x(1))/(x(1)+R_0))-(x_p*y_p*x(2)*x(3))/(x(2)+C_0);
    x_p*x(3)*(-1+(y_p*x(2))/(x(2)+C_0))];
end
