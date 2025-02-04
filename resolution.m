function [alpha,T_sol]=resolution(nomfile,n)
[~,~,Coorneu,~,Numtri,~]=Lecmail(nomfile);
[A,b,Temp]=genere_A_b(nomfile,n);
alpha=A\b;
%clacul de la solution T_sol
neg=(alpha<0);
alpha(neg)=0; %pas de resistance negative
T_sol=prob_direct(0.75,0.75,nomfile,300,0);
for i=1:n
    T_sol=T_sol+alpha(i)*Temp(:,i);
end

T_sol=full(T_sol);
figure
trisurf(Numtri,Coorneu(:,1),Coorneu(:,2),T_sol(:,1));
xlabel('x');
ylabel('y');
zlabel('temperature');
title('Probleme Inverse');
colorbar;

end