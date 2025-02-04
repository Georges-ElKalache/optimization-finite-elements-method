function T= prob_direct(x,y,nomfile,Td,A0)
[~,~,Coorneu,~,Numtri,~]=Lecmail(nomfile);
f=sec_membre(x,y,nomfile,A0);
[S,K]=elimination(nomfile,f,Td);
T=K\S;

figure;
T=full(T);
trisurf(Numtri,Coorneu(:,1),Coorneu(:,2),T);
xlabel('x');
ylabel('y');
zlabel('temperature');
title('Probleme Direct');
colorbar;

end
