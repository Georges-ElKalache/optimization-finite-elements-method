function [S,K] = elimination(nomfile,f,Td)
[Nbpt,~,~,Refneu,~,~]=Lecmail(nomfile);
K=construction_matrice(nomfile);
%recourt a la methode de pseudo_elimination direct pour faciliter les
%calculs
D=(Refneu==3);

%ajuster la forme du second membre
m=sum(D);
td=sparse(Nbpt,1);
td(D,1)=Td*ones(m,1);
g=K*td;

S=f-g;
S(D,1)=Td*ones(m,1);

%ajuster la forme de la matrice K
K(:,D)=sparse(Nbpt,m);
K(D,:)=sparse(m,Nbpt);
K(D,D)=speye(m);

end




