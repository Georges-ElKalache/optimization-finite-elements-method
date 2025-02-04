function R=genere_R_alea(n)
A0=1e4;
%generer une matrice aleatoire de taille n contenant dans sa 2eme et 3eme
%colonnes les coordonnees des resistances
R=rand(n,3); 
R(:,1)=A0*ones(n,1);
R(:,2)=((2*R(:,2))-1); %afin d'avoir des valeures entre -1 t 1
R(:,3)=((2*R(:,3))-1); %afin d'avoir des valeures entre -1 t 1

%voir si une resistance se trouve dans le domaine du par-choc pour
%regenerer une listes de positions de resistances valables
for i=1:n
    while (((-0.5<=R(i,2))&&(0.5>=R(i,2)))&&((-0.2<=R(i,3))&&(0.2>=R(i,3))))
        R(:,2:3)=rand(n,2);
        R(:,2)=(2*R(:,2)-1);
        R(:,3)=(2*R(:,3)-1);
    end
end

end