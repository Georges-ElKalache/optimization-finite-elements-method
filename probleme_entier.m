% Projet MO102 Roy AlSaddi - Georges El Kalache
function [alpha,T,T_sol]=probleme_entier(nomfile)
%%%%%%%%%%%%%%%%%%%%%%-Probleme Direct-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fonction de lecture de maillage
function [Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri]=Lecmail(nomfile) 
fid = fopen(sprintf('%s.amdba',nomfile),'r') ;
N = fscanf(fid,'%i') ;
Nbpt = N(1) ;
Nbtri = N(2) ;
line = fgets(fid) ;
tmp = fscanf(fid,'%f',[4,Nbpt]) ;
Coorneu = tmp(2:3,:)' ;
Refneu = tmp(4,:)' ;
tmp = fscanf(fid,'%i',[5,Nbtri]) ;
Numtri=tmp(2:4,:)' ;
Reftri=tmp(5,:)' ;
end

%Fonction pour la construction de la matrice de rigidite K
function [K] = construction_matrice(nomfile)
[Nbpt,Nbtri,Coorneu,~,Numtri,Reftri]=Lecmail(nomfile);
ka=1;kp=10; %valeurs donnee par hypothese pour l'air et le par-choc
K=sparse(Nbpt,Nbpt);
for l=1:Nbtri
    kloc=zeros(3,3);
    %coordonnees des trois sommets du triangle l
    x1=Coorneu(Numtri(l,1),1);y1=Coorneu(Numtri(l,1),2);
    x2=Coorneu(Numtri(l,2),1);y2=Coorneu(Numtri(l,2),2);
    x3=Coorneu(Numtri(l,3),1);y3=Coorneu(Numtri(l,3),2);
    delta=(x2-x3)*(y3-y1)-(x3-x1)*(y2-y3);
    aire=abs(delta)/2;
    %voir ou se trouve triangle
    if Reftri(l)==1
        k=ka;
    end
    if Reftri(l)==2
        k=kp;
    end
    %calcul de la matrice Kloc intermediaire avec les coordonnees
    %barycentriques
    kloc(1,1)=k*(aire/(delta*delta))*((y2-y3)*(y2-y3)+(x2-x3)*(x2-x3));
    kloc(1,2)=k*(aire/(delta*delta))*((y2-y3)*(y3-y1)+(x2-x3)*(x3-x1)); 
    kloc(1,3)=k*(aire/(delta*delta))*((y2-y3)*(y1-y2)+(x2-x3)*(x1-x2));  
    kloc(2,2)=k*(aire/(delta*delta))*((y3-y1)*(y3-y1)+(x3-x1)*(x3-x1)); 
    kloc(2,3)=k*(aire/(delta*delta))*((y3-y1)*(y1-y2)+(x3-x1)*(x1-x2));
    kloc(3,3)=k*(aire/(delta*delta))*((y1-y2)*(y1-y2)+(x1-x2)*(x1-x2));
    kloc(2,1)=kloc(1,2);kloc(3,1)=kloc(1,3);kloc(3,2)=kloc(2,3);
    %construction de la matrice K par l'intermediaire de Kloc
    for k=1:3
        for m=1:3
            i=Numtri(l,k);
            j=Numtri(l,m);
            K(i,j)=K(i,j)+kloc(k,m);
        end
    end
end
end

%Fonction pour construire le second membre f
function f = sec_membre(x,y,~,A0)
f=sparse(Nbpt,1);
for l=1:Nbtri
    %Coordonnees des trois sommets du triangle l
    x1=Coorneu(Numtri(l,1),1);y1=Coorneu(Numtri(l,1),2);
    x2=Coorneu(Numtri(l,2),1);y2=Coorneu(Numtri(l,2),2);
    x3=Coorneu(Numtri(l,3),1);y3=Coorneu(Numtri(l,3),2);
    delta=(x2-x3)*(y3-y1)-(x3-x1)*(y2-y3);
    %calcul des trois coordonnees barycentriques
    lam1=(1/delta)*((y2-y3)*(x-x3)-(x2-x3)*(y-y3));
    lam2=(1/delta)*((y3-y1)*(x-x1)-(x3-x1)*(y-y1));
    lam3=(1/delta)*((y1-y2)*(x-x2)-(x1-x2)*(y-y2));
    %voir ou se trouve la resistance pour utiliser la bonne formule
    if (0<=lam1) && (0<=lam2) && (0<=lam3) && (1>=lam1) && (1>=lam2) && (1>=lam3)
        for k=1:3
            if k==1
                f(Numtri(l,k))=A0*lam1;
            end
            if k==2
                f(Numtri(l,k))=A0*lam2;
            end
            if k==3
                f(Numtri(l,k))=A0*lam3;
            end
        end
    end
end
end


%Recourt a la methode de pseudo_elimination direct pour faciliter les
%calculs et obtenir des matrices creuses et symetriques
function [S,K] = elimination(nomfile,f,Td)
K=construction_matrice(nomfile);
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

%Fonction pour resoudre le probleme direct
function T= prob_direct(x,y,nomfile,Td,A0)
f=sec_membre(x,y,nomfile,A0);
[S,K]=elimination(nomfile,f,Td);
T=K\S;
end

%probleme direct:
A0=1e4;
Td=300;
x=0.75;y=0.75;
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri]=Lecmail(nomfile);
K=construction_matrice(nomfile);
f=sec_membre(x,y,nomfile,A0);
[S,K]=elimination(nomfile,f,Td);
T=prob_direct(x,y,nomfile,Td,A0);
%code pour visualiser la carte de temperature
figure;
T=full(T);
trisurf(Numtri,Coorneu(:,1),Coorneu(:,2),T);
xlabel('x');
ylabel('y');
zlabel('temperature');
title('Probleme Direct');
colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-Probleme Inverse-%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Fonction pour generer une matrice aleatoire de taille n contenant dans sa
% 2eme et 3eme colonnes les coordonnees des resistances
function R=genere_R_alea(n) 
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

%Fonction pour generer Temp contenant les vecteurs de temperatures en
%colonnes solution du probleme direct pour chaque resistance
function Temp=genere_temp(nomfile,n)
R=genere_R_alea(n);
Temp=sparse(Nbpt,n);
for i=1:n
    Temp(:,i)=prob_direct(R(i,2),R(i,3),nomfile,0,1e4);
end
end


%Fonction pour construire la matrice A ainsi que le vecteur b du probleme
%inverse
function [A,b,Temp]=genere_A_b(nomfile,n)
Topt=500;beta=1e-4; %valeures donnees par hypothese
%construction de la matrice A
A=beta*eye(n,n);
%generer le matrice Temp contenant les valeures t(j,r)
Temp=genere_temp(nomfile,n);
for l=1:Nbtri
    if Reftri(l)==1
        x1=Coorneu(Numtri(l,1),1);y1=Coorneu(Numtri(l,1),2);
        x2=Coorneu(Numtri(l,2),1);y2=Coorneu(Numtri(l,2),2);
        x3=Coorneu(Numtri(l,3),1);y3=Coorneu(Numtri(l,3),2);
        delta=(x2-x3)*(y3-y1)-(x3-x1)*(y2-y3);
        aire=abs(delta)/2;
        %utilisation des coordonnees barycentriques pour la construction de
        %M(3*3) intermediaire
        M=sparse(3,3);
        for i=1:3
            for j=1:3
                if i==j
                    M(i,j)=aire/6;
                else 
                    M(i,j)=aire/12;
                end
            end 
        end
        %construction de la matrice A
        for r=1:n
            for s=1:n
                for i=1:3
                    for j=1:3
                        A(r,s)=A(r,s) + M(i,j)*Temp(Numtri(l,i),r)*Temp(Numtri(l,j),s);
                    end
                end
            end
        end
    end
end
%generer la matrice T0 
T0=prob_direct(0.75,0.75,nomfile,300,0);
%calcul du second membre le vecteur B
b=sparse(n,1);
for l=1:Nbtri
      x1=Coorneu(Numtri(l,1),1);y1=Coorneu(Numtri(l,1),2);
      x2=Coorneu(Numtri(l,2),1);y2=Coorneu(Numtri(l,2),2);
      x3=Coorneu(Numtri(l,3),1);y3=Coorneu(Numtri(l,3),2);
      delta=(x2-x3)*(y3-y1)-(x3-x1)*(y2-y3);
      if Reftri(l)==1
          for r=1:n
              for i=1:3
                  b(r,1)=b(r,1)+(delta/6)*Temp(Numtri(l,i),r)*(Topt-T0(Numtri(l,i)));
              end
          end
      end
end
end

%Fonction pour resoudre le probleme imverse et visualiser la carte de
%temperature
function [alpha,T_sol]=resolution(nomfile,n)
[A,b,Temp]=genere_A_b(nomfile,n);
alpha=A\b;
%calcul de la solution T_sol
neg=(alpha<0);
alpha(neg)=0; %pas de resistance negative
T_sol=prob_direct(0.75,0.75,nomfile,300,0);
for i=1:n
    T_sol=T_sol+alpha(i)*Temp(:,i);
end
%visualisaton de la carte de temperature
T_sol=full(T_sol);
figure
trisurf(Numtri,Coorneu(:,1),Coorneu(:,2),T_sol(:,1));
xlabel('x');
ylabel('y');
zlabel('temperature');
title('Probleme Inverse');
colorbar;
end

%Probleme Inverse:
n=6;
R=genere_R_alea(n);
Temp=genere_temp(nomfile,n);
[A,b]=genere_A_b(nomfile,n);
[alpha,T_sol]=resolution(nomfile,n);

end