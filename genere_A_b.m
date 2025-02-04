function [A,b,Temp]=genere_A_b(nomfile,n)
[~,Nbtri,Coorneu,~,Numtri,Reftri]=Lecmail(nomfile);
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
    
