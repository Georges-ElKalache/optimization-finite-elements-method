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

    
