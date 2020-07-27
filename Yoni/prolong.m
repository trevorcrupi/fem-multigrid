%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Minah Oh
%Date: 07/02/2020
%
%This program generates the prolongation matrix between two Fourier finite
%element spaces Z_{H,Fk} and Z_{h,Fk} that were first described in the paper 
%by Lacoste titled Solution of Maxwell equation in axisymmetric geometry by Fourier
%series decompostion and by use of H (rot) conforming finite element,
%Numerische Mathematik volume 84, pages577–609(2000).
%
%Z_{H,Fk} is the Fourier-FEM space corresponding to mesh_level L.
%Z_{h,Fk} is the Fourier-FEM space corresponding to mesh_level L+1.
%Mesh level L+1 is obtained by connecting the midpoints of all edges in
%mesh level L, so Z_{H,Fk} is a subspace of Z_{h,Fk}. 
%
%Input: 
%L: mesh_level
%meshgeometry.m: saves the "gd, sf, ns" variables that represents the
%domain of problem.
%
%Output: Prolongation matrix P that extends vectors in Z_{H,Fk} to
%Z_{h,Fk}, i.e., if X is a vector representing a function in Z_{H,Fk}, 
%Y=P*X is a longer vector that represents the same function in Z_{h,Fk}.
%The prolongation matrix P is independent of Fk.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function P=prolong(L) 

% % model=createpde(1);
% % [a,b,c]=meshgeometry(1);
% % g=decsg(a,b,c);
% % geometryFromEdges(model,g);
% % [p,e,t]=initmesh(g,'hmax',inf);
% % %pdemesh(p,e,t)
% % 
% % for ii=1:L-1
% %    [p,e,t]=refinemesh(g,p,e,t,'regular'); 
% % end
% Replace with:
mystr             = ['PETForYoni/PETForYoni' num2str(L) '.mat'];
load(mystr);

[~,N_node]=size(p);
node=p';
[~,N_ele]=size(t);
ele=t(1:3,1:N_ele);
ele=ele';
TR=triangulation(ele,node);
edge=edges(TR);
[N_edge,~]=size(edge);
load(['newEle/new_ele',num2str(L),'.mat']);

% % [p2,~,t2]=refinemesh(g,p,e,t,'regular');
% Replace with:
temp1 = load(['PETForYoni/PETForYoni' num2str(L+1) '.mat']);
temp2 = struct2cell(temp1);
p2 = temp2{2};
t2 = temp2{3};

[~,N_node2]=size(p2);
node2=p2';
[~,N_ele2]=size(t2);
ele2=t2(1:3,1:N_ele2);
ele2=ele2';
TR2=triangulation(ele2,node2);
edge2=edges(TR2);
[N_edge2,~]=size(edge2);
temp=load(['newEle/new_ele',num2str(L+1),'.mat']);
new_ele2=cell2mat(struct2cell(temp));

NedelecBasis=zeros(3,3,N_ele);
P1Basis=zeros(3,3,N_ele);

A1=sparse(N_edge2,N_edge);
A2=sparse(N_node2,N_node);

for k2=1:N_ele2
       
        k=mod(k2,N_ele); %big trianlge number that include's the kk2-th small triangle. 
       if k==0
          k=N_ele; 
       end
       

    P1Basis(:,:,k)=[node(ele(k,1),1), node(ele(k,1),2), 1; ...
                    node(ele(k,2),1), node(ele(k,2),2), 1;...
                    node(ele(k,3),1), node(ele(k,3),2), 1]\[1 0 0; 0 1 0; 0 0 1];

    %Each column saves the NedelecBasis function formula for each edge.
    NedelecBasis(:,:,k)=[node(edge(new_ele(k,1),1),1)*node(edge(new_ele(k,1),2),2)-node(edge(new_ele(k,1),2),1)*node(edge(new_ele(k,1),1),2), ...
           node(edge(new_ele(k,1),2),1)-node(edge(new_ele(k,1),1),1),...
           node(edge(new_ele(k,1),2),2)-node(edge(new_ele(k,1),1),2);...
           node(edge(new_ele(k,2),1),1)*node(edge(new_ele(k,2),2),2)-node(edge(new_ele(k,2),2),1)*node(edge(new_ele(k,2),1),2), ...
           node(edge(new_ele(k,2),2),1)-node(edge(new_ele(k,2),1),1),...
           node(edge(new_ele(k,2),2),2)-node(edge(new_ele(k,2),1),2);...
           node(edge(new_ele(k,3),1),1)*node(edge(new_ele(k,3),2),2)-node(edge(new_ele(k,3),2),1)*node(edge(new_ele(k,3),1),2),...
           node(edge(new_ele(k,3),2),1)-node(edge(new_ele(k,3),1),1),...
           node(edge(new_ele(k,3),2),2)-node(edge(new_ele(k,3),1),2);...
          ]\[1 0 0; 0 1 0; 0 0 1];
    
      for s2=1:3
      for s=1:3
   A1(new_ele2(k2,s2),new_ele(k,s))=...
       (node2(edge2(new_ele2(k2,s2),1),1)*node2(edge2(new_ele2(k2,s2),2),2)-node2(edge2(new_ele2(k2,s2),2),1)*node2(edge2(new_ele2(k2,s2),1),2)).*NedelecBasis(1,s,k)+(node2(edge2(new_ele2(k2,s2),2),1)-node2(edge2(new_ele2(k2,s2),1),1)).*NedelecBasis(2,s,k)+(node2(edge2(new_ele2(k2,s2),2),2)-node2(edge2(new_ele2(k2,s2),1),2)).*NedelecBasis(3,s,k);
 
   A2(ele2(k2,s2),ele(k,s))=... 
        P1Basis(1,s,k).*node2(ele2(k2,s2),1) + P1Basis(2,s,k).*node2(ele2(k2,s2),2) + P1Basis(3,s,k);
      end
      end
       
end

P=blkdiag(A1,A2);

%size(P)


%Checking if prolongation matrix is working correctly. Done.
%{
Fk=1;
w=rand(N_edge+N_node,1);
c=P*w;

error=0;

for k=1:N_ele2
    
    [X,Y,Wx,Wy]=triquad(8,[node2(ele2(k,1),1),node2(ele2(k,1),2); node2(ele2(k,2),1),node2(ele2(k,2),2); node2(ele2(k,3),1),node2(ele2(k,3),2)]);
       

    P1Basis2(:,:,k)=[node2(ele2(k,1),1), node2(ele2(k,1),2), 1; ...
                     node2(ele2(k,2),1), node2(ele2(k,2),2), 1;...
                     node2(ele2(k,3),1), node2(ele2(k,3),2), 1]\[1 0 0; 0 1 0; 0 0 1];

    %Each column saves the Nedele2cBasis function formula for each edge.
    NedelecBasis2(:,:,k)=...
          [node2(edge2(new_ele2(k,1),1),1)*node2(edge2(new_ele2(k,1),2),2)-node2(edge2(new_ele2(k,1),2),1)*node2(edge2(new_ele2(k,1),1),2), ...
           node2(edge2(new_ele2(k,1),2),1)-node2(edge2(new_ele2(k,1),1),1),...
           node2(edge2(new_ele2(k,1),2),2)-node2(edge2(new_ele2(k,1),1),2);...
           node2(edge2(new_ele2(k,2),1),1)*node2(edge2(new_ele2(k,2),2),2)-node2(edge2(new_ele2(k,2),2),1)*node2(edge2(new_ele2(k,2),1),2), ...
           node2(edge2(new_ele2(k,2),2),1)-node2(edge2(new_ele2(k,2),1),1),...
           node2(edge2(new_ele2(k,2),2),2)-node2(edge2(new_ele2(k,2),1),2);...
           node2(edge2(new_ele2(k,3),1),1)*node2(edge2(new_ele2(k,3),2),2)-node2(edge2(new_ele2(k,3),2),1)*node2(edge2(new_ele2(k,3),1),2),...
           node2(edge2(new_ele2(k,3),2),1)-node2(edge2(new_ele2(k,3),1),1),...
           node2(edge2(new_ele2(k,3),2),2)-node2(edge2(new_ele2(k,3),1),2);...
          ]\[1 0 0; 0 1 0; 0 0 1];
      
      EBx2=@(x,y,t) NedelecBasis2(2,t,k)./Fk.*x-NedelecBasis2(1,t,k)./Fk.*x.*y;
      EBy2=@(x,y,t) NedelecBasis2(3,t,k)./Fk.*x+NedelecBasis2(1,t,k)./Fk.*x.^2;
     
      
      VBx2=@(x,y,t) -1./Fk.*P1Basis2(3,t,k)-P1Basis2(1,t,k)./Fk.*x-1./Fk.*P1Basis2(2,t,k).*y;
      VBth2=@(x,y,t) P1Basis2(3,t,k)+P1Basis2(1,t,k).*x+P1Basis2(2,t,k).*y;

      
   uhx2=@(x,y) c(new_ele2(k,1)).*EBx2(x,y,1)+c(new_ele2(k,2)).*EBx2(x,y,2)+c(new_ele2(k,3)).*EBx2(x,y,3)...
               +c(N_edge2+ele2(k,1)).*VBx2(x,y,1)+c(N_edge2+ele2(k,2)).*VBx2(x,y,2)+c(N_edge2+ele2(k,3)).*VBx2(x,y,3);
   uhth2=@(x,y) c(N_edge2+ele2(k,1)).*VBth2(x,y,1)+c(N_edge2+ele2(k,2)).*VBth2(x,y,2)+c(N_edge2+ele2(k,3)).*VBth2(x,y,3); 
   uhy2=@(x,y) c(new_ele2(k,1)).*EBy2(x,y,1)+c(new_ele2(k,2)).*EBy2(x,y,2)+c(new_ele2(k,3)).*EBy2(x,y,3);
          
       k0=mod(k,N_ele); %big trianlge number that include's the kk2-th small triangle. 
       if k0==0
          k0=N_ele; 
       end
       
      EBx=@(x,y,t) NedelecBasis(2,t,k0)./Fk.*x-NedelecBasis(1,t,k0)./Fk.*x.*y;
      EBy=@(x,y,t) NedelecBasis(3,t,k0)./Fk.*x+NedelecBasis(1,t,k0)./Fk.*x.^2;
     
      
      VBx=@(x,y,t) -1./Fk.*P1Basis(3,t,k0)-P1Basis(1,t,k0)./Fk.*x-1./Fk.*P1Basis(2,t,k0).*y;
      VBth=@(x,y,t) P1Basis(3,t,k0)+P1Basis(1,t,k0).*x+P1Basis(2,t,k0).*y;
      
      
   uhx=@(x,y) w(new_ele(k0,1)).*EBx(x,y,1)+w(new_ele(k0,2)).*EBx(x,y,2)+w(new_ele(k0,3)).*EBx(x,y,3)...
             +w(N_edge+ele(k0,1)).*VBx(x,y,1)+w(N_edge+ele(k0,2)).*VBx(x,y,2)+w(N_edge+ele(k0,3)).*VBx(x,y,3);
   uhth=@(x,y) w(N_edge+ele(k0,1)).*VBth(x,y,1)+w(N_edge+ele(k0,2)).*VBth(x,y,2)+w(N_edge+ele(k0,3)).*VBth(x,y,3); 
   uhy=@(x,y) w(new_ele(k0,1)).*EBy(x,y,1)+w(new_ele(k0,2)).*EBy(x,y,2)+w(new_ele(k0,3)).*EBy(x,y,3);
       
   
   f=@(x,y) ((uhx(x,y)-uhx2(x,y)).^2+(uhth(x,y)-uhth2(x,y)).^2+(uhy(x,y)-uhy2(x,y)).^2).*x;  
   out=Wx'*feval(f,X,Y)*Wy;
   error=error+out;

end

   error
%}

end
