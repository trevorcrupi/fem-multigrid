%Input:
%mesh_level
%See meshgeometry.m to change your doamin.


mesh_level=2;

model=createpde(1);
[a, b, c]=meshgeometry(1);
g=decsg(a,b,c);
geometryFromEdges(model,g);

[p,e,t] = initmesh(g,'hmax',inf);
 
for ii=1:mesh_level-1

  [p,e,t] = refinemesh(g,p,e,t,'regular');
    %pdemesh(p,e,t)
    %figure
    %pdeplot(p,e,t,'ElementLabels','on')
end


[~,N_node]=size(p); %N_node: number of nodes
node=p';  %node: the N_node by 2 matrix that has the x,y-coordinates of each node.
[~,N_ele]=size(t); %N_ele is the number of triangles.
ele=t(1:3,1:N_ele);
ele=ele'; %ele: the N_ele by 3 matrix that has the three vertex numbers for each triangle.
TR=triangulation(ele,node);
edge=edges(TR); %edge: the N_edge by 2 matrix that has the two vertex numbers for making each edge.
[N_edge,E2]=size(edge);  %N_edge is the number of edges.


%new_ele is a N_ele by 3 matrix that saves the edge numbers of each
%triangle. This will be saved in a "new_ele(mesh_level).mat") after you run
%this program.



new_ele=zeros(N_ele,3); 


for k=1:N_ele
  
    node_vector1=[ele(k,1),ele(k,2)];
    node_vector2=[ele(k,2),ele(k,3)];
    node_vector3=[ele(k,3),ele(k,1)];
    
    count=0;
    
    for k_edge=1:N_edge
        if isempty(setdiff(node_vector1,edge(k_edge,:)))==1
            new_ele(k,1)=k_edge;
            count=count+1;
        elseif isempty(setdiff(node_vector2,edge(k_edge,:)))==1
            new_ele(k,2)=k_edge;
            count=count+1;
        elseif isempty(setdiff(node_vector3,edge(k_edge,:)))==1
            new_ele(k,3)=k_edge;
            count=count+1;
        end
        
        if count==3
           break; 
        end
    end

end

save(['new_ele', num2str(mesh_level),'.mat'],'new_ele');


