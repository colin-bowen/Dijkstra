function [shortest_path_cost] = Dijkstra(principal,terminal,blockersx, blockersy, begindex, endex, polycost,resolution,radius)
%Implementation of Dijkstra's shortest path algorithm in MATLAB.
%Considering blocking polygons "blockers," with penalized weights to go
%through them, find the shortest path between two nodes. Begin by defining
%the spatial distribution of the graph, where the nodes are in opposite
%corners with the polygon in the middle. %principal is the location, in
%latitude then longitude, of the principal node. terminal is the location
%in LONGITUDE(X) then LATITUDE(Y) of the terminal node. blockers is an n-by-max
%input matrix where n is the number of blocking polygons and max is the max
%boundary input for any 1 polygon. resolution is a scalar that indicates in
%km the distance between mesh spacings in both the x and y directions. The
%mesh therefore will not necessarily be square which might skullfuck
%everything but we can see :)
%blockersx_find = blockersx(:);
%blockersy_find = blockersy(:);

%convert values from degrees to kilometers

 
%define the boundary of the mesh
far_left = min([principal(1); terminal(1); blockersx]);
far_right = max([principal(1); terminal(1); blockersx]);
bottom = min([principal(2); terminal(2); blockersy]);
top = max([principal(2); terminal(2); blockersy]);
resolution = km2deg(resolution, radius);
graph_x = far_left:resolution:far_right;%linspace(far_left-resolution,far_right+resolution,1/resolution*(abs(far_right-far_left+2*resolution)))'; %why use linspace when you can just space equally 
graph_y = bottom:resolution:top;%linspace(bottom-resolution,top+resolution,1/resolution*(abs(top-bottom+2*resolution)))'; 

[X,Y] = meshgrid(graph_x,graph_y); %creates a mesh of each of the x points and the y points
%[m,n] = size(X);
nodes = zeros(numel(X),2);
[q,~] = size(nodes);
for i = 1:q
    nodes(i,:) = [X(i) Y(i)]; %defining the node positions
end
%assign in the mesh the principal and terminal nodes
distance_from_mesh_principal = 100000;
distance_from_mesh_terminal = 100000;
for i = 1:numel(X) %mapping the actual node location onto the new node map
    dist_from_mesh_p = sqrt((X(i)-principal(1))^2 + (Y(i)-principal(2))^2);
    dist_from_mesh_t = sqrt((X(i)-terminal(1))^2 + (Y(i)-terminal(2))^2);
    if dist_from_mesh_p < distance_from_mesh_principal
        distance_from_mesh_principal = dist_from_mesh_p;
        principe_pio = [X(i) Y(i)];
        pp = i;
    end
    if dist_from_mesh_t < distance_from_mesh_terminal
        distance_from_mesh_terminal = dist_from_mesh_t;
        tincipe_tio = [X(i) Y(i)];
        tt = i;
    end
end

i = 1; %counter for node connection. 
j = 1; %counter for node number
s = zeros(length(nodes)*4,1); %each node has a maximum of four connections (node from)
t = zeros(length(nodes)*4,1); %node to. 
while 1
    if j > numel(X)
        break
    end
    if j <= numel(X)-length(graph_x) %if it's not the right hand side of the graph 
            t(i) = j+length(graph_x); %node directly to the side
            s(i) = j;
            if mod(j,length(graph_x)) ~= 0 %if it's not at the bottom of the graph
                t(i+1) = j+length(graph_x)+1; %node diagonally adjacent
                s(i+1) = j;
            end
            if mod(j,length(graph_x)) ~= 1 %if it's not at the top of the graph
                t(i+3) = j+length(graph_x)-1;
                s(i+3) = j;
            end
                
    %else
     %       t(i) = nan;
      %      s(i) = nan;
       %     t(i+1) = nan;
        %    s(i+1) = nan;
    end
    if mod(j,length(graph_x)) ~= 0 %if it's not the lower bound of the graph
            t(i+2) = j+1; %node directly below
            s(i+2) = j; 
    %else
     %       t(i+2) = nan;
      %      s(i+2) = nan;
    end
    i = i+4;
    j = j+1; 
end

s = sparse(s); %squeeze out the zeros where connections weren't made because they were on the edges of the graph
t = sparse(t);
connect_from = find(s); 
connect_to = find(t); %returns the indices in the connections, want to convert to node indices via modular arithmetic
%need to create weights of connections. evaluate based on geographical
%distance 
geo_dist = zeros(length(s),1);
for i = 1:length(connect_from) %for each of the connections, determine their Euclidean distance
    long1 = nodes(s(connect_from(i),1));
    long2 = nodes(t(connect_to(i),1));
    lat1 = nodes(s(connect_from(i),2));
    lat2 = nodes(t(connect_to(i),2));
    diff_long = long2 - long1;
    diff_lat = lat2 - lat1;
    a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
    c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
    d = c * radius; %step 3 of haversine function
    geo_dist(connect_from(i)) = d;
end
weights = geo_dist; %renaming it because i like using more memory i guess
%now consider if inpolygon 
[p,~] = size(blockersx);
inpoly = zeros(q,p);
%maybe can do a smaller mesh search because there polygon is more likely to
%be a problem in the smaller ones 
%nodes_to_check = nodes((ceil(q/5)):(ceil(4*q/5)), :); %check only the middle 3/5ish of nodes 
%[v,~] = size(nodes_to_check); actually can't do this because then you lose
%the node index
for i = 1:q %for the total number of nodes in the mesh to be checked
    for j = 1:p %for the blocking polys
        if inpolygon(nodes(i,1),nodes(i,2), blockersx(begindex(j):endex(j)), blockersy(begindex(j):endex(j)))
            inpoly(i,j) = 1;
        end
    end
end


inpoly = sparse(inpoly);
[R,S] = find(inpoly); %the indices in inpoly of the set R of nodes inside blockers S. 
%now need to penalize the connections that correspond to this node
connections_from = zeros(length(s),4);  
connections_to = zeros(length(t),4);
connections = zeros(length(s),8);
for i = 1:length(R) %for each of the nodes inside a polygon
    find_in_s = find(s==R(i))';
    find_in_t = find(t==R(i))';
    connections_from(i,1:length(find_in_s)) = find_in_s; %find in the principal node set the values of this node for outbound edges
    connections_to(i,1:length(find_in_t)) = find_in_t; %find in the terminal node set the values of this node for incident edges
    connections(i,:) = [connections_from(i,:), connections_to(i,:)];
    j = 1;
    while 1
        if j > length(connections(i,:))
            break
        end
        while connections(i,j) == 0
            j = j+1;
            if j > length(connections(i,:))
                break
            end
        end
        if j > length(connections(i,:))
                break
        end
        if j <=4 %if it's a connection from the node we're inspecting, then we want to check if the node TO is in the same poly
            if ismember(t(connections(i,j)),R) && inpolygon(nodes(t(connections(i,j)),1), nodes(t(connections(i,j)),2), blockersx(begindex(S(i)):endex(S(i))), blockersy(begindex(S(i)):endex(S(i))))     %geodist(connections(i,j)) < sqrt(deg2km(polyarea(blockersx(S(i),:),blockersy(S(i),:))))%if the adjoining node in the connection is also in the same polygon, determined by first finding it in the set of nodes inpolys and then finding whether same poly
                long1 = nodes(s(connections(i,j)),1);
                long2 = nodes(t(connections(i,j)),1);
                lat1 = nodes(s(connections(i,j)),2);
                lat2 = nodes(t(connections(i,j)),2);
                diff_long = long2 - long1;
                diff_lat = lat2 - lat1;
                a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                d = c * radius; %step 3 of haversine function
                dist_in_poly = d*(polycost(S(i))-1);
                
                
                %dist_in_poly = deg2km(sqrt((nodes(s(connections(i,j)),1)-nodes(t(connections(i,j)),1))^2 + ...
                %(nodes(s(connections(i,j)),2)-nodes(t(connections(i,j)),2))^2))*(polycost(S(i))-1);
            else %the node_to is outside the polygon and so there is an intersection point, which will correspond to the weight
                [x_poly_intersection, y_poly_intersection] = polyxpoly([nodes(s(connections(i,j)),1) nodes(t(connections(i,j)),1)], ...
                    [nodes(s(connections(i,j)),2), nodes(t(connections(i,j)),2)], blockersx(begindex(S(i)):endex(S(i))), blockersy(begindex(S(i)):endex(S(i))));
                
                long1 = nodes(s(connections(i,j)),1);
                long2 = x_poly_intersection(1);
                lat1 = nodes(s(connection(i,j)),1);
                lat2 = y_poly_intersection(1);
                diff_long = long2 - long1;
                diff_lat = lat2 - lat1;
                a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                d = c * radius; %step 3 of haversine function
                dist_in_poly = d*(polycost(S(i))-1);
                
                %dist_in_poly = deg2km(sqrt((nodes(s(connections(i,j)),1)-x_poly_intersection(1))^2 + ...
                 %   (nodes(s(connections(i,j)),2) - y_poly_intersection(1))^2))*(polycost(S(i))-1); %the distance in the poly of the node i and one of its (up to 6) j connections
                
            end %has a polycost corresponding to the node's home polygon
        else %if it's a connection TO the node we're inspecting, then we want to check if the node FROM is in the same poly
            if ismember(s(connections(i,j)),R) && inpolygon(nodes(s(connections(i,j)),1), nodes(s(connections(i,j)),2), blockersx(begindex(S(i)):endex(S(i))), blockersy(begindex(S(i)):endex(S(i)))) %checking if node FROM is also in poly 
                long1 = nodes(s(connections(i,j)),1);
                long2 = nodes(t(connections(i,j)),1);
                lat1 = nodes(s(connections(i,j)),2);
                lat2 = nodes(t(connections(i,j)),2);
                diff_long = long2 - long1;
                diff_lat = lat2 - lat1;
                a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                d = c * radius; %step 3 of haversine function
                dist_in_poly = d*(polycost(S(i))-1);
                %dist_in_poly = deg2km(sqrt((nodes(s(connections(i,j)),1)-nodes(t(connections(i,j)),1))^2 + ... if so, then it is the distance between them times the polycost
                %    (nodes(s(connections(i,j)),2)-nodes(t(connections(i,j)),2))^2))*(polycost(S(i))-1);
            else
                [x_poly_intersection, y_poly_intersection] = polyxpoly([nodes(s(connections(i,j)),1) nodes(t(connections(i,j)),1)], ...
                    [nodes(s(connections(i,j)),2), nodes(t(connections(i,j)),2)], blockersx(begindex(S(i)):endex(S(i))), blockersy(begindex(S(i)):endex(S(i))));
                
                long1 = nodes(s(connections(i,j)),1);
                long2 = x_poly_intersection(1);
                lat1 = nodes(s(connections(i,j)),2);
                lat2 = y_poly_intersection(1);
                diff_long = long2 - long1;
                diff_lat = lat2 - lat1;
                a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                d = c * radius; %step 3 of haversine function
                dist_in_poly = d*(polycost(S(i))-1);
                %dist_in_poly = deg2km(sqrt((nodes(s(connections(i,j)),1)-x_poly_intersection(1))^2 + ...
                %    (nodes(s(connections(i,j)),2) - y_poly_intersection(1))^2))*(polycost(S(i))-1); %the distance in the poly of the node i and one of its (up to 6) j connections
            end
        end
        if j<=4
            weights(connections_from(i,j)) = weights(connections_from(i,j)) + dist_in_poly; %add the penalized weight to the calculated weight. will include distance already, so we multiply by polycost-1 since the edge weight by distance for that portion was already calculated
        else
            weights(connections_to(i,j-4)) = weights(connections_to(i,j-4)) + dist_in_poly;
        end
        j = j+1;
    end
end
weights = sparse(weights);
s = s(connect_from);
t = t(connect_to);
weights = weights(connect_from);
s = full(s);
t = full(t);
weights = full(weights);
G = graph(s,t, weights); %make a directed graph with edge weights specified by the weight calculations from above
[shortest_path, shortest_path_cost] = shortestpath(G,pp,tt); %shortest path between the first node and the last node in the node set



end
