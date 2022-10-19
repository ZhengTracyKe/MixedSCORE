load('TradeInService.mat');
Values = Values14 + Values15 + Values16 + Values17 + Values18;


%%%%%%%%%%%%%Summaries%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ExAll = Values(:,203);  %%%% reported total service to World
[~, sortind] = sort(ExAll, 'descend');
table(Economy(sortind), ExAll(sortind))


%%%%%%%%%%%%%Constructing the network%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Values = (Values + Values')/2;
%%%%%% Remove 'European Union', 'Extra EU Trade', 'World'
Inds = [1:62, 64:201];    
Names = Economy(Inds);
Values = Values(Inds,Inds);
Total = sum(Values,2);


[~, ind] = sort(Total, 'descend');
table(Names(ind), Total(ind))

Temp = (Total*Total').^(-0.5);
Weights = Values.*Temp;
A = double(Weights>0.02);


%%%%%%%%%%%%SCORE projection%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = graph(A, Names);
bins = conncomp(G)';
inds = find(bins==2);

A = A(inds,inds);
Countries = Names(inds);
[U, V] = eigs(A);
v=diag(V);
[~, ind] = sort(abs(v), 'descend');
v=v(ind);
U=U(:,ind);
R = U(:,2:3)./repmat(U(:,1),1,2);




%%%%%%%%%%%%Implementing Mixed-SCORE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Vertex hunting
n = size(R,1);
temp = find(R(:,2)==min(R(:,2)));  %remove the outlier 'Montenegro'
data = R([1:(temp-1), (temp+1):n],:);
[~, clusters] = kmeans(data, 25, 'replicates', 20); 
%hold on
%scatter(clusters(:,1), clusters(:,2), 'red', 'filled');
%%% We are supposed to run SVS to select 3 clusters, but for this example,
%%% vertex hunting is so obvious that we use this greedy approach below.
%%% The SVS algorithm gives the same output but takes longer time. 
vertices = zeros([3,2]);
vertices(1,:) = clusters(clusters(:,2)==max(clusters(:,2)),:);
vertices(2,:) = clusters(clusters(:,1)==min(clusters(:,1)),:);
vertices(3,:) = clusters(clusters(:,2)==min(clusters(:,2)),:);

%%%%% Membership estimation
tempPi = [ones([n,1]), R]/[ones([3,1]), vertices];
invb1 = (v(1) + diag(vertices*diag(v(2:3))*vertices')).^0.5;
tempPi = tempPi*diag(invb1);
tempPi = max(tempPi, 0);
Pi = tempPi./repmat(sum(tempPi,2),1,3);



%%%%% Plotting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%%%%% Plot the simplex first 
hold on
plot(vertices(1:2,1), vertices(1:2,2), 'k--');
plot(vertices([1,3],1), vertices([1,3],2), 'k--');
plot(vertices(2:3,1), vertices(2:3,2), 'k--');

%%%%% Plot the data cloud
scatter(R(:,1), R(:,2),'blue', 'filled')
axis([-7,6,-11,4]);  %%% change this if the eigenvector has a sign flip


%%%%% Select a few economies to show in the figure
hold on
select = [16, 102, 115, 66, 57, 55, 23, 93, 50, 20, 109, 110,...
    67, 39, 54, 97, 62, 104, 42, 105, 85, 88, 89, 56,...
    95, 107, 60, 111, 92, 14, 17, 63, 72, 1, 37];
Locations = R(select,:);
Text = Countries(select);
%%%%% Slightly adjust the positions of some text labels
Locations(select==57,1) = Locations(select==57,1)-0.2; %%% adjust 'Korea'
Locations(select==57,2) = Locations(select==57,2)+0.3;
Locations(select==23,2) = Locations(select==23,2)+0.1; %%% adjust 'China'
Locations(select==55,1) = Locations(select==55,1)-0.1; %%% adjust 'Japan'
Locations(select==55,2) = Locations(select==55,2)-0.2; 
Locations(select==20,1) = Locations(select==20,1)-0.3; %%% adjust'Canada'
Locations(select==20,2) = Locations(select==20,2)+0.25; 
Locations(select==109,2) = Locations(select==109,2)-0.1; %%% adjust'UK' 
Locations(select==54,2) = Locations(select==54,2)-0.1; %%% adjust 'Italy'
Locations(select==105,1) = Locations(select==105,1)-0.2; %%% adjust 'Turkey'
Locations(select==105,2) = Locations(select==105,2)+0.3; 
Locations(select==14,1) = Locations(select==14,1)-0.8; %%% adjust 'Bosnia'
Locations(select==14,2) = Locations(select==14,2)-0.15;
% Locations(select==84,1) = Locations(select==84,1)-0.5; %%%%% adjust 'Philippines'
% Locations(select==84,2) = Locations(select==84,2)+0.15;
Locations(select==56,1) = Locations(select==56,1)-0.3; %%%%% adjust 'Kazakhstan'
Locations(select==56,2) = Locations(select==56,2)-0.15;
%%%%% Shorten the names of some economies 
Text{select==57}='Korean';
Text{select==109}='UK';
Text{select==110}='USA';
Text{select==14}='Bosnia-Herzegovina';
Text{select==111}='Kosovo';
Text{select==89}='Russian';
scatter(R(select,1), R(select,2), 'green', 'filled')
%%%%% Add text labels showing economy names
text(Locations(:,1), Locations(:,2), Text, 'fontsize', 8, 'color', 'red');
hold off






