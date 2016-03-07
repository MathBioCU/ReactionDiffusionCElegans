%%% Code for calculating eigenfuctions and values of the laplacian
%%% on a 2D rectangular domain and a 3D cylindrical domain
%%%
%%% Jacqui Wentz
%%% 3/7/2016

%% Set up parameter values (see pdeeig doc for what these values mean)
c = 1;
a = 0;
d = 1;

%% Calculate eigenfunctions/values on rectangular domain
numberOfPDE = 1;
model2D = createpde(numberOfPDE);

% Create domain.
% For rectangular solid: Row 1 contains 3, Row 2 contains # of edges
R = [3,4,0,10,10,0,1,1,0,0]'; 
g = decsg(R);
geometryFromEdges(model2D,g);

% Plot domain with edge labels
pdegplot(model2D, 'edgeLabels', 'on');
xlim([-1 11]);
axis equal
title 'Geometry With Edge Labels Displayed';
xlabel x
ylabel y

% Specify BCs:
% Apply Neumann boundary conditions to edges (i.e., derivative=0 at
% boundary)
applyBoundaryCondition(model2D,'Edge',[1,2,3,4],'g',0);

% Generate mesh over with to perform fem
generateMesh(model2D,'Hmax',0.2); % Hmax: Target maximum mesh edge length
figure
pdemesh(model2D);
xlim([-1 11]);
axis equal
xlabel x
ylabel y

% Solve for Laplacian eigenfunctions/values
[V,L] = pdeeig(model2D,c,a,d,[-inf,45]); %Solves: ???(c?u)+au=?du on ?

% Plot Results
for i=1:size(L,1)
    pdeplot(model2D,'xydata',V(:,i)) 
    title(['Eigenfunction for mode ' num2str(i) '. Eigenvalue: ' num2str(L(i))]);
    pause
end

%% Calculate eigenfunctions/values on cylindrical 3D domain
numberOfPDE = 1;
model3D = createpde(numberOfPDE);

% Manually create domain.
nxyPoints = 100; %Number of points to use along a circle in the xy plane
nzPoints = 100; %Number of points to use for the length of the cylinder
radius1 = 1; %Radius of the inner circle of the cylinder 
radius2 = 1.1; %Radius of the outter circule of the cylinder
length = 10; %Length of the cylinder
theta = linspace(0,2*pi,nxyPoints);
z = linspace(0,length,nzPoints);
x1 = radius1*sin(theta);
y1 = radius1*cos(theta);
x2 = radius2*sin(theta);
y2 = radius2*cos(theta);
vertices1 = [repmat(x1,1,nzPoints)', repmat(y1,1,nzPoints)', kron(z',ones(nxyPoints,1))];
vertices2 = [repmat(x2,1,nzPoints)', repmat(y2,1,nzPoints)', kron(z',ones(nxyPoints,1))];
vertices = [vertices1; vertices2]; %Matrix of vertcies for create shape

% Add manually created domain to model
shape = alphaShape(vertices); 
shape.Alpha = 1; % Change Alpha to make sure that there are no jagged edges, if this value is too large the ends of the cylinder will be filled in
plot(shape) % Verify shape
[elements,nodes] = boundaryFacets(shape);
geometryFromMesh(model3D,nodes',elements'); 

% Plot domain with face labels
hc = pdegplot(model3D, 'FaceLabels', 'on');
hc(1).FaceAlpha = 0.5; %Make faces transparent so you can see the labels
axis equal
title 'Geometry With Face Labels Displayed';
xlabel x
ylabel y

% Specify BCs:
% Apply Neumann boundary conditions to faces (i.e., derivative=0 at
% boundary)
applyBoundaryCondition(model3D,'Face',[1,2,3,4],'g',0);

% Generate mesh over which to solve fem
generateMesh(model3D,'Hmax',.5); % Hmax: Target maximum mesh edge length
figure
pdemesh(model3D);
ylim([-1.1 1.1]);
axis equal
xlabel x
ylabel y

% Solve for eigenfunctions/values
[V,L] = pdeeig(model3D,c,a,d,[-inf,10]);

% Plot Results
for i=1:size(L,1) 
    pdeplot3D(model3D,'colormapdata',V(:,i))
    title(['Eigenfunction for mode ' num2str(i) '. Eigenvalue: ' num2str(L(i))]);
    pause
end