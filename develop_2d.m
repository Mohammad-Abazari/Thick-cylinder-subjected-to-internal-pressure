%% Preamble
clc; clear; 
% startup
set(groot,'DefaultTextInterpreter','latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(findall(gcf,'-property','FontSize'),'FontSize',14);
format compact;
% % set(groot, 'DefaultFigureRenderer', 'painters');
close all;

% % Problem statement
% Geometry description ===================================================
Ri = 40; % Inner radius Ri, mm
Ro = 80; % Outer radius Ro, mm
L = 1000; % Length L, mm
Pi = 17; % Inside pressure (radial pressure), N/mm2
Po = 0; % Outside pressure (radial pressure), N/mm2
meshsize_max = 8; % maximum mesh dimension, mm
meshsize_min = 1; % minimum mesh dimension, mm
mesh_order = 'linear'; % or quadratic

% Material Properties ====================================================
E = 210e3; % Modulus of elasticity E, N/mm2
nu = 0.3; % Poisson's ratio \nu
% rho = 8000; % Mass density \rho, kg/m3 irrelevent for static model

% Measurement
theta = pi/4;

%% Model creation =========================================================
model = createpde('structural','static-planestress');
% Introduce geometry =====================================================
importGeometry(model,'mesh2d.stl');

f1 = figure('Position',[100,40,1200,1000],'Renderer','painters');
subplot(2,2,1);
pdegplot(model,'EdgeLabels','on','FaceLabels','on','VertexLabels','on');
title('Cylinder with Face/Edge/Vertex Labels')

% Assign Material Values
structuralProperties(model,'YoungsModulus',E, ...
                           'PoissonsRatio',nu);
% Applying boundary conditions
structuralBC(model,'Edge',1,'Constraint','symmetric');
structuralBC(model,'Edge',4,'Constraint','symmetric');
% 
% % Apply pressure loads in and out
structuralBoundaryLoad(model,"Edge",3,"Pressure",Pi);
structuralBoundaryLoad(model,"Edge",2,"Pressure",Po);
% 
% % Mesh description, Generate mesh
mesh = generateMesh(model,'Hmax',meshsize_max,'Hmin',meshsize_min,'GeometricOrder',mesh_order);
subplot(2,2,2);
pdeplot(model);
title(['Mesh with Triangular Elements with $\Delta_m=',num2str(meshsize_max),'$mm']);

% Specify line of nodes considered for analysis
R = ((Ri):(Ro)); n_ids = zeros(size(R));
for i = 1:length(R)
    r = R(i);
    point = r*[cos(theta);sin(theta)];
    n_ids(i) = findNodes(mesh,'nearest',point);
end
subplot(2,2,3);
pdemesh(model,'NodeLabels','on'); hold on;
plot(mesh.Nodes(1,n_ids),mesh.Nodes(2,n_ids),'or','MarkerFaceColor','g')
axis off;
subplot(2,2,4);
pdemesh(model); hold on
plot(mesh.Nodes(1,n_ids),mesh.Nodes(2,n_ids),'or','MarkerFaceColor','g')
axis off;

title(['$\theta=',num2str(rad2deg(theta)),'^\circ$'])
legend('Mesh elements','Edges','Interpolation nodes');

% print('-f',['Figures/fig01m',num2str(meshsize_max)],'-dsvg')
%% Calculate solution
result = solve(model);
%% Visualize
figure('Position',[100,100,1900,600],'Renderer','painters');
subplot(1,4,1);
pdeplot(model,'XYData',result.Stress.sxx,'ColorMap','jet')
title('Stress $\sigma_{xx}$');
axis equal;
subplot(1,4,2);
pdeplot(model,'XYData',result.Stress.syy,'ColorMap','jet')
title('Stress $\sigma_{yy}$');
axis equal;
subplot(1,4,3);
pdeplot(model,'XYData',result.Stress.sxy,'ColorMap','jet')
title('Stress $\sigma_{xy}=\sigma_{yx}$');
axis equal;
subplot(1,4,4);
pdeplot(model,'XYData',result.VonMisesStress,'ColorMap','jet');
title('Von-Mises stress');
axis equal;
print('-f',['Figures/fig02m',num2str(meshsize_max)],'-dsvg')
% save(['Results2dM',num2str(meshsize_max),'m',num2str(meshsize_min),mesh_order,'.mat'],'result','model')
%% More
st = @(r) (Pi*Ri^2-Po*Ro^2)/(Ro^2-Ri^2) + ((Ri^2*Ro^2)./r.^2)*(Pi-Po)/((Ro^2-Ri^2));
sr = @(r) (Pi*Ri^2-Po*Ro^2)/(Ro^2-Ri^2) - ((Ri^2*Ro^2)./r.^2)*(Pi-Po)/((Ro^2-Ri^2));
R = ((Ri):(Ro));

S_theta = []; S_r = [];
% figure();
for r = R
    Ct = cos(theta); St = sin(theta);
    si = interpolateStress(result,r*Ct,r*St);
    S_theta = [S_theta,si.syy*Ct^2 - 2*Ct*St*si.sxy + si.sxx*St^2];
    S_r = [S_r,si.sxx*Ct^2 + 2*Ct*St*si.sxy + si.syy*St^2];
end
f5 = figure('Position',[100,80,800,600],'Renderer','painters');
subplot(1,2,1);
plot(R,S_r,'-*r','LineWidth',0.7); hold on;
plot(R,sr(R),'--k','LineWidth',0.7);
xlabel('$r,\:\mathrm{mm}$');
ylabel('$\sigma_r,\:\mathrm{MPa}$');grid on;
subplot(1,2,2);
plot(R,S_theta,'-*r','LineWidth',0.7); hold on;
plot(R,st(R),'--k','LineWidth',0.7);
xlabel('$r,\:\mathrm{mm}$');
ylabel('$\sigma_\theta,\:\mathrm{MPa}$'); grid on;
legend('FEM','Theory');
sgtitle({['Radial stress variation, $\theta=',num2str(rad2deg(theta)),'^\circ$'],...
    ['$R_i=',num2str(Ri),'\mathrm{mm},\:R_o=',num2str(Ro),...
    '\mathrm{mm},\:P_i=',num2str(Pi),'\mathrm{MPa},\:P_o=',num2str(Po),...
    '\mathrm{MPa},\:L=',num2str(L),'\mathrm{mm}$']});
% print('-f',['Figures/fig03m',num2str(meshsize_max)],'-dsvg')