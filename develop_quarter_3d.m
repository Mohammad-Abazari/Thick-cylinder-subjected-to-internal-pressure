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
meshsize_max = 5; % maximum mesh dimension, mm
meshsize_min = 1; % minimum mesh dimension, mm
mesh_order = 'linear'; % or quadratic
% Material Properties ====================================================
E = 210e3; % Modulus of elasticity E, N/mm2
nu = 0.3; % Poisson's ratio \nu
% rho = 8000e-9; % Mass density \rho, kg/mm3 irrelevent for static model

% Measurement
theta = pi/4;

%% Model creation =========================================================
model = createpde('structural','static-solid');
% Introduce geometry =====================================================
importGeometry(model,'mesh3d_quarter.stl');
f1 = figure('Position',[100,40,1500,900],'Renderer','painters');
subplot(1,2,1);
pdegplot(model,'CellLabels','on','FaceLabels','on','EdgeLabels','on','VertexLabels','on','FaceAlpha',0.5)
title('Cylinder with Face/Edge/Cell/Vertex Labels')

% Assign Material Values
structuralProperties(model,'YoungsModulus',E, ...
                           'PoissonsRatio',nu);
% Applying boundary conditions
structuralBC(model,'Face',1,'Constraint','symmetric');
structuralBC(model,'Face',3,'Constraint','symmetric');
structuralBC(model,'Face',5,'Constraint','symmetric');
structuralBC(model,'Face',6,'Constraint','symmetric');

% Apply pressure loads in and out
structuralBoundaryLoad(model,"Face",4,"Pressure",Pi);
structuralBoundaryLoad(model,"Face",2,"Pressure",Po);

% Mesh description, Generate mesh
mesh = generateMesh(model,'Hmax',meshsize_max,'Hmin',meshsize_min,'GeometricOrder',mesh_order);
% subplot(1,2,2);
% pdeplot3D(model);
% title(['Mesh with Quadratic Tetrahedral Elements with $\Delta_m=',num2str(meshsize_max),'$mm']);
% f2 = figure('Position',[100,40,1500,900],'Renderer','painters');
% subplot(1,2,1);
% pdeplot3D(model,'NodeLabels','on'); title('Cross-Section'); view(0,90);
% % Specify line of nodes considered for analysis
R = ((Ri):(Ro)); n_ids = zeros(size(R));
for i = 1:length(R)
    r = R(i);
    point = r*[cos(theta);sin(theta);L/r/2];
    n_ids(i) = findNodes(mesh,'nearest',point);
end
% subplot(1,2,2);
% pdemesh(model,'NodeLabels','off'); hold on;
% plot(mesh.Nodes(1,n_ids),mesh.Nodes(2,n_ids),'or','MarkerFaceColor','g')
% axis off;
% view(0,-90)
%% Calculate solution
result = solve(model);

%% Evaluate Results
figure
pdeplot3D(model,'ColorMapData',result.VonMisesStress)
title('Von-Mises Stress, MPa')
colormap('jet')
figure
pdeplot3D(model,'ColorMapData',result.Stress.yy)
title('Von-Mises Stress, MPa')
colormap('jet')
save(['ResultsM',num2str(meshsize_max),'m',num2str(meshsize_min),mesh_order,'.mat'],'result','model')
%% Stress values by radius
% Theory
st = @(r) (Pi*Ri^2-Po*Ro^2)/(Ro^2-Ri^2) + ((Ri^2*Ro^2)./r.^2)*(Pi-Po)/((Ro^2-Ri^2));
sr = @(r) (Pi*Ri^2-Po*Ro^2)/(Ro^2-Ri^2) - ((Ri^2*Ro^2)./r.^2)*(Pi-Po)/((Ro^2-Ri^2));
% Based on cartesian to cylindrical transformation if Cos(theta) =
% Sine(theta) = 0 then Sxx = Syy and Stt = Syy as in the present case.
% http://solidmechanics.org/text/AppendixD/AppendixD.htm

% Practice and FEM
R = Ri:Ro;

S_theta = []; S_r = [];

for r = R
    si = interpolateStress(result,r,0,L/2);
    S_theta = [S_theta,si.syy];
    S_r = [S_r,si.sxx];
end
%% Longitudinal stress values
Sz_Ro = [];
Sz_Ri = [];
Sz_Rm = [];
Rm = (Ro+Ri)/2;

l = 0:meshsize_max:L;

for li = l
    si_Ro = interpolateStress(result,Ro,0,li);
    si_Ri = interpolateStress(result,Ri,0,li);
    si_Rm = interpolateStress(result,Rm,0,li);
    Sz_Ro = [Sz_Ro,si_Ro.szz];
    Sz_Ri = [Sz_Ri,si_Ri.szz];
    Sz_Rm = [Sz_Rm,si_Rm.szz];
end

%% Visualize results
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
sgtitle({'Radial stress variation',...
    ['$R_i=',num2str(Ri),'\mathrm{mm},\:R_o=',num2str(Ro),...
    '\mathrm{mm},\:P_i=',num2str(Pi),'\mathrm{MPa},\:P_o=',num2str(Po),...
    '\mathrm{MPa},\:L=',num2str(L),'\mathrm{mm}$']});
% print('-f',['Figures/fig05m',num2str(meshsize_max)],'-dsvg')
% matlab2tikz('figurehandle',f5,'filename','Figures/fig05++.tex' ,'standalone', false,'floatFormat','%.4g')
f6 = figure('Position',[100,80,800,600],'Renderer','painters');
plot(l,Sz_Ro,'-b','LineWidth',0.7); hold on;
plot(l,Sz_Ri,'--r','LineWidth',0.7); 
plot(l,Sz_Rm,'-.k','LineWidth',0.7);
xlabel('$L,\:\mathrm{mm}$');
ylabel('$\sigma,\:\mathrm{MPa}$');grid on;
title({'Longitudinal stress variation',...
    ['$R_i=',num2str(Ri),'\mathrm{mm},\:R_o=',num2str(Ro),...
    '\mathrm{mm},\:P_i=',num2str(Pi),'\mathrm{MPa},\:P_o=',num2str(Po),...
    '\mathrm{MPa},\:L=',num2str(L),'\mathrm{mm}$']});
% print('-f',['Figures/fig06m',num2str(meshsize_max)],'-dsvg')
% matlab2tikz('figurehandle',f6,'filename','Figures/fig06++.tex' ,'standalone', false,'floatFormat','%.4g')