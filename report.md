---
author:
- Mohammad Abazari
title: Thick cylinder under pressure
---

# Lame's problem-Thick cylinder subjected to internal pressure {#sec_theory}

Consider a thick cylinder of inner radius $R_i$, outer radius $R_o$,
length $L$ subjected to internal pressure $P_i$ and outer pressure
$P_o=0$. Two cases of plane-stress $\sigma_z=0$ and -strain
$\varepsilon_z=0$ are studied.

## Plane stress

Assuming that both cylinder ends are free thus $\sigma_z=0$, $$\nonumber
\frac{\partial\sigma_r}{\partial r}+\frac{\sigma_r-\sigma_\theta}{r}=0$$
$r$ is the only independent variable in this expression which could be
rewriten as below $$\label{eq1}
\frac{\mathrm{d}}{\mathrm{d}r}(r\sigma_r)-\sigma_\theta=0$$ Following
hoek's law, $$\nonumber\begin{aligned}
\varepsilon_r =& \frac{1}{E}\left(\sigma_r-\nu\sigma_\theta\right)\\
\varepsilon_\theta =& \frac{1}{E}\left(\sigma_r-\nu\sigma_\theta\right)\\
\sigma_r =& \frac{E}{1-\nu^2}\left(\varepsilon_r+\nu\varepsilon_\theta\right)\\
\sigma_\theta =& \frac{E}{1-\nu^2}\left(\varepsilon_\theta-\nu\varepsilon_r\right)
\end{aligned}$$ Substituting strain equations,
$$\label{eq2}\begin{aligned}
\sigma_r =& \frac{E}{1-\nu^2}\left(\frac{\mathrm{d}u_r}{\mathrm{d}r}+\nu\frac{u_r}{r}\right)\\
\sigma_\theta =& \frac{E}{1-\nu^2}\left(\frac{u_r}{r}+\nu\frac{\mathrm{d}u_r}{\mathrm{d}r}\right)
\end{aligned}$$ Substituting the above in
[\[eq1\]](#eq1){reference-type="ref" reference="eq1"} $$\nonumber
\begin{aligned}
\frac{\mathrm{d}}{\mathrm{d}r}\left(r\frac{\mathrm{d}u_r}{\mathrm{d}r}+vu_r\right)-\left(\frac{u_r}{r}+v\frac{\mathrm{d}u_r}{\mathrm{d}r}\right)=&0\\
\frac{\mathrm{d}u_r}{\mathrm{d}r}+r\frac{\mathrm{d}^2u_r}{\mathrm{d}r^2}+\nu\frac{\mathrm{d}u_r}{\mathrm{d}r}-\frac{u_r}{r}-\nu\frac{\mathrm{d}u_r}{\mathrm{d}r}=&0\\
\frac{\mathrm{d}^2u_r}{\mathrm{d}r^2}+\frac{1}{r}\frac{\mathrm{d}u_r}{\mathrm{d}r}-\frac{u_r}{r^2}=&0\\
\frac{\mathrm{d}}{\mathrm{d}r}\left[\frac{1}{r}\frac{\mathrm{d}}{\mathrm{d}r}(u,r)\right]=&0
\end{aligned}$$ Assuming $u_r$ follows the below function
$$\label{eqs03}
u_r=C_1r+\frac{C_2}{r}$$ Substituting the above in
[\[eq2\]](#eq2){reference-type="ref" reference="eq2"}, $$\begin{aligned}
\sigma_r =& \frac{E}{1-\nu^2}\left[C_1(1+\nu)-C_2(1-\nu)\frac{1}{r^2}\right]\\
\sigma_\theta =& \frac{E}{1-\nu^2}\left[C_1(1+\nu)+C_2(1-\nu)\frac{1}{r^2}\right]
\end{aligned}$$

Constants $C_1$ and $C_2$ are applying boundary conditions. $$\nonumber
\begin{aligned}
\sigma_r(r=R_i) = &-P_i= \frac{E}{1-\nu^2}\left[C_1(1+\nu)-C_2(1-\nu)\frac{1}{R_i^2}\right]\\
\sigma_r(r=R_o) = &-P_o=0=\frac{E}{1-\nu^2}\left[C_1(1+\nu)-C_2(1-\nu)\frac{1}{R_o^2}\right]\\
\end{aligned}$$ Solving for these two constants, $$\nonumber
\begin{aligned}
C_1=&\frac{1-\nu}{E}\frac{P_iR_i^2-P_oR_o^2}{R_o^2-R_i^2}\\
C_2=&\frac{1+\nu}{E}\frac{R_i^2-R_o^2}{R_o^2-R_i^2}\left(P_i-P_o\right)
\end{aligned}$$ Substituting these constants into the above equations
$$\label{eqs04}
\begin{aligned}
\sigma_{r}=&\frac{P_i R_i^{2}-P_o R_o^{2}}{R_o^{2}-R_i^{2}}-\frac{R_i^{2} R_o^{2}}{r^{2}} \frac{P_i-P_o}{R_o^{2}-R_i^{2}}\\
\sigma_{\theta}=&\frac{P_i R_i^{2}-P_o R_o^{2}}{R_o^{2}-R_i^{2}}+\frac{R_i^{2} R_o^{2}}{r^{2}} \frac{P_i-P_o}{R_o^{2}-R_i^{2}}
\end{aligned}$$

### Plane stress problem with internal pressure

If outside pressure $P_o=0$, then $$\begin{aligned}
\sigma_{r}=&\frac{P_i R_i^{2}}{R_o^{2}-R_i^{2}}\left[1-\frac{R_o^2}{r^2}\right]\label{eqs01}\\
\sigma_{\theta}=&\frac{P_i R_i^{2}}{R_o^{2}-R_i^{2}}\left[1+\frac{R_o^2}{r^2}\right]\label{eqs02}\\\end{aligned}$$
The above derivation shows that $\sigma_r$ is compressive through the
cylinder thickness while and $\sigma_\theta$ is tensile and positive.

## Plane strain

In the plane strain case, $\sigma_z$ is assumed constant,
[\[eq1\]](#eq1){reference-type="ref" reference="eq1"} $$\nonumber
\frac{\mathrm{d}}{\mathrm{d}r}(r\sigma_r)-\sigma_\theta=0$$ From Hoek's
law $$\nonumber
\begin{aligned}
\varepsilon_{r}=&\frac{1}{E}\left[\sigma_{r}-v\left(\sigma_{\theta}+\sigma_{z}\right)\right]\\
\varepsilon_{\theta}=&\frac{1}{E}\left[\sigma_{\theta}-v\left(\sigma_{r}+\sigma_{z}\right)\right]\\
\varepsilon_{z}=&\frac{1}{E}\left[\sigma_{r}-v\left(\sigma_{r}+\sigma_{\theta}\right)\right]
\end{aligned}$$ With $\varepsilon_z=0$ $$\nonumber
\begin{aligned}
\sigma_{z}=&v\left(\sigma_{r}+\sigma_{\theta}\right)\\
\varepsilon_{r}=&\frac{1+v}{E}\left[(1-v) \sigma_{r}-v \sigma_{\theta}\right]\\
\varepsilon_{\theta}=&\frac{1+v}{E}\left[(1-v) \sigma_{\theta}-v \sigma_{r}\right]
\end{aligned}$$ Solving for stress components $$\nonumber
\begin{aligned}
\sigma_{\theta}=&\frac{E}{(1-2 v)(1+v)}\left[v \varepsilon_{r}+(1-v) \varepsilon_{\theta}\right]\\
\sigma_{r}=&\frac{E}{(1-2 v)(1+v)}\left[(1-v) \varepsilon_{r}+v \varepsilon_{\theta}\right]
\end{aligned}$$ Substituting strain $$\label{eq4}
\begin{aligned}
\sigma_{r}=&\frac{E}{(1-2 v)(1+v)}\left[(1-v) \frac{d u_{r}}{d r}+v \frac{u_{r}}{r}\right]\\
\sigma_{\theta}=&\frac{E}{(1-2 v)(1+v)}\left[v \frac{d u_{r}}{d r}+(1-v) \frac{u_{r}}{r}\right]\\
\end{aligned}$$ Substituting the above in the equilibrum equations
$$\nonumber\begin{aligned}
\frac{\mathrm{d}}{\mathrm{d} r}\left[(1-v) r \frac{d u_{r}}{d r}+v u_{r}\right]-v \frac{d u_{r}}{d r}-(1-v) \frac{u_{r}}{r}&=0\\
\frac{d u_{r}}{d r}+r \frac{\mathrm{d}^{2} u_{r}}{\mathrm{d} r^{2}}-\frac{u_{r}}{r}&=0\\
\frac{\mathrm{d}}{\mathrm{d} r}\left(\frac{\mathrm{d} u}{\mathrm{d} r}+\frac{u_{r}}{r}\right)&=0
\end{aligned}$$ Assuming [\[eqs03\]](#eqs03){reference-type="ref"
reference="eqs03"} for $u_r$ and substituting into
[\[eq4\]](#eq4){reference-type="ref" reference="eq4"}
$$\nonumber\begin{aligned}
\sigma_{\theta}=&\frac{E}{(1-2 v)(1+v)}\left[C_{1}+(1-2 v) \frac{C_{2}}{r^{2}}\right]\\
\sigma_{r}=&\frac{E}{(1-2 v)(1+v)}\left[C_{1}-(1-2 v) \frac{C_{2}}{r^{2}}\right]
\end{aligned}$$ Again applying boundary conditions $$\nonumber
\begin{aligned}
\sigma_r(r=R_i)=&-P_i=\frac{E}{(1-2 v)(1+v)}\left[C_{1}-(1-2 v) \frac{C_{2}}{R_i^{2}}\right]\\
\sigma_r(r=R_o)=&-P_{o}=\frac{E}{(1-2 v)(1+v)}\left[C_{1}+(1-2 v) \frac{C_{2}}{R_o^{2}}\right]
\end{aligned}$$ Thus, $$\nonumber\begin{aligned}
C_{1}=&\frac{(1-2 v)(1+v)}{E} \frac{P_{o} R_o^{2}-P_{i} R_i^{2}}{R_i^{2}-R_o^{2}}\\
C_{2}=&\frac{1+v}{E} \frac{\left(P_{o}-P_{i}\right) R_i^{2} R_o^{2}}{R_i^{2}-R_o^{2}}
\end{aligned}$$ Substituting these constants into the above,
$$\nonumber\begin{aligned}
\sigma_{r}=&\frac{P_{i} R_i^{2}-P_{o} R_o^{2}}{R_o^{2}-R_i^{2}}-\frac{R_i^{2} R_o^{2}}{r^{2}} \frac{P_{i}-P_{o}}{R_o^{2}-R_i^{2}} \\
\sigma_{\theta}=&\frac{P_{i} R_i^{2}-P_{o} R_o^{2}}{R_o^{2}-R_i^{2}}+\frac{R_i^{2} R_o^{2}}{r^{2}} \frac{P_{i}-P_{o}}{R_o^{2}-R_i^{2}}
\end{aligned}$$ Which is equals [\[eqs04\]](#eqs04){reference-type="ref"
reference="eqs04"}.

# MATLAB modeling {#sec_numm}

A cylinder of inner radius of $R_i=40\mathrm{mm}$ outer radius
$R_o=80\mathrm{mm}$ length $L=1000\mathrm{mm}$ inside pressure
$17\mathrm{MPa}$ outside pressure $2$ was meshed with minimum element
dimension of $1\mathrm{mm}$ and maximum element dimension of
$10\mathrm{mm}$.

![3D model geometry.](Figures/fig01.pdf){#fig01 width="40%"}

![2D model and mesh.](Figures/fig01m8.pdf){#fig01m width="\\textwidth"}

![2D results.](Figures/fig03m8.pdf){#fig05a width="80%"}

![3D results.](Figures/fig05m10.pdf){#fig05b width="80%"}

![Stress variation in length.](Figures/fig06m10.pdf){#fig06 width="80%"}

# Conclusions

-   Both 2D and 3D numerical models ([2](#sec_numm){reference-type="ref"
    reference="sec_numm"}) agree with theoretical derivation
    ([1](#sec_theory){reference-type="ref" reference="sec_theory"})
    ([\[fig05\]](#fig05){reference-type="ref" reference="fig05"}).

-   Based on both numerical and theoretical results radial ($\sigma_r$)
    and angular stresses ($\sigma_\theta$) reduce
    radially([\[fig05\]](#fig05){reference-type="ref"
    reference="fig05"}).

-   3D numerical results shows no significant change in longitudinal
    stress ($\sigma_z$) which sattisfies the plane stress and strain
    assumptions([5](#fig06){reference-type="ref" reference="fig06"}).

-   Longitudinal stress oscilations in the cylinder length
    ([5](#fig06){reference-type="ref" reference="fig06"}) showcases
    model sensitivity to proper boundary conditions which is dealt with
    through finer mesh sizes.

# Appendices

## Stress transformation from cartesian to cylindrical coordinates

## 3D model

``` {.matlab style="Matlab-editor"}
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

% % Problem statement =====================================================
% Geometry description ===================================================
Ri = 40; % Inner radius Ri, mm
Ro = 80; % Outer radius Ro, mm
L = 1000; % Length L, mm
Pi = 17; % Inside pressure (radial pressure), N/mm2
Po = 0; % Outside pressure (radial pressure), N/mm2
meshsize_max = 20; % maximum mesh dimension, mm
meshsize_min = 5; % minimum mesh dimension, mm
mesh_order = 'linear'; % or quadratic
% Material Properties ====================================================
E = 210e3; % Modulus of elasticity E, N/mm2
nu = 0.3; % Poisson's ratio \nu
%% Model creation =========================================================
model = createpde('structural','static-solid');
% Introduce geometry =====================================================
gm = multicylinder([Ri Ro],L,'Void',[true,false]);
model.Geometry = gm;
figure
pdegplot(model,'CellLabels','on','FaceLabels','on','EdgeLabels',...
    'on','VertexLabels','on','FaceAlpha',0.5)
view(30,30);
title('Cylinder with Face/Edge/Cell/Vertex Labels')

% Assign Material Values =================================================
structuralProperties(model,'YoungsModulus',E, ...
                           'PoissonsRatio',nu);
% Applying boundary conditions ===========================================
structuralBC(model,'Face',1,'Constraint','symmetric');
structuralBC(model,'Face',2,'Constraint','symmetric');

% Apply pressure loads in and out ========================================
structuralBoundaryLoad(model,"Face",3,"Pressure",Pi);
structuralBoundaryLoad(model,"Face",4,"Pressure",Po);

% Mesh description, Generate mesh ========================================
generateMesh(model,'Hmax',meshsize_max,'Hmin',meshsize_min,...
    'GeometricOrder',mesh_order);
figure;
pdeplot3D(model);
title(['Mesh with Quadratic Tetrahedral Elements with $\Delta_m=',...
    num2str(meshsize_max),'$mm']);

% Calculate solution =====================================================
result = solve(model);

%% Evaluate Results =======================================================
figure
pdeplot3D(model,'ColorMapData',result.VonMisesStress)
title('Von-Mises Stress, MPa')
colormap('jet')
figure
pdeplot3D(model,'ColorMapData',result.Stress.yy)
title('Von-Mises Stress, MPa')
colormap('jet')
save(['ResultsM',num2str(meshsize_max),'m',...
    num2str(meshsize_min),mesh_order,'.mat'],'result','model')
%% Stress values by radius ================================================
% Theory
st = @(r) (Pi*Ri^2-Po*Ro^2)/(Ro^2-Ri^2) ...
    + ((Ri^2*Ro^2)./r.^2)*(Pi-Po)/((Ro^2-Ri^2));
sr = @(r) (Pi*Ri^2-Po*Ro^2)/(Ro^2-Ri^2) ...
    - ((Ri^2*Ro^2)./r.^2)*(Pi-Po)/((Ro^2-Ri^2));
% Based on cartesian to cylindrical transformation if Cos(theta) =
% Sine(theta) = 0 then Sxx = Syy and Stt = Syy as in the present case.
% http://solidmechanics.org/text/AppendixD/AppendixD.htm

% Practice and FEM =======================================================
R = Ri:Ro;

S_theta = []; S_r = [];

for r = R
    si = interpolateStress(result,r,0,L/2);
    S_theta = [S_theta,si.syy];
    S_r = [S_r,si.sxx];
end
%% Longitudinal stress values =============================================
figure;
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

%% Visualize results ======================================================
f5 = figure;
subplot(1,2,1);
plot(R,S_r,'b','LineWidth',0.7); hold on;
plot(R,sr(R),'--k','LineWidth',0.7);
xlabel('$r,\:\mathrm{mm}$');
ylabel('$\sigma_r,\:\mathrm{MPa}$');grid on;
subplot(1,2,2);
plot(R,S_theta,'b','LineWidth',0.7); hold on;
plot(R,st(R),'--k','LineWidth',0.7);
xlabel('$r,\:\mathrm{mm}$');
ylabel('$\sigma_\theta,\:\mathrm{MPa}$'); grid on;
legend('FEM','Theory');
sgtitle({'Radial stress variation',...
    ['$R_i',num2str(Ri),'\mathrm{mm},\:R_o=',num2str(Ro),...
    '\mathrm{mm},\:P_i=',num2str(Pi),'\mathrm{MPa},\:P_o=',num2str(Po),...
    '\mathrm{MPa},\:L=',num2str(L),'\mathrm{mm}$']});
f6 = figure;
plot(l,Sz_Ro,'-b','LineWidth',0.7); hold on;
plot(l,Sz_Ri,'--r','LineWidth',0.7); 
plot(l,Sz_Rm,'-.k','LineWidth',0.7);
xlabel('$L,\:\mathrm{mm}$');
ylabel('$\sigma,\:\mathrm{MPa}$');grid on;
title({'Longitudinal stress variation',...
    ['$R_i',num2str(Ri),'\mathrm{mm},\:R_o=',num2str(Ro),...
    '\mathrm{mm},\:P_i=',num2str(Pi),'\mathrm{MPa},\:P_o=',num2str(Po),...
    '\mathrm{MPa},\:L=',num2str(L),'\mathrm{mm}$']});
```

## 2D model

``` {.matlab style="Matlab-editor"}
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
generateMesh(model,'Hmax',meshsize_max,'Hmin',meshsize_min,'GeometricOrder',mesh_order);
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
save(['Results2dM',num2str(meshsize_max),'m',num2str(meshsize_min),mesh_order,'.mat'],'result','model')
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
```

## 2D geometry generation with ABAQUS
