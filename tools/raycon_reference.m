% raycon_reference.m -- dump machine-precision reference data from the original
% RAYCON MATLAB code for verifying the Julia port in src/raycon/.
%
% Run:  /Applications/MATLAB_R2026a.app/bin/matlab -batch ...
%       "addpath('.../tools/matlab_shadow'); addpath('/Users/jake/PlasmaWorkspace/raycon'); run('.../tools/raycon_reference.m')"
%
% tools/matlab_shadow shadows dmagnetic.m with the FD step matched to the
% Julia port (see tools/matlab_shadow/README.md).
%
% Output: test/reference/raycon_reference.json
%
% Uses the C-Mod ICRF case (data('cmod')) with main.m's overrides:
% kant = [-31.5 -10 0], TYPE='Con', odeDim=4, NRAY=1.

global plasma rays cnst sys

plasma = initPlasma; rays = initRays; cnst = initCnst; sys = initSys;
data('cmod');
rays.TYPE = 'Con'; rays.odeDim = 4; rays.NRAY = 1;
rays.inKspace = 0; rays.time = 0;
plasma.thant = [-.5; .5];
plasma.kant  = [-31.5 -10 0.];

ref = struct();
ref.params = struct('b0',plasma.b0,'r0',plasma.r0,'q0',plasma.q0, ...
    'iaspr',plasma.iaspr,'elong',plasma.elong,'psin',plasma.psin, ...
    'freq',plasma.freq,'omega',plasma.omega,'amass',plasma.amass, ...
    'acharge',plasma.acharge,'n0',plasma.n0,'na',plasma.na,'nb',plasma.nb, ...
    't0',plasma.t0,'ta',plasma.ta,'tb',plasma.tb,'kant',plasma.kant, ...
    'sant',plasma.sant,'MODEL',plasma.MODEL);

% ---------------------------------------------------------------- solovev
disp('SECTION solovev');
% sample (rho,theta) pairs inside the plasma
rhoa = plasma.elong*plasma.iaspr*plasma.r0;
rhos   = [0.03 0.08 0.15 0.22 0.30]*rhoa/0.35;   % spread of minor radii [m]
thetas = [0.001 0.8 2.0 -1.2 3.14];
sol = {};
for k = 1:numel(rhos)
  [sflx,dsdr,dsdz,dsdr2,dsdrz,dsdz2,dsdr3,dsdr2z,dsdrz2,dsdz3] = ...
      solovev(rhos(k),thetas(k),plasma.r0,plasma.iaspr,plasma.elong,0);
  sol{end+1} = struct('rho',rhos(k),'theta',thetas(k),'sflx',sflx, ...
    'dsdr',dsdr,'dsdz',dsdz,'dsdr2',dsdr2,'dsdrz',dsdrz,'dsdz2',dsdz2, ...
    'dsdr3',dsdr3,'dsdr2z',dsdr2z,'dsdrz2',dsdrz2,'dsdz3',dsdz3);
end
ref.solovev = [sol{:}];

% ---------------------------------------------------------------- magnetic
disp('SECTION magnetic');
mag = {};
for k = 1:numel(rhos)
  [b,dbds,dbdt,bp,sflx,dsdr,dsdz,dtdr,dtdz, ...
   ener,enef,enez,eber,ebef,ebez,eper,epef,epez]=magnetic(rhos(k),thetas(k));
  m = struct('rho',rhos(k),'theta',thetas(k),'b',b,'dbds',dbds,'dbdt',dbdt, ...
    'bp',bp,'sflx',sflx,'dsdr',dsdr,'dsdz',dsdz,'dtdr',dtdr,'dtdz',dtdz, ...
    'ener',ener,'enef',enef,'enez',enez,'eber',eber,'ebef',ebef,'ebez',ebez, ...
    'eper',eper,'epef',epef,'epez',epez);
  % second-level outputs
  [b,dbds,dbdt,bp,sflx,dsdr,dsdz,dtdr,dtdz, ...
   ener,enef,enez,eber,ebef,ebez,eper,epef,epez, ...
   dbds2,dbdst,dbdt2,dsdr2,dsdrz,dsdz2,dsdr3,dsdr2z,dsdrz2,dsdz3, ...
   dtdr2,dtdrz,dtdz2]=magnetic(rhos(k),thetas(k));
  m.dbds2=dbds2; m.dbdst=dbdst; m.dbdt2=dbdt2;
  % basis-vector first derivatives (49-output form)
  [b,dbds,dbdt,bp,sflx,dsdr,dsdz,dtdr,dtdz, ...
   ener,enef,enez,eber,ebef,ebez,eper,epef,epez, ...
   dbds2,dbdst,dbdt2,dsdr2,dsdrz,dsdz2,dsdr3,dsdr2z,dsdrz2,dsdz3, ...
   dtdr2,dtdrz,dtdz2, ...
   denerdr,denefdr,denezdr,denerdz,denefdz,denezdz, ...
   deberdr,debefdr,debezdr,deberdz,debefdz,debezdz, ...
   deperdr,depefdr,depezdr,deperdz,depefdz,depezdz]=magnetic(rhos(k),thetas(k));
  m.denerdr=denerdr; m.denezdr=denezdr; m.denerdz=denerdz; m.denezdz=denezdz;
  m.deberdr=deberdr; m.debefdr=debefdr; m.debezdr=debezdr;
  m.deberdz=deberdz; m.debefdz=debefdz; m.debezdz=debezdz;
  m.deperdr=deperdr; m.depefdr=depefdr; m.depezdr=depezdr;
  m.deperdz=deperdz; m.depefdz=depefdz; m.depezdz=depezdz;
  mag{end+1} = m;
end
ref.magnetic = [mag{:}];

% ---------------------------------------------------------------- map flux
disp('SECTION mapflux');
% (r,z) of chosen (s,theta) launch points via modern-fzero equivalent
maps = {};
svals = [0.15 0.4 0.7 0.9]; tvals = [0.001 0.8 2.0 -1.2];
for k = 1:numel(svals)
  f = @(rho) solovev(rho,tvals(k),plasma.r0,plasma.iaspr,plasma.elong,svals(k));
  rho = fzero(f,[1e-12 1.5*rhoa]);
  maps{end+1} = struct('s',svals(k),'theta',tvals(k),'rho',rho, ...
      'r',rho*cos(tvals(k))+plasma.r0,'z',-rho*sin(tvals(k)));
end
ref.mapflux = [maps{:}];

% ---------------------------------------------------------------- disp_eig
disp('SECTION disp_eig');
% phase-space sample points y = [r z kr kz]
ys = [ 0.80  0.05 -31.5  5.0;
       0.75 -0.10 -20.0  0.0;
       0.62  0.02  10.0  3.0;
       0.55  0.12 -31.5 10.0;
       0.70  0.00 -31.5  0.0 ];
de = {};
for k = 1:size(ys,1)
  y4 = ys(k,:);
  d = struct('y',y4);
  d.U    = disp_eig(y4.','Dsp');
  d.mon  = disp_eig(y4.','Mon');
  d.trj  = disp_eig(y4.','Trj').';
  p      = disp_eig(y4.','Pol');
  d.pol_re = real(p).'; d.pol_im = imag(p).';
  y20  = [y4 reshape(eye(4),1,16)].';
  t20  = disp_eig(y20,'Trj');
  d.trj20 = t20.';
  de{end+1} = d;
end
ref.disp_eig = [de{:}];

% also cld3x3 pointwise (same points, switched model)
plasma.MODEL = 'cld3x3';
de3 = {};
for k = 1:size(ys,1)
  y4 = ys(k,:);
  d = struct('y',y4);
  d.U   = disp_eig(y4.','Dsp');
  d.mon = disp_eig(y4.','Mon');
  d.trj = disp_eig(y4.','Trj').';
  y20 = [y4 reshape(eye(4),1,16)].';
  d.trj20 = disp_eig(y20,'Trj').';
  de3{end+1} = d;
end
ref.disp_eig_3x3 = [de3{:}];
plasma.MODEL = 'cld2x2';

% ---------------------------------------------------------------- dispertok
disp('SECTION dispertok');
% pointwise: Dsp, Sgn, Trj, Frq, Msw with synthetic zdot/zddot for eval2nd ops
dt = {};
zdot  = [0.5; -0.3; 2000; 1000];      % synthetic but fixed: same in Julia tests
zddot = [0.1; 0.2; -500; 300];
for k = 1:size(ys,1)
  y4 = ys(k,:).';
  d = struct('y',ys(k,:));
  d.Dsp = dispertok(0.,y4,y4,y4,0,'Dsp');
  d.Sgn = dispertok(0.,y4,y4,y4,0,'Sgn');
  d.Trj = dispertok(0.,y4,y4,y4,0,'Trj').';
  d.Msw = dispertok(0.,y4,y4,y4,0,'Msw');
  d.Frq = dispertok(0.,y4,y4,y4,0,'Frq').';
  d.Mon = dispertok(0.,y4,zdot,zddot,0,'Mon');   % [mon1 mon2 eta2 yg1 yg3]
  d.Sdl = dispertok(0.,y4,zdot,zddot,0,'Sdl');
  dt{end+1} = d;
end
ref.dispertok = [dt{:}];
ref.dispertok_zdot = zdot.'; ref.dispertok_zddot = zddot.';

% ---------------------------------------------------------------- cgamma
disp('SECTION cgamma');
zs = [0 .5+.5i -.5+.5i -.5-.5i .5-.5i 1 1+1i 1i -1+1i -1 -1-1i -1i 1-1i].';
zs = [zs; -1i*0.25; -1i*1.0; -1i*2.7];
% NB upstream cgamma is scalar-only in modern MATLAB (its shift loop
% `for m=1:floor(9-x)` uses the whole vector x; old MATLAB silently took the
% first element, so vector calls were also historically wrong): map per element.
cg = arrayfun(@(zz) cgamma(zz), zs);
ref.cgamma = struct('z_re',real(zs).','z_im',imag(zs).', ...
                    'g_re',real(cg).','g_im',imag(cg).');

% ------------------------------------------------- full conversion pipeline
disp('SECTION conversion');
% Reproduce ray.m convert_list order at a plausible near-conversion point:
% Mon (stores yalf0+eta2est) -> Sdl iterate -> Trs -> Cnv (tau/beta).
try
  y0 = [0.62 0.02 10.0 3.0].';   % near the mode-conversion layer for cmod
  rays.yalf0 = []; rays.eta2est = [];
  mon = dispertok(0,y0,zdot,zddot,0,'Mon');
  zst = dispertok(0,y0,zdot,zddot,0,'Sdl');
  for it = 1:30
    zst_new = dispertok(0,zst.',zdot,zddot,0,'Sdl');
    if norm((zst_new-zst)./zst_new,inf) < 1e-4, zst = zst_new; break; end
    zst = zst_new;
  end
  conv = struct('y0',y0.','mon',mon,'zst',zst);
  try
    ytrs = dispertok(0,zst.',zdot,zddot,0,'Trs');
    conv.ytrs = ytrs.';
  catch etrs
    conv.ytrs_error = etrs.message;
  end
  rays.tau=[]; rays.beta=[];
  ycnv = dispertok(0,zst.',zdot,zddot,0,'Cnv');
  conv.ycnv = ycnv.';
  conv.tau = rays.tau;
  conv.beta_re = real(rays.beta); conv.beta_im = imag(rays.beta);
  ref.conversion = conv;
catch econv
  ref.conversion_error = econv.message;
end

% ------------------------------------------------- end-to-end cmod ray
disp('SECTION ray');
% Try the legacy ODE stack (may fail on R2026a); best effort.
try
  rays.timespan = 5e-2; rays.timeintv = rays.timespan;
  rays.initialstep = 1e-7*rays.timespan;
  th0 = 0.001;
  rays.sray0 = plasma.sant; rays.thray0 = th0;
  f = @(rho) solovev(rho,th0,plasma.r0,plasma.iaspr,plasma.elong,rays.sray0);
  rho = fzero(f,[1e-12 1.5*rhoa]);
  r = rho*cos(th0)+plasma.r0; z = -rho*sin(th0);
  y = adjust_disp_m([r,z,plasma.kant(1),plasma.kant(3)],0);
  ref.launch = struct('r',r,'z',z,'y0',y.');
  rays.y = y; rays.MON=[]; rays.TR=[]; rays.YR=[]; rays.stp=0;
  rays.tspan = [0 rays.timespan];
  [tspan,yy,odeOptions] = inittok(rays.timespan,y);
  [tr,yr] = ode45('trajectory',tspan,yy,odeOptions);
  ref.ray = struct('t',tr.','y',yr);
catch eode
  ref.ray_error = eode.message;
end

% ---------------------------------------------------------------- write
txt = jsonencode(ref);
fid = fopen('/Users/jake/PlasmaWorkspace/HybridPlasmaPIC.jl/test/reference/raycon_reference.json','w');
fwrite(fid,txt); fclose(fid);
disp('raycon_reference.json written');
