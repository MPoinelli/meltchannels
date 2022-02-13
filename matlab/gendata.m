clear all, close all, clc

% Choose rift lenght [0,1,2,3]

rift = 0; % km

% Dimensions of grid
nx=200;
ny=320;
nz=60;
deltaZ = 10;

dx = 250;
dy = 250;

% Nominal depth of model (meters)
H = -600;
Hmin = -250;
Hcav = -200;

% RIFT
Rift_1 = 150;
Rift_2 = 153; % 2km width

%% Draw Geometry
icetopo_profile = zeros(nx,1);
icetopo_profile(1:70)=linspace(Hmin,Hcav,70);
icetopo_profile(71:200)= Hcav;

icetopo_profile(Rift_1:Rift_2) = -10;

icetopo = zeros(nx, ny);
for jj = 1:ny % y levels are the same 
    icetopo(:,jj)=icetopo_profile;
end

% % closing partial south-north border
% icetopo(1:344,1:20)=0;
% icetopo(1:344,end-19:end)=0;

% frontdepth = zeros(nx, ny);
% frontdepth(341, 21:end-20) = -icetopo(340, 21:end-20);
% frontcirc = zeros(nx, ny);
% frontcirc(341, 21:end-20) = deltaZ./dy;
% 
% if rift ~= 0 
%     frontdepth(Rift_1, 21:end-20) = -icetopo(Rift_1 - 1, 21:end-20)-10;
%     frontdepth(Rift_2, 21:end-20) = -icetopo(Rift_2 + 1, 21:end-20)-10;
%     frontcirc(Rift_1, 21:end-20) = deltaZ./dy;
%     frontcirc(Rift_2, 21:end-20) = deltaZ./dy;
% end

bathy = -600.* ones(nx,ny);
% bathy_profile(1:200)=linspace(Hmin-20,-600, 200);
% 
% bathy = H .* ones(nx,ny);
% for jj = 1:ny % y levels are the same 
%     bathy(:,jj)=bathy_profile;
% end

% closing partial East-south-north border
bathy(:, 1:2)=0;
bathy(:, end-1:end)=0;

bathy(1:2,:) = 0 ;

writebin('bathy.bin',bathy);
writebin('draft.bin',icetopo);
% writebin(strcat('frontdepth_rift',num2str(rift),'_400x800.bin'),frontdepth);
% writebin(strcat('frontcirc_rift',num2str(rift),'_400x800.bin'),frontcirc);

%% Compute Pressure Anomaly
% After modifying the code in calc_phi_hyd.F on Apr26,2012 this is the
% consistent way of computing phi0surf. For this, we need the grid
% information (hFacC's). For convenience, it's taken from a previous model
% run.
%
% The way of computing phi0surf consistent with code prior to Apr26,2012
% is recovered by setting drloc*dphi=0

% % Read hfacc
% if rift == 0
%     hfac=readbin('/nobackupp11/mpoinell/testCaseRift/run_testCaseRift00a.1/hFacC.data',[nx ny nz]); 
% else
%     hfac=readbin(strcat('/nobackupp11/mpoinell/testCaseRift/run_testCaseRift00c.',num2str(rift),'/hFacC.data'),[nx ny nz]);
% end
% 
% % read temp_salinity profiles
% theta = readbin('THETA_400x800x60.bin',[nx,ny,nz]); t = squeeze(theta(1,1,:));
% salt  = readbin('SALT_400x800x60.bin',[nx,ny,nz]);  s = squeeze(salt(1,1,:));
% 
% 
% % compute potential field underneath ice shelf
% gravity=9.81;
% rhoConst = 1030;
% 
% dz = deltaZ*ones(1,nz);
% zgp1 = [0,cumsum(dz)];
% zc = .5*(zgp1(1:end-1)+zgp1(2:end));
% zg = zgp1(1:end-1);
% dz = diff(zgp1);
% 
% %eos = 'linear';
% eos = 'jmd95z';
% %eos = 'mdjwf';
% 
% talpha = 0;
% sbeta  = 1e-4;
% 
% %tref = -1.9*ones(nz,1);
% %t = tref;
% %sref = 34.4*ones(nz,1);
% %s = sref;
% 
% k=1;
% dzm = abs([zg(1)-zc(1) .5*diff(zc)]);
% dzp = abs([.5*diff(zc) zc(end)-zg(end)]);
% p = abs(zc)*gravity*rhoConst*1e-4;
% dp = p;
% kp = 0;
% 
% while std(dp) > 1e-13
%   phiHydF(k) = 0;
%   p0 = p;
%   kp = kp+1
%   for k = 1:nz
%     switch eos
%      case 'linear'
%       drho = rhoConst*(1-talpha*(t(k)-tref(k))+sbeta*(s(k)-sref(k)))-rhoConst;
%      case 'jmd95z'
%       drho = densjmd95(s(k),t(k),p(k))-rhoConst;
%      case 'mdjwf'
%       drho = densmdjwf(s(k),t(k),p(k))-rhoConst;
%      otherwise
%       error(sprintf('unknown EOS: %s',eos))
%     end
%     phiHydC(k)   = phiHydF(k) + dzm(k)*gravity*drho/rhoConst;
%     phiHydF(k+1) = phiHydC(k) + dzp(k)*gravity*drho/rhoConst;
%   end
%   switch eos
%    case 'mdjwf'
%     p = (gravity*rhoConst*abs(zc) + phiHydC*rhoConst)/gravity/rhoConst;
%   end
%   dp = p-p0;
% end
% 
% msk=sum(hfac,3); msk(msk>0)=1;
% phi0surf = zeros(nx,ny);
% 
% for ix=1:nx
%   for iy=1:ny
%     k=max(find(abs(zg)<abs(icetopo(ix,iy))));
%     if isempty(k)
%       k=0;
%     end
%     if k>0
%       kp1=min(k+1,nz);
%       drloc=1-hfac(ix,iy,k);
%       %drloc=(abs(icetopo(ix,iy))-abs(zg(k)))/dz(k);
%       dphi = phiHydF(kp1)-phiHydF(k);
%       phi0surf(ix,iy) = (phiHydF(k)+drloc*dphi)*rhoConst*msk(ix,iy);
%     end
%   end
% end
% 
% figure, mypcolor(bathy'),colorbar
% figure, mypcolor(icetopo'),colorbar
% figure, mypcolor(frontcirc'),colorbar
% figure, mypcolor(frontdepth'),colorbar
% figure, mypcolor(phi0surf'),colorbar
% 
% writebin(strcat('phi0surf_rift',num2str