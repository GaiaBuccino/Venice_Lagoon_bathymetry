clear all
close all

% grigliato regolare: mesh completa con profondità zero nei punti terra
filein='/share/data/squerin/COBRA/storage/home/squerin/SPILL/Batimetria/RUOTATORE_Adriatic_HV/Batimetria.dat';
load(filein,'-ascii')

% ADRIA02
bat=Batimetria;
lat=bat(1:1921:end,2);ly=size(lat,1);
lon=bat(1:1921,1);lx=size(lon,1);
% figure;scatter(bat(:,1),bat(:,2),1,bat(:,3));

% punti terra
% bat(bat<=0.0)=0.0;

% nuova griglia a 1/128
nx=792; ny=424;
res=1/128;

xc0=12.22265625; xce=xc0+(nx-1)*res;
yc0=42.47265625; yce=yc0+(ny-1)*res;

mask=-2000;

% 1/128
[Y,X]=meshgrid(yc0:res:yce,xc0:res:xce);

% ADRIA02
[LAT,LON]=meshgrid(lat,lon);
tmp=reshape(-bat(:,3),lx,ly);
bat01=interp2(LAT,LON,tmp,Y,X);
% bat01=interp2(LAT,LON,tmp,Y,X,'nearest');

% bat01=griddata(lat,lon,bat(:,3),X,Y);
% bat02=griddata(bat(:,1),bat(:,2),bat(:,3),X,Y,'nearest');
% figure;pcolor(bat01');shading flat;colorbar
% tanto per sicurezza...
bat01(bat01>=0.0)=0.0;
% elimino zone basse
bat01(bat01>=-0.5)=0.0;
% abbasso il resto a -3.001
bat01(bat01>=-3.001 & bat01<-0.5)=-3.001;
% spiana tutto a -278.423 (bottom del livello 30)
maxdepth1=min(min(bat01));
bat01(bat01<=-278.423)=-278.423;

% Mascherone ad Est della punta dell'Istria
for i=205:225
    for j=200:ny
        if (bat01(i,j)>=-17.503 && bat01(i,j)<0.0)
            bat01(i,j)=-17.503;  % mask
        elseif (bat01(i,j)>=-21.842 && bat01(i,j)<-17.503)
            bat01(i,j)=-21.842;  % mask
        end
    end
end
for i=226:nx
    for j=200:ny
        if (bat01(i,j)>=-26.507 && bat01(i,j)<0.0)
            bat01(i,j)=-26.507;  % mask
        elseif (bat01(i,j)>=-31.519 && bat01(i,j)<-26.507)
            bat01(i,j)=-31.519;  % mask
        end
    end
end
for i=300:nx
    for j=1:ny
        if (bat01(i,j)>=-26.507 && bat01(i,j)<0.0)
            bat01(i,j)=-26.507;  % mask
        elseif (bat01(i,j)>=-31.519 && bat01(i,j)<-26.507)
            bat01(i,j)=-31.519;  % mask
        end
    end
end

% Maschera sul resto della costa istriana
for i=150:204
    for j=300:388
        if (bat01(i,j)>=-17.503 && bat01(i,j)<0.0)
            bat01(i,j)=-17.503;  % mask
        elseif (bat01(i,j)>=-21.842 && bat01(i,j)<-17.503)
            bat01(i,j)=-21.842;  % mask
        end
    end
end

% Maschera da Grignano a "Pirano"
for i=170:210
    for j=389:405
        if (bat01(i,j)>=-13.468 && bat01(i,j)<0.0)
            bat01(i,j)=-13.468;  % mask
        elseif (bat01(i,j)>=-17.503 && bat01(i,j)<-13.468)
            bat01(i,j)=-17.503;  % mask
        end
    end
end
for i=185:210
    for j=406:410
        if (bat01(i,j)>=-13.468 && bat01(i,j)<0.0)
            bat01(i,j)=-13.468;  % mask
        elseif (bat01(i,j)>=-17.503 && bat01(i,j)<-13.468)
            bat01(i,j)=-17.503;  % mask
        end
    end
end

% Pulizia lagune e dintorni...

% laguna di Grado-Marano
bat01(110:161,415:424)=0.0;
bat01(143:160,413:414)=0.0;
bat01(144:158,412)=0.0;
bat01(110:122,414)=0.0;

% laguna Venezia
bat01(1:10,340:400)=0.0;
bat01(11:41,385:400)=0.0;
bat01(11:28,379:384)=0.0;
bat01(11:15,367:380)=0.0;
for i=10:26
    for j=360:378
        if (bat01(i,j)>=-4.0 && bat01(i,j)<0.0)
            bat01(i,j)=0.0;  % mask
        end
    end
end

% ULTIMI SKAMUFFI

% punti "singolari" e canali da aprire/chiudere

% elimino punti isolati
for i=2:nx-1
    for j=2:ny-1
        if (bat01(i+1,j)==0.0 && bat01(i-1,j)==0.0...
                && bat01(i,j+1)==0.0 && bat01(i,j-1)==0.0)
            bat01(i,j)=0.0;  % mask
        end
    end
end

% Veglia
bat01(299,355)=-26.507;
bat01(300,355)=-26.507;
bat01(304,354)=-26.507;
bat01(305,354)=-26.507;
bat01(303,354)=-26.507;
bat01(303,355)=-26.507;

% zone isolate da eliminare
bat01(182:187,341:342)=0.0;
bat01(217:218,300:302)=0.0;
bat01(234:240,321:328)=0.0;
bat01(297:303,361:363)=0.0;
bat01(306:308,344:345)=0.0;
bat01(307:308,328:329)=0.0;
bat01(415:438,214:232)=0.0;
bat01(435:438,171:172)=0.0;
bat01(467:471,161:165)=0.0;
bat01(667,75:76)=0.0;
bat01(700:706,48:51)=0.0;
bat01(162,415:416)=0.0;

% baia di Panzano
bat01(169,423)=0.0;
bat01(169:170,424)=0.0;
bat01(174:175,424)=-3.001;

% FIUMI

% -2 - HPP Dub/Kup + Ombla - NEW!
bat01(750,24:29)=-3.001;
% -1 - Neretva
bat01(667:672,73)=-3.001;
%  0 - Cetina + Jadro - NEW!
bat01(572,125:130)=-3.001;
%  1 - Krka - NEW!
bat01(465,161:166)=-3.001;
%  2 - Zrmanja - NEW!
bat01(415:420,231)=-3.001;
%  3 - HPP Senj + Crikvenica - NEW!
bat01(320:325,345)=-3.001;
%  4 - Bakarac + Rjecina - NEW!
bat01(287,366:371)=-3.001;
%  5 - Rasa - NEW!
bat01(237,320:325)=-3.001;
%  6 - Mirna - NEW!
bat01(177:182,365)=-3.001;
%  7 - Dragonja - NEW!
bat01(177:182,386)=-3.001;
%  8 - Rizana - NEW!
bat01(195:200,396)=-3.001;
%  9 - Timavo
bat01((184:189)-5,421+3)=-3.001;
% 10 - Isonzo
bat01(166:171,418)=-3.001;
bat01(166,418:422)=-3.001;
% 11 - Tagliamento
bat01(113,406:411)=-3.001;
% 12 - Livenza
bat01(93-10,(404:409)-4)=-3.001;
% 13 - Piave
bat01(60+5,(391:396)+1)=-3.001;
% 14 - Sile
bat01(48-2,(387:392)-1)=-3.001;
% 15 - Brenta
bat01((6:11)+2,350-2)=-3.001;
% 16 - Adige
bat01((9:14)+1,340+4)=-3.001;
% 17 - Po
bat01(4:44,318:319)=-6.235;
bat01(25,297:319)=-6.235;
% 18 - Reno + Lamone
bat01((4:9)-1,259+12)=-3.001;
% 19 - Foglia
bat01(82:87,186)=-3.001;

% FINE FIUMI

% OBCs

% Adattamento per il nesting col modello ADRI
% IDENTICO al boundary + 2 celle di "nudging"

load('comp_bathy_y_70_mod.txt','-ascii')
z_OBCs=comp_bathy_y_70_mod;
% Attenzione al cambio di segno!
z_OBCs(z_OBCs>278.423 )=278.423;

for i=1:(nx/4-1)
    bat01(4*(i-1)+1,1:2)=-z_OBCs(i);
    bat01(4*(i-1)+2,1:2)=-z_OBCs(i);
    bat01(4*(i-1)+3,1:2)=-z_OBCs(i);
    bat01(4*(i-1)+4,1:2)=-z_OBCs(i);
end

for i=1:(nx-8)
    bat01(i,3)=bat01(i,2)-1*(bat01(i,2)-bat01(i,5))/3;
    bat01(i,4)=bat01(i,2)-2*(bat01(i,2)-bat01(i,5))/3;
end

% FINE OBCs

nxN=494;
bat02=bat01(1:nxN,129:end);
bat03=bat02;
% OBCs (bat02)

% Adattamento per il nesting col modello Copernicus
% IDENTICO al boundary + 2 celle di "nudging"
% (maxdepth -floor of bottom cell- of the grid: -218.364 m)

load('OBS_depw_214.txt','-ascii')
z_OBCs=OBS_depw_214;
% Attenzione al cambio di segno!

for i=1:floor(nxN/8)
    for j=1:8
        bat02(8*(i-1)+j,1:2)=-z_OBCs(i);
    end
end

for i=1:(nxN-8)
    bat02(i,3)=bat02(i,2)-1*(bat02(i,2)-bat02(i,5))/3;
    bat02(i,4)=bat02(i,2)-2*(bat02(i,2)-bat02(i,5))/3;
end

% aggiustamento linea di costa
bat02(181:182,3:4)=0.0;

% FINE OBCs (bat02)

% Chiusura bordi a nord e est (bat02)!
bat02(488:end,:)=0.0;
bat02(:,297:300)=0.0;

% Chiusura bordi a nord e est (bat01)!
bat01(nx,:)=0.0;
bat01(:,ny)=0.0;

% Trova la profondità max
maxdepth2=min(min(bat01));
maxdepth2_NORTH=min(min(bat02));

nomefileout1='bathy_ADRI_CADEAU_ADRIA02_rev2.bin';
fidout1=fopen(nomefileout1,'w','l');
fwrite(fidout1,bat01,'float32');
fclose(fidout1);

nomefileout2='bathy_ADRI_CADEAU_ADRIA02_N2.bin';
fidout2=fopen(nomefileout2,'w','l');
fwrite(fidout2,bat02,'float32');
fclose(fidout2);

bat01(bat01==0.0)=NaN;
figure;pcolor(bat01');colorbar;shading flat;
caxis([-270 0])
set(gcf,'renderer','zbuffer','PaperPositionMode','auto')
print -dpng -r600 bathy_ADRI_CADEAU_ADRIA02_rev2

bat02(bat02==0.0)=NaN;
figure;pcolor(bat02');colorbar;shading flat;
caxis([-270 0])
set(gcf,'renderer','zbuffer','PaperPositionMode','auto')
print -dpng -r600 bathy_ADRI_CADEAU_ADRIA02_N2

bat03(bat03==0.0)=NaN;
figure;pcolor(bat03');colorbar;shading flat;
caxis([-270 0])
set(gcf,'renderer','zbuffer','PaperPositionMode','auto')
print -dpng -r600 bathy_ADRI_CADEAU_ADRIA02_N2_noCopernicus

disp('Fatto!')
