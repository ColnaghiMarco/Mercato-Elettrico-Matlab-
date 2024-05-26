%% Stima con il polinomio kernel
% residui dovrebbero rappresentare serie de-trendizzata

%% Importazione dati da file excel 2021-2024

% Importazione dati di Train
m2021 = readmatrix("Anno 2021_12.xlsx","Sheet","Prezzi-Prices","Range","O2:O8761");
m2022 = readmatrix("Anno 2022_12.xlsx","Sheet","Prezzi-Prices","Range","O2:O8761");
m2023 = readmatrix("Anno 2023_12.xlsx","Sheet","Prezzi-Prices","Range","O2:O8761");

y = vertcat(m2021,m2022,m2023);
n = length(y);

t1 = datetime(2021, 1, 1, 1, 0, 0); % Prima ora del 1 gennaio 2021
t2 = datetime(2023, 12, 31, 24, 0, 0); % Ultima ora del 31 dicembre 2023
t = (t1:hours(1):t2)'; % Creazione del vettore delle date/orari

%% grafico della serie
% g = figure('Name','Serie Temperature Globali');
plot(t, y, LineWidth=2);     xlim([t1, t2]);
set(gca,'TickLabelInterpreter','latex');
xlabel('t', Interpreter='latex')
ylabel('Gradi centigradi','Interpreter','latex',Rotation=90)
grid on; box on;
% exportgraphics(g,'gTemperature.pdf')

%% plot kernel 

ch = 20;
vh = (-ch:ch)';

vw_g = exp(((vh).^2)/(-2*ch)) ;
vw_g = vw_g/(sum(vw_g)); % kernel gaussiano

vw_t = (1-(abs(vh)/ch).^3).^3; 
vw_t = vw_t/(sum(vw_t));  % kernel tricube

subplot(1,2,1);
bar(vh,vw_g,0.5,'r')
title('Gaussian ','Interpreter','latex');    
ylim([0 0.15]); 
set(gca,'TickLabelInterpreter','latex');
subplot(1,2,2);
bar(vh, vw_t, 0.5, 'r'); 
title('Tricube ','Interpreter','latex');    
ylim([0 0.15]); 
set(gca,'TickLabelInterpreter','latex');

%% Estensione della serie con h previsioni 
% mediante trend lineare adattato alle prime e ultime m  osservazioni
chmax = 40;   
m = 30; 
% Modello stimato in base alle prime cm osservazioni
y0   = flip(y(1:m)); 
X0 = (1:m)'; 
mdl0  = fitlm(X0,y0);

yhatprima  = predict(mdl0, (m+1:m+chmax)'); 
yhatprima  = flip(yhatprima);

y1   = y(n-m+1:n); 
X1 = (1:m)'; 
mdl1  = fitlm(X1,y1);
yhatdopo  = predict(mdl1, (m+1:m+chmax)'); 

vyext = [yhatprima; y; yhatdopo];   % serie estesa 
%plot(vyext)

%% Kernel Gaussiano
chmax = 25;
chmin = 3; 
h=0;

Yfitted=zeros(n,chmax-chmin+1);
vCVscore = NaN(chmax-chmin,1);
for h = chmin:chmax
    vyext = [yhatprima(end-h+1:end); y; yhatdopo(1:h)];
    vh = (-h:h)';

    % kernel gaussiano
    vw_g = exp(((vh).^2)/(-2*ch)) ;
    vw_g = vw_g/(sum(vw_g));
   
    dw0  = vw_g(h+1);
    vyg = filter(vw_g, 1, vyext );
    vyhat = vyg(2*h+1:end);
    Yfitted(:,h-chmin+1) = vyhat;
    vCVscore(h-chmin+1) = sum((y-vyhat).^2)/((1-dw0)^2);
end

[min_val_Gauss, idx_min_Gauss] = min(vCVscore);
h_opt = chmin + idx_min_Gauss - 1;
vyext = [yhatprima(end-h+1:end); y; yhatdopo(1:h)];
vh = (-h:h)';
vw_g = exp(((vh).^2)/(-2*ch)) ;
vw_g = vw_g/(sum(vw_g));
dw0  = vw_g(h+1);
vyg = filter(vw_g, 1, vyext );
vyhat_gs = vyg(2*h+1:end);
plot(t, [y, vyhat_gs]); 

%% Kernel Tricube
chmax = 25;
chmin = 3; 
h=0;

Yfitted=zeros(n,chmax-chmin+1);
vCVscore = NaN(chmax-chmin,1);
for h = chmin:chmax
    vyext = [yhatprima(end-h+1:end); y; yhatdopo(1:h)];
    vh = (-h:h)';

    % kernel tricube
    vw_t = (1-(abs(vh)/h).^3).^3; 
        vw_t = vw_t/(sum(vw_t));
    dw0  = vw_t(h+1);
    vyt = filter(vw_t, 1, vyext );
    vyhat = vyt(2*h+1:end);
    Yfitted(:,h-chmin+1) = vyhat;
    vCVscore(h-chmin+1) = sum((y-vyhat).^2)/((1-dw0)^2);
end

[min_val_tricube, idx_min_tricube] = min(vCVscore);
h_opt = chmin + idx_min_tricube - 1;
vyext = [yhatprima(end-h+1:end); y; yhatdopo(1:h)];
vh = (-h:h)';
vw_t = (1-(abs(vh)/h).^3).^3; 
vw_t = vw_t/(sum(vw_t));
dw0  = vw_t(h+1);
vyt = filter(vw_t, 1, vyext );
vyhat_tc = vyt(2*h+1:end);
plot(t, [y, vyhat_tc]);

%% Confronto tra i due kernel gaussiano-tricube

subplot(2,1,1);
plot(t, [y, vyhat_gs, vyhat_tc]);

subplot(2,1,2);
b = [min_val_Gauss min_val_tricube]
bar(b)
set(gca, 'XTickLabel', {'Gaussiano', 'Tricube'});
xlabel('Tipologia Kernel');
ylabel('Minimo');
title('Grafico a Barre dei Valori Minimi');

%% Grafico dei residui best kernel 3 anni
residui = y - vyhat_tc;

% Plotta i dati originali e i valori adattati
figure;
subplot(2,1,1);
plot(t, y, 'b', t, vyhat_tc, 'r');
legend({'Dati Originali', 'Dati Adattati (Tricube)'});
xlabel('Tempo');
ylabel('Valori');
title('Dati Originali e Adattati');

% Plotta i residui
subplot(2,1,2);
scatter(t, residui, '.k');
xlabel('Tempo');
ylabel('Residui');
title('Grafico dei Residui');

%% Grafico dei residui best kernel 3 mesi

t1 = datetime(2021, 1, 1, 1, 0, 0); 
t2 = datetime(2021, 03, 31, 23, 0, 0); 
tm = (t1:hours(1):t2)'; % Creazione del vettore delle date/orari

y_tm = y(1:length(tm));
vyhat_tc_tm = vyhat_tc(1:length(tm));
residui_tm = residui(1:length(tm));

figure;
subplot(2,1,1);
plot(tm, y_tm, 'b', tm, vyhat_tc_tm, 'r');
legend({'Dati Originali', 'Dati Adattati (Tricube)'});
xlabel('Tempo (mesi)');
ylabel('Valori');
title('Dati Originali e Adattati su 3 mesi');

% Plotta i residui su tm
subplot(2,1,2);
plot(tm, residui_tm, '.k');
xlabel('Tempo (mesi)');
ylabel('Residui');
title('Grafico dei Residui su 3 mesi');

%% Grafico dei residui best kernel 3 settimane

t1 = datetime(2021, 1, 1, 1, 0, 0); 
t2 = datetime(2021, 01, 21, 23, 0, 0); 
ts = (t1:hours(1):t2)'; % Creazione del vettore delle date/orari

y_ts = y(1:length(ts));
vyhat_tc_ts = vyhat_tc(1:length(ts));
residui_ts = residui(1:length(ts));

figure;
subplot(2,1,1);
plot(ts, y_ts, 'b', ts, vyhat_tc_ts, 'r');
legend({'Dati Originali', 'Dati Adattati (Tricube)'});
xlabel('Tempo (mesi)');
ylabel('Valori');
title('Dati Originali e Adattati su 3 sett.');

% Plotta i residui su ts
subplot(2,1,2);
plot(ts, residui_ts, '.k');
xlabel('Tempo (mesi)');
ylabel('Residui');
title('Grafico dei Residui su 3 sett.');