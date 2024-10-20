%% 1-Caricamento dati di train

% Importazione dati di Train
m2021 = readtable("Anno 2021_12.xlsx","Sheet","Prezzi-Prices",VariableNamingRule='preserve');
m2022 = readtable("Anno 2022_12.xlsx","Sheet","Prezzi-Prices",VariableNamingRule='preserve');
m2023 = readtable("Anno 2023_12.xlsx","Sheet","Prezzi-Prices",VariableNamingRule='preserve');

dati = vertcat(m2021,m2022,m2023);

n = size(dati,1); % n.ro delle osservazioni
y = dati.PUN; % selezione della serie dei prezzi elettrici
nDays = n/24;
   
matrice_prezzi=reshape(y,[24,nDays])';
yh24 = array2table(matrice_prezzi);

nomivar = "h"+(1:24);
yh24.Properties.VariableNames = nomivar;

t1 = datetime(2021, 1, 1, 1, 0, 0); % Prima ora del 1 gennaio 2021 - data iniziale
n = size(yh24,1); % n.ro delle osservazioni
t = t1 + caldays(0:n-1)';

y_PUN = yh24.h7;

% Creazione della sequenza temporale con ritardo
tX = t1 - 7;
tX = tX + caldays(0:n-1+7)';

% Creazione delle variabili dummy
mD = dummyvar(month(tX));       % 12 variabili dummy per i mesi
dD = dummyvar(weekday(tX));     % 7 variabili dummy per i giorni della settimana

% Imposto dicembre come combinazione lineare degli altri 11 mesi
december_indices = month(tX) == 12;  % Indici delle righe relative a dicembre
mD(december_indices, 1:11) = mD(december_indices, 1:11) - 1/11;
mD(:, 12) = [];  % Rimuove dicembre per avere solo 11 colonne

% Imposto Domenica come combinazione lineare degli altri giorni
sunday_indices = weekday(tX) == 7;  % Indici delle righe relative alla domenica
dD(sunday_indices, 1:6) = dD(sunday_indices, 1:6) - 1/6;
dD(:, 7) = [];  % Rimuove la domenica per avere solo 6 colonne

% Creazione della tabella dei predittori
X_dummy = array2table([mD, dD]);
X_dummy = table2array(X_dummy);

%% 1.5 - Caricamento dati di previsione

% Importazione dati di Test
mPrev = readtable("Anno 2024.xlsx","Sheet","Prezzi-Prices",VariableNamingRule='preserve');

np = size(mPrev,1); % n.ro delle osservazioni
yp = mPrev.NORD; % selezione della serie dei prezzi elettrici
nDaysP = np/24;
   
matrice_prezzi_previsione=reshape(yp,[24,nDaysP])';
yh24P = array2table(matrice_prezzi_previsione);

nomivar = "h"+(1:24);
yh24P.Properties.VariableNames = nomivar;

tp1 = datetime(2024, 1, 1, 1, 0, 0); % Prima ora del 1 gennaio 2021 - data iniziale
np = size(yh24P,1); % n.ro delle osservazioni
tp = tp1 + caldays(0:np-1)';

y_PUN_2024 = yh24P.h7;

% Creazione della sequenza temporale con ritardo
tpX = tp1 - 7;
tpX = tpX + caldays(0:np-1+7)';

% Creazione delle variabili dummy
mD = dummyvar(month(tpX));       % 12 variabili dummy per i mesi
dD = dummyvar(weekday(tpX));     % 7 variabili dummy per i giorni della settimana

% Imposto dicembre come combinazione lineare degli altri 11 mesi
december_indices = month(tpX) == 12;  % Indici delle righe relative a dicembre
mD(december_indices, 1:11) = mD(december_indices, 1:11) - 1/11;
mD(:, 12) = [];  % Rimuove dicembre per avere solo 11 colonne

% Imposto Domenica come combinazione lineare degli altri giorni
sunday_indices = weekday(tpX) == 1;  % Indici delle righe relative alla domenica
dD(sunday_indices, 1:6) = dD(sunday_indices, 1:6) - 1/6;
dD(:, 7) = [];  % Rimuove la domenica per avere solo 6 colonne

% Creazione della tabella dei predittori
X_dummy_Prev = array2table([mD, dD]);
X_dummy_Prev = table2array(X_dummy_Prev);

%% 2-Applicazione dello smoother (Wavelet)

% Definizione dei parametri delle wavelet
wname = 'db4'; % Tipo di wavelet (Daubechies 4)
level = 5;     % Livello di decomposizione
t_array = days(t - t(1));

% Inizializza la tabella per i residui e per i valori smussati (trend_data)
residues = zeros(size(y_PUN));
trend_data = zeros(size(y_PUN));

spline_data = y_PUN; % Estrai i dati 
    
% Applica il smoothing con spline (spaps)
p = 0.9999; % Parametro di smoothing, puoi regolarlo tra 0 e 1
[sp, smoothed_data_col] = spaps(t_array, spline_data, p);
smoothed_data = smoothed_data_col'; % Traspongo l'array dei dati smoothed

residues = spline_data - smoothed_data;

subplot(2,1,1);
plot(t, y_PUN, 'b-'); % Linea blu continua
hold on;
plot(t, smoothed_data, 'r--'); % Linea rossa tratteggiata
hold off;
title('Serie originale e smussata');
xlabel('Tempo');
ylabel('Valori');
legend('Serie originale', 'Serie smussata'); % Aggiungi legenda se desiderato

subplot(2,1,2);
plot(t, residues, 'k-'); % Residui in nero continuo
title('Residui');
xlabel('Tempo');
ylabel('Valori');


%% 3 - Stima stagionalità modello regARIMA(7;0;7) - Stima dei parametri

% Valutazione del parametro d tramite test ADF

% Creazione di una tabella per memorizzare i risultati del test ADF
adf_results = table('Size', [1, 5], ... % Cambia la dimensione per una singola variabile
                    'VariableTypes', {'string', 'double', 'double', 'double', 'double'}, ...
                    'VariableNames', {'Variable', 'HypothesisRejected', 'pValue', 'TestStatistic', 'CriticalValue'});

% Esegui il test ADF
[h, pValue, stat, cValue] = adftest(residues);

% Memorizza i risultati nella tabella
adf_results.Variable = {'Residues'}; % Usa un nome generico
adf_results.HypothesisRejected = h;
adf_results.pValue = pValue;
adf_results.TestStatistic = stat;
adf_results.CriticalValue = cValue;

% Visualizza i risultati
disp('ADF Test Results:');
disp(adf_results);

% Calcola la PACF
[pacf, lags] = parcorr(residues); % calcola PACF con 28 lag

% Plot dei residui
figure;
subplot(2, 1, 1); % Poiché hai solo una colonna, usa solo due sotto-grafici
plot(residues);
title('Residui serie storica');

% Plot della PACF
subplot(2, 1, 2);
stem(lags, pacf); % Utilizza stem per visualizzare correttamente la PACF
title('PACF - Residui');
xlabel('Lag');
ylabel('PACF');

% Opzionale: aggiungi una griglia per migliorare la leggibilità
grid on;

% Plot dei residui
figure;
subplot(2, 1, 1); % Solo due sotto-grafici
plot(residues);
title('Residui serie storica');

% Calcola e plotta la funzione di autocorrelazione (ACF)
subplot(2, 1, 2);
autocorr(residues, 28); % calcola ACF con 28 lag
title('Autocorrelazione (ACF) - Residui');

% Opzionale: aggiungi una griglia per migliorare la leggibilità
grid on;

%% 3.5 - Stima stagionalità modello regARIMA(7;0;7) - Stima del modello

% Creazione di una tabella per memorizzare i valori stimati
fitted_values = zeros(size(residues));

Y_Arima = residues;

Mdl = regARIMA('ARLags', 1:7, 'D', 0,'MALags' , 1:7);
estimatedModel = estimate(Mdl, Y_Arima, 'X', X_dummy);

 % Salva il modello stimato
fitModel = estimatedModel;
% Applicare il modello ARIMA stimato ai predittori X
Y_fit = simulate(estimatedModel, length(Y_Arima), 'X', X_dummy);

% Salvare i valori stimati nella tabella
fitted_values= Y_fit;
    
subplot(3,1,1);
plot(t,residues);
title('residui serie');
subplot(3,1,2);
plot(t,fitted_values);
title('stima modello regARIMA')
subplot(3,1,3);
plot(t,residues);
hold on;
plot(t,fitted_values);
hold off;
title('Comparazione residui-stima')

%% 4 - Previsone modello additivo

n_predictions = np;
window_size = nDays; 
X0 = X_dummy;
XF = X_dummy_Prev;

% Creazione di una tabella per memorizzare i valori previsti con ARIMA
rolling_predicted_values = array2table(zeros(n_predictions, 1));
residues_forecast = zeros(n_predictions,1);

% Creo array dove salvare le future previsioni
Y_forecast_2024 = array2table(zeros(n_predictions, 1));

% stima ARIMA primo valore
estimateModel = fitModel;

%dati storici
Y_dati_storici = yh24{:,7};  %Cambiare in base all'ora in esame
Y_residui_storici = residues;
smoothed_forecast = trend_data;

t_array = days(t - t(1));

store_smoothed = zeros(n_predictions,1);

for i=1:n_predictions
    %Previsione i-esimo elemento
    Y_pred = forecast(estimateModel, 1, 'Y0', Y_residui_storici, 'X0', X0, 'XF', XF);
    rolling_predictions(i) = Y_pred;
    %Modello additivo trend+stagionalità
    Y_forecast_2024{i,1} = Y_pred + smoothed_forecast(nDays);

    % Store dei valori di smoothing
    store_smoothed (i) = smoothed_forecast(nDays);

    % Avanzo di 1 nella rolling windows
    %Y_dati_storici = [Y_dati_storici(i+1:end,:); Y_forecast_2024{1:i,:}];   
    Y_dati_storici = [Y_dati_storici(2:end,:); yh24P.h7(i,:)];
   
    % Applica il smoothing con spline (spaps)
    [sp, smoothed_data_col] = spaps(t_array, Y_dati_storici, p);
    smoothed_forecast = smoothed_data_col'; % Traspongo l'array dei dati smoothed

     % Aggiorno anche le variabili esogene X0 e XF
    X0 = [X0(2:end, :); XF(1,:)]; % Aggiorna X0 con la nuova variabile esogena
    XF = XF(2:end, :); % Aggiorna XF per i passi successivi

    residues_forecast(i, 1) = Y_dati_storici(end, 1) - smoothed_forecast(end, 1);
    Y_residui_storici = [Y_residui_storici(2:end,:);residues_forecast(i,:)];
end

% Grafico previsioni-serie
figure;
plot(tp, yh24P.h7);              % Primo grafico
hold on;
plot(tp, Y_forecast_2024{:,1});  % Secondo grafico
hold off;
ylim([0 inf]);                   % Imposta l'asse delle ordinate per partire da 0
title('Comparazione serie e previsione')

%% 5 - Valutazione delle previsioni tramite metriche MAE, RMSE, MAPE

% Calcola le metriche di valutazione
actual = yh24P.h7; % Valori reali
predicted = Y_forecast_2024{:,1}; % Valori previsti

% Calcolo delle metriche
MAE = mean(abs(actual - predicted)); % Mean Absolute Error
RMSE = sqrt(mean((actual - predicted).^2)); % Root Mean Squared Error
MAPE = mean(abs((actual - predicted) ./ actual)) * 100; % Mean Absolute Percentage Error

% Visualizza i risultati
fprintf('Metriche di valutazione delle previsioni:\n');
fprintf('MAE: %.4f\n', MAE);
fprintf('RMSE: %.4f\n', RMSE);
fprintf('MAPE: %.4f%%\n', MAPE);

