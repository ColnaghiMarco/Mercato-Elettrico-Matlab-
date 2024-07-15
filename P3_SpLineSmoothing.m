%% Analisi Prezzo mercato elettrico zona NORD Italia
% Analisi delle fasce orarie 3.00-7.00-12.00-15.00-21.00
% 0-Caricamento dati 2021-2023
% 1-Estrapolazione elle 5 serie storiche
% 2-Eliminazione del rumore con Wavelet e creazione serie storiche
% con i residui del wavelet
% 3-Modello ARIMA (7,0,0) sulle serie sotriche dei residui Kernel

%% 0-Caricamento dei dati

% Importazione dati di Train
m2021 = readtable("Anno 2021_12.xlsx","Sheet","Prezzi-Prices",VariableNamingRule='preserve');
m2022 = readtable("Anno 2022_12.xlsx","Sheet","Prezzi-Prices",VariableNamingRule='preserve');
m2023 = readtable("Anno 2023_12.xlsx","Sheet","Prezzi-Prices",VariableNamingRule='preserve');

dati = vertcat(m2021,m2022,m2023);

n = size(dati,1); % n.ro delle osservazioni
y = dati.NORD; % selezione della serie dei prezzi elettrici
nDays = n/24;
   
matrice_prezzi=reshape(y,[24,nDays])';
yh24 = array2table(matrice_prezzi);

nomivar = "h"+(1:24);
yh24.Properties.VariableNames = nomivar;

t1 = datetime(2021, 1, 1, 1, 0, 0); % Prima ora del 1 gennaio 2021 - data iniziale
n = size(yh24,1); % n.ro delle osservazioni
t = t1 + caldays(0:n-1)';

%% 1- Creazione delle 5 serie storice

yAverage = zeros(1,24);

for i = 1:24
    yAverage(i) = mean(yh24{i,:});
end

% Istogramma
bar (1:24,yAverage);
% Aggiunta di etichette e titoli
xlabel('Fasce Orarie');
ylabel('Media');
title('Istogramma delle Medie delle Fasce Orarie');

y3 = yh24.h3; % Da sistemare poi

y_Ore = array2table([yh24.h3,yh24.h7,yh24.h12,yh24.h15,yh24.h21],"VariableNames",{'h3','h7','h12','h15','h21'});

%% 2 - Smoothing con Spline

% Creazione della tabella per i residui
residues = table('Size', size(y_Ore), 'VariableTypes', repmat({'double'}, 1, width(y_Ore)), 'VariableNames', y_Ore.Properties.VariableNames);
t_array = days(t - t(1));

% Ciclo per ogni colonna di y_Ore
for col = 1:width(y_Ore)
    spline_data = y_Ore{:, col}; % Estrai i dati della colonna corrente
    
    % Applica il smoothing con spline (spaps)
    p = 0.9999; % Parametro di smoothing, puoi regolarlo tra 0 e 1
    [sp, smoothed_data] = spaps(t_array, spline_data, p);

    residues{:, col} = spline_data - smoothed_data';
    
    % Grafici
    figure;
    subplot(4,1,1);
    plot(t, spline_data);
    title(['Serie Temporale Originale con Rumore (', y_Ore.Properties.VariableNames{col}, ')']);
    xlabel('Tempo');
    ylabel('Valori');
    
    subplot(4,1,2);
    plot(t, smoothed_data);
    title(['Serie Temporale Denoised con Smoothing Spline (', y_Ore.Properties.VariableNames{col}, ')']);
    xlabel('Tempo');
    ylabel('Valori');
    
    subplot(4,1,3);
    plot(t, [spline_data, smoothed_data']);
    title(['Confronto Serie e Approssimazione Spline (', y_Ore.Properties.VariableNames{col}, ')']);
    xlabel('Tempo');
    ylabel('Valori');

    subplot(4,1,4);
    plot(t, residues{:, col});
    title(['Residui (', y_Ore.Properties.VariableNames{col}, ')']);
    xlabel('Tempo');
    ylabel('Valori');
end

%% 3.1- Modello regARIMA(7;0;0) - Verifica parametri p,d,q

%% Verifico che la serie residui sia stazionaria (d=0)

% Creazione di una tabella per memorizzare i risultati del test ADF
adf_results = table('Size', [width(residues), 5], ...
                    'VariableTypes', {'string', 'double', 'double', 'double', 'double'}, ...
                    'VariableNames', {'Variable', 'HypothesisRejected', 'pValue', 'TestStatistic', 'CriticalValue'});

% Esegui il test di Dickey-Fuller aumentato per ogni colonna
for col = 1:width(residues)
    [h, pValue, stat, cValue, reg] = adftest(residues{:, col});
    
    % Memorizza i risultati nella tabella
    adf_results.Variable(col) = residues.Properties.VariableNames{col};
    adf_results.HypothesisRejected(col) = h;
    adf_results.pValue(col) = pValue;
    adf_results.TestStatistic(col) = stat;
    adf_results.CriticalValue(col) = cValue;
end

% Visualizza i risultati
disp('ADF Test Results:');
disp(adf_results);

%% Verifica del parametro p=7

figure;
for col = 1:width(residues)
    [pacf, lags] = parcorr(residues{:, col});
    
    % Plot
    subplot(width(residues), 2, 2*col-1);
    plot(residues{:, col});
    title(['Residui serie storica - ', residues.Properties.VariableNames{col}]);
    
    subplot(width(residues), 2, 2*col);
    parcorr(residues{:, col}, 28); % dal grafico imposterei p=3
    title(['PACF - ', residues.Properties.VariableNames{col}]);
end

%% Verifica parametro q=0

figure;
for col = 1:width(residues)
    % Plot
    subplot(width(residues), 2, 2*col-1);
    plot(residues{:, col});
    title(['Residui serie storica - ', residues.Properties.VariableNames{col}]);
    
    % Calcola e plotta la funzione di autocorrelazione (ACF)
    subplot(width(residues), 2, 2*col);
    autocorr(residues{:, col}, 28); % dal grafico imposterei q=7;
    title(['Autocorrelazione (ACF) - ', residues.Properties.VariableNames{col}]);
end


%% 3.2- Modello regARIMA(7;0;0) - Stima del modello 
% Devo includere i 7 elementi precedenti al 1 gennaio alla tab predittori dato che applico un modello con ritardo 7
tX = t1 - 7;
tX = tX + caldays(0:n-1+7)';

mD = dummyvar(month(tX));       % 12 dummy mensili
dD = dummyvar(weekday(tX));     % 7 dummy per giorno della settimana

% Tabella coi predittori
X_dummy = array2table([mD, dD]);
X_dummy = table2array(X_dummy);

% Creazione di una tabella per memorizzare i valori stimati
fitted_values = array2table(zeros(height(residues), width(residues)), ...
                            'VariableNames', residues.Properties.VariableNames);

% Ciclo per ogni colonna di residues
for col = 1:width(residues)
    Y_Arima = residues{:, col};
    
    Mdl = regARIMA('ARLags', 7, 'D', 0, 'MALags', [], 'Intercept', 0);
    estimatedModel = estimate(Mdl, Y_Arima, 'X', X_dummy);
    
    % Visualizzare i risultati
    disp(['Estimated Model for ', residues.Properties.VariableNames{col}, ':']);
    disp(estimatedModel);
    
    % Applicare il modello ARIMA stimato ai predittori X
    Y_fit = simulate(estimatedModel, length(Y_Arima), 'X', X_dummy);
    
    % Salvare i valori stimati nella tabella
    fitted_values{:, col} = Y_fit;
    
    % Plot dei risultati
    figure;
    t = 1:length(Y_Arima);
    plot(t, Y_Arima, 'b', 'LineWidth', 1.5); % Serie temporale originale
    hold on;
    plot(t, Y_fit, 'r--', 'LineWidth', 1.5); % Serie ottenuta dal modello ARIMA con predittori
    hold off;
    
    legend('Original Series', 'ARIMA Model Fitted with Exogenous Predictors');
    title(['ARIMA(7,0,0) Model Fitted with Exogenous Predictors - ', residues.Properties.VariableNames{col}]);
    xlabel('Time');
    ylabel('Value');
end

% Visualizzare i valori stimati
disp('Fitted Values:');
disp(fitted_values);

%% 3.3 - Valutazione modello regARIMA

% Creazione di una tabella per memorizzare i risultati
evaluation_results = table('Size', [width(residues), 6], ...
                           'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double'}, ...
                           'VariableNames', {'Variable', 'ADF_HypothesisRejected', 'ADF_pValue', 'ADF_TestStatistic', ...
                                             'MSE', 'MAE'});

% Creazione di una singola figura per tutti i grafici
figure;

% Ciclo per ogni colonna di residues
for col = 1:width(residues)
    Y_Arima = residues{:, col};
    Y_fit = fitted_values{:, col};
    
    % Analisi dei residui
    residuals_ARIMA = Y_Arima - Y_fit;

    % Plot dell'autocorrelazione dei residui
    subplot(3, 2, col);  % Creazione del subplot
    autocorr(residuals_ARIMA);
    title(['Autocorrelazione dei Residui - ', residues.Properties.VariableNames{col}]);

    % Test di normalità dei residui
    [~, pValue, ~] = adtest(residuals_ARIMA);
    disp(['Test di normalità dei residui per ', residues.Properties.VariableNames{col}, ' - p-Value: ', num2str(pValue)]);
    
    % Test di eteroschedasticità dei residui
    [h, pValue, stat] = lbqtest(residuals_ARIMA);
    disp(['Test di eteroschedasticità dei residui per ', residues.Properties.VariableNames{col}, ' - p-Value: ', num2str(pValue)]);

    % Calcolo di metriche di performance
    mse = mean(residuals_ARIMA.^2);
    mae = mean(abs(residuals_ARIMA));
    disp(['Errore Quadratico Medio (MSE) per ', residues.Properties.VariableNames{col}, ': ', num2str(mse)]);
    disp(['Errore Assoluto Medio (MAE) per ', residues.Properties.VariableNames{col}, ': ', num2str(mae)]);

    % Calcolo del R quadro
    RSS = sum(residuals_ARIMA.^2);
    TSS = sum((Y_fit - mean(Y_fit)).^2);
    R2 = 1 - (RSS / TSS);
    disp(['Coefficiente di Determinazione (R^2) per ', residues.Properties.VariableNames{col}, ': ', num2str(R2)]);
    
    % Salvataggio dei risultati nella tabella
    evaluation_results.Variable(col) = residues.Properties.VariableNames{col};
    evaluation_results.ADF_HypothesisRejected(col) = h;
    evaluation_results.ADF_pValue(col) = pValue;
    evaluation_results.ADF_TestStatistic(col) = stat;
    evaluation_results.MSE(col) = mse;
    evaluation_results.MAE(col) = mae;
end

% Visualizzare i risultati complessivi
disp('Risultati della Valutazione del Modello:');
disp(evaluation_results);


%% 3.4 - Test di Ljung-Box

% Creazione di una tabella per memorizzare i risultati del test di Ljung-Box
lbq_results = table('Size', [width(residues), 4], ...
                    'VariableTypes', {'string', 'logical', 'double', 'double'}, ...
                    'VariableNames', {'Variable', 'LBQ_HypothesisRejected', 'LBQ_pValue', 'LBQ_TestStatistic'});

% Ciclo per ogni colonna di residues per il test di Ljung-Box
for col = 1:width(residues)
    Y_Arima = residues{:, col};
    Y_fit = fitted_values{:, col};
    
    % Analisi dei residui
    residuals_ARIMA = Y_Arima - Y_fit;

    % Test di Ljung-Box per verificare l'autocorrelazione residua
    lag = 10; % Imposta il numero di ritardi per il test di Ljung-Box
    [h_lbq, pValue_lbq, stat_lbq] = lbqtest(residuals_ARIMA, 'Lags', lag);
    disp(['Test di Ljung-Box per ', residues.Properties.VariableNames{col}, ' - p-Value: ', num2str(pValue_lbq)]);

    % Salvataggio dei risultati nella tabella
    lbq_results.Variable(col) = residues.Properties.VariableNames{col};
    lbq_results.LBQ_HypothesisRejected(col) = h_lbq;
    lbq_results.LBQ_pValue(col) = pValue_lbq;
    lbq_results.LBQ_TestStatistic(col) = stat_lbq;
end

% Visualizzare i risultati del test di Ljung-Box
disp('Risultati del Test di Ljung-Box:');
disp(lbq_results);