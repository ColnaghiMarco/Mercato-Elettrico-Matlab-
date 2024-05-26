%% Stimare il modello di regressione con dummy giornaliere e mensili
% Residuo rappresenta serie storica destagionalizzata

%% Importazione dati da file excel 2021-2024

% Importazione dati di Train
m2021 = readmatrix("Anno 2021_12.xlsx","Sheet","Prezzi-Prices","Range","O2:O8761");
m2022 = readmatrix("Anno 2022_12.xlsx","Sheet","Prezzi-Prices","Range","O2:O8761");
m2023 = readmatrix("Anno 2023_12.xlsx","Sheet","Prezzi-Prices","Range","O2:O8761");

mTrain = vertcat(m2021,m2022,m2023);

% Importazione dati di test
mTest = readmatrix("Anno 2024_04.xlsx","Sheet","Prezzi-Prices","Range","O2:O2905");

% Creazione timetable

t1 = datetime(2021, 1, 1, 1, 0, 0); % Prima ora del 1 gennaio 2021
t2 = datetime(2023, 12, 31, 24, 0, 0); % Ultima ora del 31 dicembre 2023
t_train = (t1:hours(1):t2)'; % Creazione del vettore delle date/orari

t3 = datetime(2024, 1, 1, 1, 0, 0); % Prima ora del 1 gennaio 2024
t4 = datetime(2024, 04, 30, 24, 0, 0); % Ultima ora del 30 aprile 2024
t_test = (t3:hours(1):t4)'; % Creazione del vettore delle date/orari

% mTrain = array2table(mTrain);
% mTrain = array2timetable(mTrain,"Rowtimes",t_train);
% mTest = table(mTest,t_test); %sistemare ora legale

%% Importazione dati da file excel 2016-2019

% Importazione dati di Train
m2016 = readmatrix("Anno 2016_12.xlsx","Sheet","Prezzi-Prices","Range","O2:O8785");
m2017 = readmatrix("Anno 2017_12.xlsx","Sheet","Prezzi-Prices","Range","O2:O8761");
m2018 = readmatrix("Anno 2018_12.xlsx","Sheet","Prezzi-Prices","Range","O2:O8761");

mTrain = vertcat(m2016,m2017,m2018);

% Importazione dati di test
mTest = readmatrix("Anno 2019_12.xlsx","Sheet","Prezzi-Prices","Range","N2:N2881");

% Creazione timetable

t1 = datetime(2016, 1, 1, 1, 0, 0); % Prima ora del 1 gennaio 2021
t2 = datetime(2018, 12, 31, 24, 0, 0); % Ultima ora del 31 dicembre 2023
t_train = (t1:hours(1):t2)'; % Creazione del vettore delle date/orari

t3 = datetime(2019, 1, 1, 1, 0, 0); % Prima ora del 1 gennaio 2024
t4 = datetime(2019, 04, 30, 24, 0, 0); % Ultima ora del 30 aprile 2024
t_test = (t3:hours(1):t4)'; % Creazione del vettore delle date/orari

%% Creazione Matrici per Train per fasce orarie 

% Estraggo le ore dall'array di date/orari
hours_of_day = hour(t_train);
day_of_week = day(t_train,'dayofweek')

% Creo variabili dummy per le ore del giorno
TD = dummyvar(hours_of_day + 1); % Aggiungo 1 per convertire l'intervallo [0,23] a [1,24]

%TDtt=array2table(TD)
% TDtt=array2timetable(TD,'Rowtimes',t_train);

%Creo variability dummy i giorni della settimana
TS = dummyvar(day_of_week)
TDtt = array2timetable([mTrain TS],'RowTimes',t_train);

% Viene creata  TS^*_1 (lun-dom), ..., G^*_6 (sab-dom)
%TDtt{:,1:6}=TDtt{:,2:7}-TDtt{:,1};
%TDtt = TDtt{:,1:6}

TS(:,1:6)=TS(:,2:7)-TS(:,1);
TS = TS(:,1:6);


TStt=array2timetable([mTrain TD TS],'RowTimes',t_train); % TimeTable con variabile di risposta (prezzo) e le 24 variabili dummies
TStt.Properties.VariableNames = {'Price','h1','h2','h3','h4','h5','h6','h7','h8','h9','h10','h11','h12','h13','h14','h15','h16','h17','h18','h19','h20','h21','h22','h23','h24','d1','d2','d3','d4','d5','d6'};

%% Creazione Matrici per Test per fasce orarie 

% Estraggo le ore dall'array di date/orari
hours_of_day = hour(t_test);
day_of_week = day(t_test,'dayofweek')

% Creo variabili dummy per le ore del giorno
TD = dummyvar(hours_of_day + 1); % Aggiungo 1 per convertire l'intervallo [0,23] a [1,24]
TS = dummyvar(day_of_week)

TS(:,1:6)=TS(:,2:7)-TS(:,1);
TS = TS(:,1:6);


TNtt=array2timetable([mTest TD TS], 'RowTimes',t_test);
TNtt.Properties.VariableNames = {'Price','h1','h2','h3','h4','h5','h6','h7','h8','h9','h10','h11','h12','h13','h14','h15','h16','h17','h18','h19','h20','h21','h22','h23','h24','d1','d2','d3','d4','d5','d6'};

%% Trend grafico delle fasce orarie (3.00am;7.00am;12.00am;15.00am;21.00am)
    
figure;
fasce_orarie = [3,7,12,15,21]

% Ciclo per creare i grafici
for i = 1:width(fasce_orarie) % Scorre tutte le colonne tranne la prima
    % Seleziona la variabile di ore corrente
    j=fasce_orarie(i);
    oreCorrente = TStt{:, j+1}; % Considera la colonna corrente (escludi la variabile di risposta)
    
    % Trova gli indici delle ore in cui la variabile di risposta è 1
    indici = find(oreCorrente == 1);
    
    % Se ci sono ore in cui la variabile di risposta è 1, crea un grafico
    if ~isempty(indici)
        % Crea un nuovo subplot
        subplot(5, 1, i); % 6 righe, 4 colonne, i-esimo subplot
        
        % Plot variabile di risposta rispetto alle ore correnti in cui è presente
        plot(t_train(indici), TStt{indici, 1}, '.'); % Accedi alla variabile di risposta con indici e visualizza l'orario
        
        % Aggiungi etichette
        xlabel('Day');
        ylabel('Price');
        title(['Grafico per Ora ', num2str(j)]);
        
        % Ruota etichette sull'asse x per una migliore leggibilità
        xtickangle(45);
        
        % Formatta l'asse x per visualizzare l'orario nel formato desiderato
        datetick('x', 'HH:MM');
    end
end


%% Istogramma delle medie (3.00am;7.00am;12.00am;15.00am;21.00am)

fasce_orarie = [3,7,12,15,21]
medie = zeros(size(fasce_orarie))

% Ciclo per creare istogramma
for i = 1:width(fasce_orarie) % Scorre tutte le colonne tranne la prima
    % Seleziona la variabile di ore corrente
    j=fasce_orarie(i);
    oreCorrente = TStt{:, j+1}; % Considera la colonna corrente (escludi la variabile di risposta)
    
    % Trova gli indici delle ore in cui la variabile di risposta è 1
    indici = find(oreCorrente == 1);
    
    medie(i) = mean(TStt{indici,1})
    
end

figure;

% Istogramma
bar (fasce_orarie,medie)
% Aggiunta di etichette e titoli
xlabel('Fasce Orarie');
ylabel('Media');
title('Istogramma delle Medie delle Fasce Orarie');

%% Regressione fasce orarie (3.00am;7.00am;12.00am;15.00am;21.00am)

X = TStt{:,2:end}
Y = TStt{:,1}

mdl = fitlm(X,Y)
% mdl = fitlm(TStt{:,2:end}, TStt{:,1})

%% Regressione manuale (usata per verifica dei risultati)

X = TStt{:,2:end}
Y = TStt{:,1}

n = size(Y,1)

Beta = (X'*X)^(-1)*X'*Y;  
H=X*(X'*X)^(-1)*X';
M=eye(n)-H;
yhat=H*y;
e=M*y;

%% Stima dei componenti

%trend
vT =  mdl.Coefficients.Estimate(1)+X(:,1)*mdl.Coefficients.Estimate(2);

% componente oraria
vdelt0 = mdl.Coefficients.Estimate(2:25);
vS = X(:,1:24)*vdelt0;

% componente giornaliera
vC = X(:,25:end)* mdl.Coefficients.Estimate(26:end)


%% Grafici

% Residui
subplot(2,1,1); 
plot(t_train, Y-vT-vS-vC);  
title("Residui", Interpreter="latex"); set(gca,'TickLabelInterpreter','latex');
line([t_train(1) t_train(end)], [0 0], 'LineStyle','--','Color', 'r');

% Autocorrelazione nei residui
subplot(2,1,2);   
autocorr(mdl.Residuals.Raw); title({'Autocorrelazione dei residui'}, Interpreter="latex");
set(gca,'TickLabelInterpreter','latex'); 
xlabel(' ', Interpreter='latex')
ylabel(' ')

%% Previsione primi 4 mesi 2024

[ypred, yci] = predict(mdl,TNtt{:,2:end},'Alpha',0.001);
% ypred: variabile di risposta
% yci: intervallo di confidenza

%% Rappresentazione previsione e intervallo di confidenza

times = TNtt.Properties.RowTimes;

% Traccia le previsioni
figure;
hold on;
plot(times, ypred, 'b-', 'LineWidth', 2); % Previsioni
plot(times, yci(:,1), 'r--', 'LineWidth', 1.5); % Limite inferiore dell'intervallo di confidenza
plot(times, yci(:,2), 'r--', 'LineWidth', 1.5); % Limite superiore dell'intervallo di confidenza
xlabel('Time');
ylabel('Response Variable');
title('Predictions with 99.9% Confidence Intervals');
legend('Predictions', 'Lower Confidence Bound', 'Upper Confidence Bound');
hold off;

%% Grafico comparazione Forecast e valori registrati

times = TNtt.Properties.RowTimes;

% Creazione del grafico di confronto
figure;
hold on;
plot(times, mTest, 'b-', 'LineWidth', 2); % Valori osservati
plot(times, ypred, 'r--', 'LineWidth', 2); % Previsioni
xlabel('Time');
ylabel('Response Variable');
title('Confronto tra Valori Osservati e Previsioni');
legend('Valori Osservati', 'Previsioni');
hold off;
