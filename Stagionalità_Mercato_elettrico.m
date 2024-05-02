%% Estrazione dati
mY = readmatrix("2022_2023.xlsx", "Sheet", "Foglio1", "Range", "O2:O17521")';
y = log(mY(:,:)'); % Trasformazione logaritmica del dato

t1 = datetime(2022, 1, 1, 0, 0, 0); % Prima ora del 1 gennaio 2022
t2 = datetime(2023, 12, 31, 23, 0, 0); % Ultima ora del 31 dicembre 2023
t = (t1:hours(1):t2)'; % Creazione del vettore delle date/orari

n=length(y);

%creazione estrattore TR
TR = timerange(t1-30,t2)

%% creiamo le variabili di calendario fino al 2050
t1c = datetime(1900, 1, 1);
t2c = datetime(2050,12,31);
vdate = (t1c:hours(1):t2c)';     % giorni di calendario

vdow = day(vdate, 'dayofweek');  % giorni della settimana (1=domenica)
TD = dummyvar(vdow);

TDtt=array2timetable(TD,'Rowtimes',vdate);
Gtt = retime(TDtt,'monthly','sum'); 

GBtt=Gtt(TR,:);

%% Effetti mensili

% Variabili dummies mensili (effetti stagionali)
tempi=GBtt.Properties.RowTimes;
Seas=dummyvar(tempi.Month);

%parte da eliminare nel programma completo
mese_anno = month(t)*100 + year(t);
mesi_unici = unique(mese_anno);

prezzo_medio_mese = zeros(size(mesi_unici));

% Loop sui mesi unici
for i = 1:length(mesi_unici)
    indici_mese = find(mese_anno == mesi_unici(i));
    prezzi_mese = y(indici_mese);
    media_mese = mean(prezzi_mese);
    prezzo_medio_mese(i) = media_mese;
end


% Fine parte da eliminare

trend=(1:24)';
TStt=array2timetable([prezzo_medio_mese trend Seas],'RowTimes',tempi);

%% Unisco effetti giornalieri + mensili

yXtt=[TStt GBtt];
yXtt.Properties.VariableNames(1:14)=["Prezzo Medio Mese" "trend" "Seas"+(1:12)];

%% MOdello di previsione

yX=timetable2table(yXtt,'ConvertRowTimes',false);

% utilizziamo fino al 31 dic 2022
ff=(12);
yXt=yX(1:ff,:);
t=tempi(1:ff);
%%

X=yXt{:,2:end};
y=yXt{:,1};
% I dati da gen 2023 fino a dic 2023 sono utilizzati per il test
yXtnew=yX(ff+1:ff+12,:);
tnew=tempi(ff+1:ff+12);
ynew=yX{ff+1:ff+12,1};

%% Stima del modello
mdl=fitlm(yXt,'ResponseVar','y');
disp(mdl)



