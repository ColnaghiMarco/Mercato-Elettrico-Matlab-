%% Importazione dati da tabella Excel contenente il prezzo dell'energia elettrica 

mY = readmatrix("2022_2023.xlsx", "Sheet", "Foglio1", "Range", "O2:O17521");

t1 = datetime(2022, 1, 1, 1, 0, 0); % Prima ora del 1 gennaio 2022
t2 = datetime(2023, 12, 31, 24, 0, 0); % Ultima ora del 31 dicembre 2023
t = (t1:hours(1):t2)'; % Creazione del vettore delle date/orari

%%mY.Time=t

%% Grafico Trend

y = log(mY(:,:)'); % Trasformazione logaritmica del dato

subplot(2,1,1)
plot(t, y)
title('Logaritmo')
subplot(2,1,2)
plot(t, mY(:,:)')
title('NON log')

%% Istogrammi

% Inizializzo vettori nulli dove memorizzare le medie
medie_ore = zeros(1,24);
ore_giorno = hour(t);

medie_giorni=zeros(1,7);
giorno_settimana = day(t);

medie_mesi = zeros(1, 12);
mese_anno = month(t);

for ora=1:24
    indici_ora = find(ore_giorno == ora);
    
    prezzi_ora = mY(indici_ora,:);
    
    media_ora = mean(prezzi_ora);
  
    medie_ore(ora) = media_ora;
end


%Ciclo giorni della settimana
for giorno=1:7
    indici_giorno = find(giorno_settimana == giorno);
    
    prezzi_giorno = mY(indici_giorno,:);
    
    media_giorno = mean(prezzi_giorno);
    
    medie_giorni(giorno) = media_giorno;
end

% Ciclo mesi dell'anno
for mese = 1:12
    indici_mese = find(mese_anno == mese);
    
    prezzi_mese = mY(indici_mese,:);
    
    media_mese = mean(prezzi_mese);
   
    medie_mesi(mese) = media_mese;
end



% MOstra a videp l'istogramma delle medie dei prezzi 
figure;
subplot(3,1,1)

bar(1:24, medie_ore);
xlabel('Ora del Giorno');
ylabel('Media Prezzo');
title('Media Prezzo per Ora del Giorno');
xticks(1:24); 
xticklabels({'00:00', '01:00', '02:00', '03:00', '04:00', '05:00', '06:00', '07:00', '08:00', '09:00', '10:00', '11:00', '12:00', '13:00', '14:00', '15:00', '16:00', '17:00', '18:00', '19:00', '20:00', '21:00', '22:00', '23:00'}); 

subplot(3,1,2)
bar(1:7, medie_giorni);
xticks(1:7);
xticklabels({'Lunedì', 'Martedì', 'Mercoledì', 'Giovedì', 'Venerdì', 'Sabato', 'Domenica'});
xlabel('Giorno della Settimana');
ylabel('Media Prezzo');
title('Media Prezzo per Giorno della Settimana');

subplot(3,1,3)
bar(1:12, medie_mesi);
xticks(1:12);
xticklabels({'Gennaio','Febbraio','Marzo','Aprile','Maggio','Giugno','Luglio','Agosto','Settembre','Ottobre','Novembre','Dicembre'});
xlabel('Mese');
ylabel('Media Prezzo');
title('Media Prezzo per Mese');


