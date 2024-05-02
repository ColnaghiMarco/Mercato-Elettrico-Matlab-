%% Analisi ARMA Electric price, Italian Market  NORTH ZONE

%% Caricamento dati
y = readmatrix("2022_2023.xlsx", "Sheet", "Foglio1", "Range", "O2:O17521");

t1 = datetime(2022, 1, 1, 0, 0, 0); % Prima ora del 1 gennaio 2022
t2 = datetime(2023, 12, 31, 23, 0, 0); % Ultima ora del 31 dicembre 2023
t = (t1:hours(1):t2)'; % Creazione del vettore delle date/orari

mY = table(y, t)
mY.Properties.VariableNames = {'Price', 'Date'} % Rinomina le colonne

all_price = mY.Price
all_data = mY.Date

[n,~] = size (y)
n2 = n*0.20 % ?da capire la dimensione per la finsetra della media mobile

price2023=y(n2:n,1)

%% Media mobile (metodo iterativo=troppo dispendioso)
hm = zeros(n,1)
hm = all_price

for i=1:n2
    hm(i+n2) = mean(all_price(n2:n2+i))
end

MSE = mean((hm(n2+1:n)-all_price(n2+1:n)).^2);

%% Media mobile funz movmean
windows = 24*30*12
hm = movmean(y,windows)

subplot(2,1,1)
plot(t,y)
xlabel('Tempo')
ylabel('Price')
title('Price value')

subplot(2,1,2)
plot(t,hm)
xlabel('Tempo')
ylabel('Media Mobile')
title('Media mobile con finestra mobile')

%% Autocorrelazione vs Par. Autocorrelazione
figure
subplot(2,1,1)
autocorr(log(all_price))
title("Autocorr")
subplot(2,1,2)
parcorr(log(all_price))
title("Part. Autocorr")

%% Calcolo AIC e BIC 
cpmax=5
cqmax=2

mAIC = NaN(cpmax+1,cqmax+1)
mBIC= NaN(cpmax+1,cqmax+1)

for p=0:cpmax
    for q=0:cqmax
        Mdl = arima(p,0,q)
        EstMdl = estimate(Mdl,all_price)
        mAIC(p+1,q+1) = EstMdl.summarize.AIC
        mBIC(p+1,q+1) = EstMdl.summarize.BIC
    end
end
%%
%Selezione con AIC
[dminAIC,iInd] = min(mAIC(:))
[ir, ic] = ind2sub(size(mAIC),iInd)
cpsel = ir-1
cqsel = ic-1

% Selezione mediante BIC 
[dminBIC iInd] = min(mBIC(:))
[ir, ic] = ind2sub(size(mBIC),iInd)
cpsel = ir-1
cqsel = ic-1
% stima del modello scelto
Mdl = arima(cpsel,0,cqsel)
Mdl.Constant = mean(y) %applico la media come costante del modello
EstMdl = estimate(Mdl, y)
summarize(EstMdl)

%% Calcolo RMS
residual = infer(EstMdl, y)
RMSE = sqrt(sum(residual.^2))/length(residual)

Previsione = all_price+residual

%% Forecasting 2023
[Final_Forecast,ym] = forecast(EstMdl, n2,'Y0', price2023);
lower_prac = Final_Forecast - 1.96*sqrt(ym);
upper_prac = Final_Forecast + 1.96*sqrt(ym);

figure;
h1 = plot(all_data,price2023,'Color','g', 'LineWidth',1.5);
hold on;
h2 = plot(all_data(n2+1:n+n2),Final_Prac_Forecast,'LineWidth',2, 'Color','r');
h3 = plot(all_data(n2+1:n+n2),all_price(n2+1:n+n2));
h4 = plot(all_data(n2+1:n+n2),lower_prac,'r:','LineWidth',2);
h5 = plot(all_data(n2+1:n+n2),upper_prac,'r:','LineWidth',2);
h6 = plot(all_data(n2+1:n+n2),hm(n2+1:n+n2),'LineWidth',2, 'Color','m', 'LineStyle','--');
legend([h1,h2,h3,h4, h6],"Original Data","ARIMA","Original","95% Confidence Interval", "Historical Mean","Location","Best");
grid on;
hold off;
