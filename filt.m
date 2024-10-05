% Parametry filtru
Fs = 96000;          % Częstotliwość próbkowania (w Hz)
Fpass = 20000;       % Częstotliwość przepustu (w Hz)
Fstop = 25000;       % Częstotliwość zatrzymania (w Hz)
Rp = 1;              % Tłumienie w paśmie przepustu (w dB)
Rs = 60;             % Tłumienie w paśmie zaporowym (w dB)

% Projektowanie filtru
[N, Wn] = cheb1ord(Fpass/(Fs/2), Fstop/(Fs/2), Rp, Rs); % Obliczenie parametrów filtru
[b, a] = cheby1(N, Rp, Wn); % Projektowanie filtru Chebysheva

% Zapisywanie współczynników do pliku
save('filter.mat', 'b', 'a');

% Jeśli chcesz sprawdzić charakterystykę filtru, możesz użyć:
fvtool(b, a);