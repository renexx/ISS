# Created by Octave 4.0.0, Sat Dec 15 18:49:39 2018 CET René Bolf xbolfr00@stud.fit.vutbr.cz
#####################################################################################################################################################################
printf('\nTOTO JE 1. PRIKLAD.\n');
[sn, fs] = audioread('xbolfr00.wav');
sn = sn'; % potrebujeme radkovy vektor
N = length(sn);
t = N/fs;
pocet_bin_sym = N/16;
printf('\nVZORKOVACIA FREKVENCIA SINGNALU JE : %f [Hz].\n', fs); %16000
printf('DLZKA SIGNALU VO VZORKOCH : %f, \n V SEKUNDACH JE : %f [s].\n', N, t); %32000  , 2000s
printf('POCET REPREZENTOVANYCH BINARNYCH SYMBOLOV  : %f, \n ', pocet_bin_sym); %20000
#####################################################################################################################################################################
printf('\nTOTO JE 2. PRIKLAD.\n');
x = 8;
for i = 1:2000                                                                                                                                                                                                                            
   if sn(x) > 0
        dekod(i) = 1;
   else sn(x) < 0
        dekod(i) = 0;
   endif
   x = x + 16;
endfor  
figure('2');
hold on;
plot(sn(1:320),'b');
stem([8:16:320],dekod(1:20),'r');
xlabel('t [s]');
ylabel('sn[n],symbols');
set(gca, 'XTick', 0:32:320);
set(gca, 'XTickLabel', 0:0.002:0.02);
set(gca, 'YTick', -1:0.2:1);
grid;
title("Uloha 2");
legend("audio signal","dekod symbol")
hold off;
print('-dpng', 'picture/dvojka.png');
#####################################################################################################################################################################
printf('\nTOTO JE 3. PRIKLAD.\n');
b = [0.0192 -0.0185 -0.0185 0.0192];
a = [1.0000 -2.8870 2.7997 -0.9113];
p = roots(a); 
if (isempty(p) | abs(p) < 1) 
  printf('filter je stabilni\n');
else
  printf('NESTABILNI !!!\n');
endif
figure('3');
zplane(b, a);
xlabel('Realna cast');
ylabel('Imaginarna cast');
title("Uloha 3");
legend("zplane","poly","nuly");
print('-dpng', 'picture/3.png');
#####################################################################################################################################################################
printf('\nTOTO JE 4. PRIKLAD.\n');
figure('4');
f1 =(0:255) / 256 * fs / 2;
f1
[H1,H2] = freqz(b,a,256);
[minval,column] = min(min(H1(1:32),[],2));
printf("mezna je %f",f1(column));
plot(f1, abs(H1));
xlabel('f [Hz]');
ylabel('|H(f)|');
title("Uloha 4");
legend('Modul kmitoctovej charakteristiky');
grid;
print('-dpng', 'picture/4.png');
printf('\nObrazok je v priecinku picture.\n');
printf('\nMALA BY TO BYT DOLNA PRIEPUST.\n');
#####################################################################################################################################################################
#####################################################################################################################################################################
printf('\nTOTO JE 5. PRIKLAD.\n');
filtered_signal = (filter(b, a, sn));
figure('5');
hold on;
plot(sn(1:320),'b');
stem([8:16:320],dekod(1:20),'m');
plot(filtered_signal(1:320),'r');
xlabel('t [s]');
ylabel('sn[n],ssn[n],ss_symbols');
set(gca, 'XTick', 0:32:320);
set(gca, 'XTickLabel', 0:0.002:0.02);
set(gca, 'YTick', -1:0.2:1);
grid;
title("Uloha 5");
hold off;
legend("audio signal","dekod symbol","filtered signal");
print('-dpng', 'picture/khhkkh.png');
#####################################################################################################################################################################
printf('\nTOTO JE 6. PRIKLAD.\n');
filtered_signal = (filter(b, a, sn));
signalshift = shift(filtered_signal,-17);
xsh = 8;
for ish = 1:2000                                                                                                                                              #fclose(FILE);                                                                                
   if signalshift(xsh) > 0
        dekodsh(ish) = 1;
   else signalshift(xsh) < 0
        dekodsh(ish) = 0;
   endif
   xsh = xsh + 16;
endfor  
figure('6');
hold on;
plot(sn(1:320),'b');
stem([8:16:320],dekodsh(1:20),'m');
plot(filtered_signal(1:320),'r');
plot(signalshift(1:320),'g');
xlabel('t [s]');
ylabel('sn[n],ss_{shifted}[n],symbols');
set(gca, 'XTick', 0:32:320);
set(gca, 'XTickLabel', 0:0.002:0.02);
set(gca, 'YTick', -1:0.2:1);
grid;
title("Uloha 6");
legend("audio signal","dekod sym shift sn","filtered signal","signal shifted");
hold off;
print('-dpng', 'picture/6.png');
#####################################################################################################################################################################
fprintf('\nTOTO JE 7. PRIKLAD.\n');
pocet_chybnych = 0;
for k = 1:2000
  pocet_chybnych = pocet_chybnych + xor(dekod(k), dekodsh(k)) ;
endfor 
vysledok = pocet_chybnych / 2000 * 100;
printf("pocet chyb je %f\n",pocet_chybnych);
printf("chyba je v percentach %f\n",vysledok);
#####################################################################################################################################################################
fprintf('\nTOTO JE 8. PRIKLAD.\n');                                            
signalorig = abs(fft(sn));
signalfiltered = abs(fft(filtered_signal));
signalorig = signalorig (1:N/2);
signalfiltered = signalfiltered (1:N/2);
freq =(0:N/2-1)/N * fs; 
figure("8");
hold on;
plot(freq,signalorig,"b");
plot(freq,signalfiltered,"r");
xlabel('f [Hz]');
grid;
hold off;
legend("signal original"," filtered signal")
print('-dpng', 'picture/8.png');
################################################################################################################################################################
printf('\nTOTO JE 9. PRIKLAD.\n');  
xmin = min(min(sn));
xmax = max(max(sn));
value = 50;
x = linspace(xmin,xmax,value);
vzdialenostbodov = abs(xmin - xmax);
probability = hist(sn,x) / N / vzdialenostbodov;
figure("9");
plot(x,probability);
title("Uloha 9");
xlabel("s[n]");
ylabel("Pravdepodobnost");
legend("funckia hustoty rozdelenia pravdepodobnosti");
set(gca, 'XTick', -1:0.2:1);
grid;
integralp = sum(probability * vzdialenostbodov);
print('-dpng', 'picture/9.png');
printf("Integrál: %f\n",integralp);
#####################################################################################################################################################################
fprintf('\nTOTO JE 10. PRIKLAD.\n');
k = (-50 : 50);
R = xcorr(sn,"biased");
R = R(k + N)
figure(10);
plot(k, R);
xlabel('k');
ylabel('R[k]');
legend('Autokorelacny koeficient');
title("Uloha 10");
grid;
print('-dpng', 'picture/10.png');
#####################################################################################################################################################################
% 11) Napište hodnotu koeficientu R[0],R[1],R[16]
printf('\nTOTO JE 11. PRIKLAD\n');
printf('Hodnota koeficientu R[0] je %f.\n', R(51));
printf('Hodnota koeficientu R[1] je %f.\n', R(52));
printf('Hodnota koeficientu R[16] je %f.\n', R(67));
#####################################################################################################################################################################
printf('\nTOTO JE 12. PRIKLAD)\n');
L = 50;
value12 = 50;
min12 = min(sn);
max12 = max(sn);
x = linspace(min12,max12, value12);
h = zeros(L, L);
xcol = x(:);
ycol = sn(:)';
bigy = repmat(ycol, L, 1); 
bigx = repmat(xcol, 1, N);
[dummy, ind1] = min(abs(bigy - bigx));
ind2 = ind1(1 + 1 : N);
for i = 1 : N - 1,
	d1 = ind1(i);
	d2 = ind2(i);
	h(d1, d2) = h(d1, d2) + 1;
end
surf = (x(2) - x(1)) ^ 2;
p = h / N / surf;
figure("12");
imagesc(x, x, p);
axis xy;
colorbar;
xlabel('x2');
ylabel('x1');
title("Uloha 12");
printf(" Vygeneroval sa obazok\n");
print('-dpng', 'picture/11.png');
#####################################################################################################################################################################
printf('\nPRIKLAD 13.)\n');
check = sum(sum (p)) * surf;
printf('Overenie, ze sa jedna o spravnu zdruzenu funkciu hustoty pravdepodobnosti teda, ze vysledok je  %d.\n',check);
#####################################################################################################################################################################
printf('\nTOTO JE 14. PRIKLAD)\n');
% make col vector out of x and clone it L times. 
x = x(:); X1 = repmat(x,1,L);
% make row vector out of x and clone it L times. 
x = x'; X2 = repmat(x,L,1); 
% now put it together, don't forget to multipl by tile surface
r = sum(sum (X1 .* X2 .* p)) * surf;
printf('Hodnota koeficientu R[1] je %f.\n',r);
