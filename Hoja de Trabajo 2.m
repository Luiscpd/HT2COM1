cd 'COM'
M = dlmread('Hoja No 1 - Entrada.txt', '\t',1,0)
size(M)
t = M(:,1);
V = M(:,2);
N = length(t)
N = length(V)
Ts = (max(t)-min(t))/(N-1)
N = length(t)
Ts = (max(t)-min(t))/(N-1)
Fs = 1/Ts
delta = t(2)-t(1);
EsConstante = true(N-1,1);
for i = 2:N
  if abs((t(i)-t(i-1)-delta)/delta)*100 > 0.01
  EsConstante(i-1) = false;
  end
end
q = sign(V); 
q(q==0) = 1;
adq = abs(diff(q));
adq(242) = 0;
Np2 = sum(adq>0);
Np = Np2/2;
fc = Np/(max(t)-min(t))
Tc = 1/fc
s = sin(2*pi*fc*t);
pkg load signal
[R Ntau] = xcorr(s,V);
[Rmax i] = max(R);
tau = Ntau(i)*delta
P1A = mean(V.^2)
P1B=mean(V(1:end-1).^2)
P1C = mean(V(1:end-2).^2)
c = cos(2*pi*fc*t);
A = V'*c/(c'*c);
P2 = A^2/2;

##### Grafica #####

subplot(3,1,1); p1 = plot(t,V);
xlabel('t (s)');
ylabel('V(t)')

##### Grafica #####

subplot(3,1,2); p2 = plot(Ntau*delta,R);
title('Correlacion Cruzada');
xlabel('\tau (s)');
ylabel('R(\tau)')
grid;


###### Pregunta 5 #####

gi = A*c;
Nt = Tc/Ts;
Po = 1/Nt * (gi'*gi);
error = ((P2-Po)/Po)*100;


##### Pregunta 6 #####

subplot(3,1,3); p3 = plot(t, V, 'b', t, adq, 'r');
xlabel('t(s)');
legend('V(t)','adq(t)')
grid on
grid minor
