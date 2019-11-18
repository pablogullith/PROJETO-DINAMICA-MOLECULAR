#Autores: 
#Pablo Gullith de Melo Dantas
#Gabriel Zuza Diniz
#Davi de Araújo Fernandes
from random import random
from numpy import zeros, copy, sqrt
import pylab as pl
import sys

# Define arquivo de saida
saida = open('MD2.xyz','w')

# Função que escreve a trajetória, que deve ser lida usando o programa VMD
def dump():
    saida.write('{0:d} \n'.format(N))
    saida.write('\n')
    for i in range(N):
        saida.write( '1 {0:.3f} {1:.3f} 0.000 \n'.format( x[i], y[i]) )

# Função que calcula a força entre duas partículas. O primeiro 'if' implementa o raio de cut-off. O 'if' dentro do 'else' é útil em casos onde na configuração inicial existem átomos muito próximos, e evita que átomos saiam da caixa. Recomendo remover este segundo if após ter certeza que a configuração a ser usada é estável.
def f(r):
    if r>rcut:
        return 0.0
    else:
        f = 24.0*eps*( 2.0*(1/r)**(13) -  (1.0/r)**(7))
        if (f<100000):
            return f
        else:
            return 100000.0

# Define as posições iniciais
def init_pos():
    g = 1
    for i in range(N):
        x[i] =  L/2 + g*L/2/(Nx) * (i % (Nx)) #+ L/(Nx) * (i % (Nx)) + random()*rfac*(L/Nx) 
        divy = N/Ny
        y[i] = L/2 + g*L/2/Ny * int((i / divy)) #+ L/Ny * int((i / divy)) + random()*rfac*(L/Ny) 

# Define as velocidades iniciais
def init_vel():
    for i in range(N):
         vx[i] = 2.0*v0*(random()-0.5)
         vy[i] = 2.0*v0*(random()-0.5)

# Inicializa as posições anteriores
def init_prev():
    for i in range(N):
         xprev[i] = x[i]-vx[i]*dt
         yprev[i] = y[i]-vy[i]*dt

# Programa que determina a força resultante sobre cada átomo
def calcforces():
    for j in range(N):
        if (j!=i):
            xji = x[j]-x[i]
            if (xji>0.5*L):
                xji = xji - L
            elif (xji <-0.5*L):
                xji = xji + L
                
            yji = y[j]-y[i]
            if (yji>0.5*L):
                yji = yji - L
            elif (yji <-0.5*L):
                yji = yji + L
            rji = sqrt (xji**2 + yji**2)
            fji = f(rji)
            coseno = -xji/rji
            seno = -yji/rji
            fx[i]+= fji*coseno
            fy[i]+= fji*seno

def distancia_media():
    d = 0
    for i in range(3):
        d += sqrt( (x[i]-x[i+1])**2 + (y[i]-y[i+1])**2)
    d += sqrt( (x[3]-x[1])**2 + (y[3]-y[1])**2)
    d += sqrt( (x[1]-x[3])**2 + (y[1]-y[3])**2)
    d += sqrt( (x[2]-x[0])**2 + (y[2]-y[0])**2)
    d = d/6
    return d

eps = 1.0
######################### Parâmetros da simulação ########################
# Número de átomos na simulação
N = 4
# Nx e Ny definem como será o arranjo 2D de átomos. N = Nx * Ny
Nx = 2
Ny = 2
# A velocidade inicial dos átomos é aleatória, e cada componente se situa no intervalo [-v0,v0]
v0 = 1e-4
# A variável abaixo determina quão aleatória pode ser a posição inicial. Valores altos podem gerar configurações iniciais com átomos muito próximos. 
rfac = 0.4
# A caixa é assumida quadrada, logo basta fornecer um valor de L
L = 5

# Passo da simulação (time-step)
dt = 0.005
# Número de passos
Npassos = 20000
# Raio de cut-off
rcut = 3.0
tsave = 0.1
R = 2
P = []
########################### Fim dos parâmetros  ##########################

if (Nx*Ny != N):
    print ("O arranjo escolhido de átomos não é consistente com o número total de átomos.")
    sys.exit()
    
x = zeros(N,float)
y = zeros(N,float)
vx = zeros(N,float)
vy = zeros(N,float)
xprev = zeros(N,float)
yprev = zeros(N,float)
xnew = zeros(N,float)
ynew = zeros(N, float)
fx = zeros(N,float)
fy = zeros(N,float)
T =[]

init_pos()
init_vel()
init_prev()
    
dump()
t = 0
t1 = 1.0e-6
t2 = 1.0e-4
TF = 20000*0.005
P.append(distancia_media())
T.append(t)

def Energia_Cinetica():
    K = 0
    for i in range(N):
        K += (vx[i]**2 +  vy[i]**2)*0.5
    return K

def Energia_Potencial():
    U = 0
    for i in range(N):
        for j in range(i+1, N):
            if(i!=j):
                a = x[j]-x[i]
                if (a>0.5*L):
                    a = a - L
                elif (a <-0.5*L):
                    a = a + L  
                b = y[j]-y[i]
                if (b>0.5*L):
                    b = b - L
                elif (b <-0.5*L):
                    b = b + L
                r = ((a)**2 + (b)**2)**0.5
                U += 4*((1/r**12) - (1/r**6))
    return U

while (t < TF):

    t+=dt
    for i in range(N):
        fx[i] = 0.0
        fy[i] = 0.0
        calcforces()

        if (int(t/dt)%6000 == 0):
            xprev[i] = x[i] - R*(x[i] - xprev[i])
            yprev[i] = y[i] - R*(y[i] - yprev[i])
            
        xnew[i] = 2.0*x[i] - xprev[i] + fx[i]*dt**2
        if xnew[i]<0.0:
            xnew[i]+=L
            x[i]+=L 
            xprev[i]+=L
        elif xnew[i] > L:
            xnew[i]-=L
            x[i]-=L
            xprev[i]-=L

        ynew[i] = 2.0*y[i] - yprev[i] + fy[i]*dt**2
        if ynew[i]<0.0:
            ynew[i]+=L
            y[i]+=L
            yprev[i]+=L
        elif ynew[i] > L:
            ynew[i]-=L
            y[i]-=L
            yprev[i]-=L

        vx[i] = (xnew[i]-xprev[i])/(2*dt)
        vy[i] = (ynew[i]-yprev[i])/(2*dt)

    xprev = copy(x)
    yprev = copy(y)
    x = copy(xnew)
    y = copy(ynew)
    #dump()
    P.append(distancia_media())
    T.append(t)
    if ((t+1)%tsave < t2):
        #P.append(distancia_media())
        print(t/dt/Npassos)
        #T.append(t)
        #for k in range(N):
        #    print("(%f, %f)" %(x[k], y[k]), end = " ")
        print()
        dump()

saida.close()

pl.plot(T,P)
pl.title(r"$Figura\ 2$")
pl.xlabel(r"$Tempo$")
pl.ylabel(r"$Distância$")
pl.legend()
pl.show()

"""
Analisando o gráfico obitido (Figura 2), podmos observar como inicialmente a distancia entre
as particulas é quase constante possuindo pouca variabilidade, indicando que o etádo físico 
formado por esse conjunto de particulas é o sólido.Isso ocorre até um certo ponto onde a distância
começa a variar bruscamente de maneira estocástica , nesse momento as partículas formam um líquido
Completando a simulação de derretimento de um sólido. Vale ressaltar que o sólido passa para o estado
líquido sem um 'meio termo' (esse meio termo seria ainda no estádo sólido porem com distancia entre
as particulas em uma constante diferente mais elevada), isso ocorre pois o processo de derretimento
ocorre em nossa simulação de forma mais brusca em relação a outra que apresente esse meio termo
"""