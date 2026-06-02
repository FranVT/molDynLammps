"""
    Script en el que se adapta el código de FORTRAN de Claudia a Julia
"""

# Usefull functions



# System information
lxn=32;
lyn=32;
lzn=32;
sigma=1;
pi2=2*pi;
ntot=5000; # Número total de partículas

# Definition of parameters

# --- Parámetros de la caja y de la rejilla en q ---
xc = lxn      # longitud de la caja en x
yc = lyn      # longitud de la caja en y
zc = lzn      # longitud de la caja en z
s = sigma     # factor de escala (sigma)
x2 = xc / 2.0 # mitad de la caja en x (no usada después, posible para centrar)
y2 = yc / 2.0 # mitad de la caja en y
z2 = zc / 2.0 # mitad de la caja en z
bin0 = 1;
qmax0 = 7;
ball=0;     # No se para que es esta variable

# parametros misteriosos
b1 = 1; # Selecciona particula de tipo 1 | central particle
b2 = 1; # Selecciona particula de tipo 2 | central particle 
b3 = 0; # Selecciona particula de tipo 3 | patch particle
b4 = 0; # Selecciona particula de tipo 4 | patch particle


# Espaciado mínimo en el espacio recíproco y ancho de bin escalado
dq0 = pi2 / xc               # Δq fundamental
bin0 = dq0 * bin0            # nuevo ancho de bin (bin0 original * dq0)
qmax = Int(floor(qmax0 / dq0)) # número entero de pasos hasta qmax0
rq0 = qmax * dq0             # valor máximo real de |q| usado
numbin = Int(floor(qmax * dq0 / bin0)) + 1  # número total de bines

# Reserva de memoria
r = zeros(ntot, 4)           # posiciones (columnas 1:3) y tipo (columna 4)
sfhis = zeros(numbin, 2)     # histograma: col1 = suma de |S(q)|^2, col2 = conteos
bc = zeros(Int, ntot)        # factor de peso (0 o 1/b1..) por partícula
kr = zeros(ntot)             # producto escalar q·r para cada partícula


# --- Lectura de posiciones y asignación de pesos ---
# Si ball==1, se usan todas las partículas con peso 1.
# Si ball==0, se usan pesos b1,b2,b3,b4 según el tipo (columna 4).
if ball == 1
    bc .= 1
    ntotav = ntot
else
    bc .= 0
    ntotav = 0
end

