Modification/update of parm99.dat (Hornak & Simmerling)
MASS
DN  14.01         0.530               sp2 nitrogen in amide groups

BOND
C -DN  490.0    1.335
CT-DN  337.0    1.449
CX-DN  337.0    1.449

ANGL
CX-C -DN    70.0      116.60    AA general  (was CT-C-N)
H1-CT-DN    50.0      109.50    AA general  changed based on NMA nmodes
H1-CX-DN    50.0      109.50    AA general  (was H1-CT-N)
CT-CT-DN    80.0      109.70    AA ala, general    (JACS 94, 2657)
CT-CX-DN    80.0      109.70    AA ala, general    (was CT-CT-N)
C -DN-CT    50.0      121.90    AA general
C -DN-CX    50.0      121.90    AA general  (was C-N-CT)
CT-DN-CX    50.0      118.00    AA pro             (DETAR JACS 99,1232) (was CT-N -CT)
C -CX-DN    63.0      110.10    AA general  (was C-CT-N)
DN-C -O     80.0      122.90    AA general

DIHE
X -C -DN-X    4   33.33        180.0             2.         AA,NMA
X -CT-DN-X    6    0.00          0.0             2.         JCC,7,(1986),230
X -CX-DN-X    6    0.00          0.0             2.         JCC,7,(1986),230 (was X -CT-N -X )
C -DN-CX-C    1    0.00          0.0            -4.         four amplitudes and
C -DN-CX-C    1    0.42          0.0            -3.         phases for phi
C -DN-CX-C    1    0.27          0.0            -2.
C -DN-CX-C    1    0.00          0.0             1.
N -CX-C -DN   1    0.00          0.0            -4.         four amplitudes and
N -CX-C -DN   1    0.55        180.0            -3.         phases for psi
N -CX-C -DN   1    1.58        180.0            -2.  
N -CX-C -DN   1    0.45        180.0             1.
DN-CX-C -N    1    0.00          0.0            -4.         four amplitudes and
DN-CX-C -N    1    0.55        180.0            -3.         phases for psi
DN-CX-C -N    1    1.58        180.0            -2.
DN-CX-C -N    1    0.45        180.0             1.
CT-CT-DN-C    1    0.00          0.0            -4.         four amplitudes and
CT-CT-DN-C    1    0.40          0.0            -3.         phases for phi'
CT-CT-DN-C    1    2.00          0.0            -2.
CT-CT-DN-C    1    2.00          0.0             1.
CT-CX-DN-C    1    0.000         0.0            -4.
CT-CX-DN-C    1    0.800         0.0            -3.
CT-CX-DN-C    1    1.800         0.0            -2.
CT-CX-DN-C    1    2.000         0.0             1.
DN-C -CX-C8   1    0.000         0.0            -4.    psi'
DN-C -CX-C8   1    0.400         0.0            -3.
DN-C -CX-C8   1    0.200         0.0            -2.
DN-C -CX-C8   1    0.200         0.0             1.
DN-C -CX-CT   1    0.000         0.0            -4.
DN-C -CX-CT   1    0.400         0.0            -3.
DN-C -CX-CT   1    0.200         0.0            -2.
DN-C -CX-CT   1    0.200         0.0             1.
DN-C -CX-2C   1    0.000         0.0            -4.
DN-C -CX-2C   1    0.400         0.0            -3.
DN-C -CX-2C   1    0.200         0.0            -2.
DN-C -CX-2C   1    0.200         0.0             1.
DN-C -CX-3C   1    0.000         0.0            -4.
DN-C -CX-3C   1    0.400         0.0            -3.
DN-C -CX-3C   1    0.200         0.0            -2.
DN-C -CX-3C   1    0.200         0.0             1.

IMPR
X -CT-DN-CT         1.0          180.          2.           JCC,7,(1986),230
X -CT-DN-CX         1.0          180.          2.           JCC,7,(1986),230 (was X -CT-N -CT)

NONB
DN          1.8240  0.1700
