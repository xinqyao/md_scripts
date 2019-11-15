kb = 0.0019872041
temp = 300

## constant chsi
dgc <- -kb*temp * log((6.022*10^-4)/(8 * pi**2))

raw_comp <- read.table('PBTOT_complex.dat')[, 1]
raw_receptor <- read.table('PBTOT_receptor.dat')[, 1]
raw_ligand <- read.table('PBTOT_ligand.dat')[, 1]

## mean energy
mean_comp <- mean(raw_comp)
mean_receptor <- mean(raw_receptor)
mean_ligand <- mean(raw_ligand)

## energy variance
var_comp <- var(raw_comp)
var_receptor <- var(raw_receptor)
var_ligand <- var(raw_ligand)

## energy 3rd order moment
tri_comp <- mean((raw_comp - mean_comp)**3)
tri_receptor <- mean((raw_receptor - mean_receptor)**3)
tri_ligand <- mean((raw_ligand - mean_ligand)**3)

## entropy change up to 2nd order
s2  <- 1/(2*kb*temp)*(var_comp - var_receptor - var_ligand)

## entropy change up to 3rd order
s3  <- 1/(3*2*(kb*temp)**2)*(tri_comp - tri_receptor - tri_ligand)

## Binding free energy up to 0th order (mean)
dg0 <- mean_comp - mean_receptor - mean_ligand

## Binding free energy up to 2nd order entropy
dg2 <- dg0 + s2 + dgc

## Binding free energy up to 3rd order entropy
dg3 <- dg2 + s3

## print results
cat("            Mean (0th order)     2nd order      3rd order\n")
cat("Delta G:    ", dg0, "      ", dg2, "       ", dg3, "\n")
cat("-T*Delta S:         -         ", s2,  "       ", s3, "\n")


cis0 <- raw_comp - mean_receptor - mean_ligand + dgc
cis2 <- cis0 + s2
cis3 <- cis2 + s3

write.table(cis0, row.names=FALSE, col.names=FALSE, file='p151_cis_first.mmpbsa')
write.table(cis2, row.names=FALSE, col.names=FALSE, file='p151_cis_second.mmpbsa')
write.table(cis3, row.names=FALSE, col.names=FALSE, file='p151_cis_third.mmpbsa')


