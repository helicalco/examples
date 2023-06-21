#!/bin/bash
# Example bash script using helical toolset to compute the partitioned genetic variance
# accounting for the proportion of genetic variance explained by markers
# Required input files (just need to be available on disk in the current folder) are:
# 1. "G0" genetic variance-covariance matrix 
# 2. "m" marker matrix with a row for every animal and a column for every marker
# Output files are
# 1. "invG0" inverse of G0
# 2. "Gepsilon" genetic variance-covariance explained by markers for non-genotyped animals
# 3. "invGepsilon" inverse of Gepsilon
# 4. "Ga" genetic variance-covariance matrix part not explained by markers
# 5. "invGa" inverse of Ga
# 6. "Galpha" genetic variance-covariance matrix attributed to each fitted locus
# 7. "invGalpha" inverse of Galpha

C=0.45
PI=0.0
PQ=0.25

cat << EOF > ss-shm-vc.expr
//#C                   # This is the proportion of genetic variance explained by markers
//#pi                  # The proportion of markers with zero effect - often 0 meaning all markers useful
//#num=cols(m)                 # The number of loci
c=${C}
pi=${PI}
pq=${PQ}

multiplier = c/(2*pq*cols(m)*(1-pi))

invG0=inv(G0)
write("invG0",invG0)
//# The genetic variance explained by markers for non-genotyped animals (ie residual polygenic effect RPE)
Gepsilon    = c*G0
invGepsilon = inv(Gepsilon)
write("Gepsilon",Gepsilon)
write("invGepsilon",invGepsilon)

//# The genetic variance not explained by markers
Ga          = (1.0-c)*G0
invGa       = inv(Ga)
write("Ga",Ga)
write("invGa",invGa)

//# The genetic variance attributed to each fitted locus c/(2*0.25*pq*num*(1-p))
Galpha      = multiplier*G0
invGalpha   = inv(Galpha)

write("Galpha",Galpha)
write("invGalpha",invGalpha)
EOF

helical euler expr ss-shm-vc.expr
