#!/bin/bash
# helical example using RDataset "longley"

mkdir -p runB
cd runB

###
#
#
#
# Second using helical tools to go through in one step with the expr file
#
#
#
###

# print regressor variable (:GNP is a symbol we can use for column names)
#longley[:, :GNP]
awk -F, 'NR>1{print $3}' ../datasets/longley 

# generate y vector (Employed)
#y = longley[:, :Employed]
awk -F, 'NR>1{print $8}' ../datasets/longley > y

# to view "y" 
# cat y

# generate X matrix
#X = [ones(size(longley, 1), 1) longley[:,:GNP]]
awk -F, 'NR>1{print 1,$3}' ../datasets/longley > X

# to view "X"
# cat X

cat << EOF > lsq.expr
# This is a simple least square example expression file.
# Expression files support comments.
XpX = X' * X
XpXi = inv(XpX)
Xpy = X' * y
bhat= XpXi * Xpy
yhat = X * bhat
ehat = y - yhat
sigma2_e_hat=(ehat' * ehat)/(rows(X) - cols(X))
var_hat_beta_hat = XpXi * sigma2_e_hat

# The alternative to writing out variables with the write("file",variable) function
# is to use the --write-all flag to output all computed variables to file.
# However, many variables are not needed again, so we write out only the few we need.
write("bhat", bhat)
write("yhat", yhat)
write("ehat", ehat)
write("se_bhat", sqrt(diag(var_hat_beta_hat)) )
EOF

helical euler expr lsq.expr
echo "bhat from brute force inversion:"
helical euler print bhat


# The alternative is to use an iterative solver like PCG.
# Solves lhs * bhat = rhs
cat << EOF > lhs
MAP
X'*X
EOF

cat << EOF > rhs
MAP
X'*y
EOF

# An alternative to solving the system with the inverse is to use an iterative solver like preconditioned conjugate gradient (pcg)
# PCG also understands expressions - meaning blocks don't need to be precomputed necessarily.
helical euler pcg lhs rhs bhat
# bhat.map now contains the full solved system
echo "L2 error between bhat from pcg with diagonal preconditioning to exact bhat:"
helical euler diff bhat.precon.map bhat
# Or we can look at the single element, bhat.X_y
# helical euler print bhat.X_y

# We can also provide a custom preconditioner to PCG.
# ...in this case the perfect preconditioner consisting of the inverse of the lhs (X'*X)
cat << EOF > precon
MAP
inv(X'*X)
EOF

helical euler pcg lhs rhs bhat.precon --preconditioner precon
echo "L2 error between bhat from pcg with full inverse preconditioner to exact bhat:"
helical euler diff bhat.precon.map bhat
