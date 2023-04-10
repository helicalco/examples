#!/bin/bash
# helical example using RDataset "longley"

mkdir -pv runA
cd runA


###
#
#
#
# First using helical tools to go through step by step
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

# Define a helper function to run `helical euler`
he() {
    helical euler "$@"
}
export -f he

# X'X Matrix
#XpX = X'X
# he is now a shortcut for "helical euler" tools
he expr "X' * X" -o XpX

# (X'X)^1 (inverse of X'X)
#XpXi = inv(XpX)
he expr "inv(XpX)" -o XpXi

# X'y
#Xpy = X'y
he expr "X' * y" -o Xpy

# solve for b-hat
#bhat = XpXi * Xpy
he expr "XpXi * Xpy" -o bhat

cat bhat 
#or
he print bhat

# y-hat = X * bhat
#yhat = X * bhat
he expr "X * bhat" -o yhat

# e-hat (y - y_hat) because y = y_hat + e
#ehat = y .- yhat
he expr "y - yhat" -o ehat

# model residual variance (2nd part is the 'degrees of freedom')
#sigma2_e_hat = sum(ehat.^2) / (size(X, 1) - size(X, 2))
he expr "(ehat' * ehat)/(rows(X) - cols(X))" -o sigma2_e_hat

# sqrt that variance
#sqrt(sigma2_e_hat)
he expr "sqrt(sigma2_e_hat)" -o sqrt_sigma2_e_hat

# (co)variance matrix of coefficients (variance in estimates)
#var_hat_beta_hat = XpXi * sigma2_e_hat
he expr "XpXi * sigma2_e_hat" -o var_hat_beta_hat

# equivalent to vcov() of the model later...

# calculate SE of coef estimates
#sqrt.(diag(var_hat_beta_hat))
he expr "sqrt(diag(var_hat_beta_hat))" -o se_bhat

# An alternative to solving the system with the inverse is to use an iterative solver like preconditioned conjugate gradient (pcg)
he pcg XpX Xpy bhat
# bhat.map now contains the full solved system
he print bhat.map
# Or we can look at the single element, bhat.Xpy
# he print bhat.Xpy

echo "L2 error between bhat from pcg to exact bhat:"
he diff bhat.map bhat
