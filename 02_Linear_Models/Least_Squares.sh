# helical example using RDataset "longley"

mkdir run
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
cat y

# generate X matrix
#X = [ones(size(longley, 1), 1) longley[:,:GNP]]
awk -F, 'NR>1{print 1,$3}' ../datasets/longley > X

# to view "X"
cat X


# X'X Matrix
#XpX = X'X
helical euler expr "X' * X" -o XpX

# (X'X)^1 (inverse of X'X)
#XpXi = inv(XpX)
helical euler expr "inv(XpX)" -o XpXi

# X'y
#Xpy = X'y
helical euler expr "X' * y" -o Xpy

# solve for b-hat
#bhat = XpXi * Xpy
helical euler expr "XpXi * Xpy" -o bhat

cat bhat 
#or
mprint -p bhat 

# y-hat = X * bhat
#yhat = X * bhat
helical euler expr "X * bhat" -o yhat

# e-hat (y - y_hat) because y = y_hat + e
#ehat = y .- yhat
helical euler expr "y - yhat" -o ehat

# model residual variance (2nd part is the 'degrees of freedom')
#sigma2_e_hat = sum(ehat.^2) / (size(X, 1) - size(X, 2))
sizeX1=$(wc X | awk '{print $1}')
sizeX2=$(awk '{print NF}' X | head -1)
helical euler expr "(ehat' * ehat)/($sizeX1 - $sizeX2)" -o sigma2_e_hat

# sqrt that variance
#sqrt(sigma2_e_hat)
helical euler expr "sqrt(sigma2_e_hat)" -o sqrt_sigma2_e_hat

# (co)variance matrix of coefficients (variance in estimates)
#var_hat_beta_hat = XpXi * sigma2_e_hat
helical euler expr "XpXi * sigma2_e_hat" -o var_hat_beta_hat

# equivalent to vcov() of the model later...

# calculate SE of coef estimates
#sqrt.(diag(var_hat_beta_hat))
helical euler expr "sqrt(diag(var_hat_beta_hat))" -o se_bhat



cd ../
mkdir runB
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
cat y

# generate X matrix
#X = [ones(size(longley, 1), 1) longley[:,:GNP]]
awk -F, 'NR>1{print 1,$3}' ../datasets/longley > X

# to view "X"
cat X

sizeX1=$(wc X | awk '{print $1}')
sizeX2=$(awk '{print NF}' X | head -1)

cat << EOF > expr-1
XpX = X' * X
XpXi = inv(XpX)
Xpy = X' * y
bhat= XpXi * Xpy
yhat = X * bhat
ehat = y - yhat
sigma2_e_hat=(ehat' * ehat)/($sizeX1 - $sizeX2)
var_hat_beta_hat = XpXi * sigma2_e_hat

write("bhat", bhat)
write("yhat", yhat)
write("ehat", ehat)
write("se_bhat", sqrt(diag(var_hat_beta_hat)) )

EOF
helical euler expr expr-1
