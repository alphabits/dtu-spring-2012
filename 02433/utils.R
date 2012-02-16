statdist = function (gamma) {
    m = dim(gamma)[1]
    one.row = rep(1, m)
    U = matrix(rep(one.row, m), nrow=m)
    Id = diag(one.row)
    inv = solve(Id - gamma + U)
    return(one.row %*% inv)
}
