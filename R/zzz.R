.onLoad <- function(libname, pkgname) {
    drule[["lbeta"]] <- drule[["w_lbeta"]] <- alist(a=dfun(a,b),
                                                    b=dfun(b,a))
    drule[["dfun"]] <- alist(x=dfun2(x,y),
                             y=dfun2(y,x))
}
