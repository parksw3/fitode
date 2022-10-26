.onLoad <- function(libname, pkgname) {
    k <- x <- a <- b <- y <- NULL ## defeat visible-binding checks
    drule[["NBconst"]] <- alist(k=dfun(k,x),
                                x=dfun(x,k)+1/x)

    drule[["lbeta"]] <- alist(a=dfun(a,b), b=dfun(b,a))

    drule[["dfun"]] <- alist(x=dfun2(x,y),
                             y=dfun2(y,x))

}
