## Data set 1 - Fig1 in Podani & Miklos (2002) Ecology 83
podani1 <- function(species = 42, sites = 40, labels = TRUE) {
    mat <- matrix(0, ncol = species, nrow = sites)
    diag(mat) <- 1
    mat[seq(sites + 1, prod(sites, species), by = sites + 1)] <- 3
    mat[seq((2 * sites + 1), prod(sites, species),
            by = sites + 1)] <- 1
    if(labels) {
        colnames(mat) <- paste("Spp", seq_len(species), sep = "")
        rownames(mat) <- paste("", seq_len(sites), sep = "")
    }
    mat
}

## Data set 2 - Fig 2 in Podani & Miklos (2002) Ecology 83
podani2 <- function(species = 6, sites = 40, labels = TRUE) {
    matrix(c((g <- rep(1:40, species / 2)),
             rep(rev(g), species / 2)),
           ncol = species, nrow = sites)
}

## Data set 3 - Fig9 in Podani & Miklos (2002) Ecology 83
podani3 <- function(species = 11, sites = 19,
                    abund = c(1,2,3,7,8,7,4,2,1),
                    labels = TRUE) {
    mat <- matrix(0, ncol = species, nrow = sites)
    ind <- c(outer(seq_len(length(abund)),
                   seq(0, by = sites + 1, length = species), "+"))
    mat[ind] <- abund
    if(labels) {
        colnames(mat) <- paste("Spp", seq_len(species), sep = "")
        rownames(mat) <- paste("", seq_len(sites), sep = "")
    }
    mat
}
