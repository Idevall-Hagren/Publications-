# Adjusted Positive Agreement (APA): A Chance-Corrected Agreement Measure for Unconstrained Pile Sorts

## Publication Information

**Title:** Adjusted Positive Agreement: A Chance-Corrected Agreement Measure for Unconstrained Pile Sorts

**Author:** Jonas Idevall Hagren

**Journal:** Field Methods

**Publication Year:** 2026

**DOI:** https://doi.org/10.1177/1525822X261464653



## Use APA on your own data

Copy this into R. Your data need one row per respondent and one column per item, with each cell giving the pile that respondent put the item in. Pile labels can be numbers or text and only need to be consistent within a row.

```r
# Co-placement vector
coplacement <- function(piles) {
  same <- outer(piles, piles, "==")
  as.integer(same[upper.tri(same)])
}

# APA between two co-placement vectors
APA <- function(A_i, A_j) {
  observed <- sum(A_i * A_j)
  s_i <- sum(A_i)
  s_j <- sum(A_j)
  n_pairs <- length(A_i)

  expected <- (s_i * s_j) / n_pairs
  maximum <- min(s_i, s_j)
  denominator <- maximum - expected

  if (!is.finite(denominator) || denominator <= 0) {
    return(NA_real_)
  }

  value <- (observed - expected) / denominator
  max(0, min(1, value))
}

# Respondent-by-respondent APA matrix
APA_matrix <- function(sorts) {
  A <- lapply(seq_len(nrow(sorts)), function(r) {
    coplacement(as.character(unlist(sorts[r, ])))
  })
  n <- length(A)
  M <- diag(1, n)
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      M[i, j] <- M[j, i] <- APA(A[[i]], A[[j]])
    }
  }
  dimnames(M) <- list(rownames(sorts), rownames(sorts))
  M
}
```

### Example

This made-up example has 25 respondents sorting 16 animals. Run the functions above first, then read in the data:

```r
sorts <- read.csv(row.names = 1, text = "
id,dog,wolf,cat,lion,cow,horse,whale,seal,eagle,sparrow,duck,swan,salmon,shark,bee,butterfly
id1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2
id2,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8
id3,1,1,1,1,1,1,1,1,2,2,2,2,3,3,4,4
id4,1,2,1,2,1,1,3,3,2,2,1,3,3,3,2,2
id5,1,1,1,1,2,2,3,3,4,4,4,4,5,5,6,6
id6,1,1,2,2,3,3,4,4,5,5,5,5,6,6,7,7
id7,1,2,1,2,3,3,4,4,2,5,3,4,4,4,5,5
id8,1,1,1,1,1,1,1,1,2,2,3,3,4,4,5,5
id9,1,1,1,1,1,1,1,1,2,2,3,3,4,4,2,2
id10,1,1,1,1,2,2,3,3,4,4,5,5,6,6,7,7
id11,1,2,1,2,3,3,4,4,2,2,3,4,4,4,2,2
id12,1,1,1,1,1,1,2,2,3,3,3,3,2,2,4,4
id13,1,1,2,2,3,3,4,4,5,5,6,6,7,7,7,7
id14,1,2,1,2,3,3,4,4,2,5,3,6,6,4,5,5
id15,1,1,1,1,1,1,1,1,2,2,2,2,3,3,3,3
id16,1,2,1,2,3,3,4,4,5,5,5,5,6,6,6,6
id17,1,1,1,1,1,1,2,2,3,3,3,3,4,4,5,5
id18,1,2,1,2,1,1,3,3,2,4,1,3,3,3,4,4
id19,1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6
id20,1,1,1,1,2,2,3,3,4,4,4,4,4,4,4,4
id21,1,2,1,2,3,3,4,4,2,2,3,5,5,4,2,2
id22,1,1,1,1,1,1,2,2,3,3,4,4,5,5,6,6
id23,1,1,1,1,2,2,3,3,4,4,5,5,3,3,6,6
id24,1,2,1,2,1,1,3,3,2,4,1,5,5,3,4,4
id25,1,1,1,1,1,1,2,2,3,3,3,3,3,3,3,3
")
```

Each row is one respondent, and each number is the pile that respondent put the animal in. For example, id2 made eight piles: dog and wolf in pile 1, cat and lion in pile 2, cow and horse in pile 3, and so on. Fourteen of the respondents sort the animals by kind (mammals, birds, fish, insects). Seven (id4, id7, id11, id14, id18, id21, id24) sort them by where they live, such as at home, on the farm, in the wild or in water. Four (id9, id12, id16, id23) mix the two.

Compute the APA matrix and summarise it:

```r
M <- APA_matrix(sorts)
round(M, 2)

# Average APA across all pairs of respondents
M_off <- M
diag(M_off) <- NA   # leave out each respondent's agreement with themselves
mean(M_off, na.rm = TRUE)

# Per respondent: number of piles, sorting density (number of item pairs
# placed together) and average APA with the other respondents
respondents <- data.frame(
  piles = apply(sorts, 1, function(x) length(unique(x))),
  density = apply(sorts, 1, function(x) sum(coplacement(x))),
  mean_APA = rowMeans(M_off, na.rm = TRUE)
)
round(respondents, 2)
```

Average APA is 0.56. Sorting density ranges from 8 for id2 (eight piles of two) to 56 for id1 (two piles of eight). The seven habitat sorters have the lowest average APA (0.32 to 0.42). Those who sort by kind range from 0.52 to 0.81.

Then decompose the matrix:

```r
ev <- eigen(M, symmetric = TRUE)

ev$vectors <- apply(ev$vectors, 2, function(v) if (sum(v) < 0) -v else v)

# Eigenvalues
eigenvalues <- data.frame(rank = seq_along(ev$values), eigenvalue = ev$values)
round(eigenvalues, 3)

# Respondents' coordinates
coordinates <- ev$vectors
dimnames(coordinates) <- list(rownames(M), paste0("EV", seq_along(ev$values)))
round(coordinates, 3)

# Scree plot
plot(ev$values, type = "b", pch = 19, xlab = "Eigenvalue rank", ylab = "Eigenvalue")
abline(h = 0, lty = 2)
```

On the first eigenvector, all respondents have positive coordinates, with the seven respondents sorting by habitat being closest to zero (0.08 to 0.10, against 0.20 to 0.28 for everyone else). On the second eigenvector, the habitat sorters (0.31 to 0.39) are clearly separated from everyone else (−0.15 to 0.10).


## License

Code: MIT License (see `LICENSE`). Article: CC BY 4.0.

## Contact

Jonas Idevall Hagren, Department of Global Public Health, Karolinska Institutet.
