library(testthat)

test_that("no CDs are added if add_cds is FALSE", {
    d <- data.frame(NAME = "A", SMILES = "CC", RT = 1)
    out <- preprocess_data(d, add_cds = FALSE, verbose = 0)
    expect_equal(out, d)
})

test_that("polynomial features are added correctly", {
    M <- data.frame(NAME = c("A","B"), SMILES = c("CC","CCC"), RT = c(1,2))
    X <- data.frame(Fsp3 = c(1,2))
    d <- cbind(M, X)
    out <- preprocess_data(d, add_cds = FALSE, degree_polynomial = 3, verbose = 0)
    expect_true(all(c("Fsp3","Fsp3^2","Fsp3^3") %in% colnames(out)))
    expect_equal(out$`Fsp3^2`, out$Fsp3^2)
})

test_that("interaction terms are added correctly", {
    M <- data.frame(NAME = c("A","B"), SMILES = c("CC","CCC"), RT = c(1,2))
    X <- data.frame(Fsp3 = c(1,2), apol = c(3,4))
    d <- cbind(M, X)
    out <- preprocess_data(d, add_cds = FALSE, interaction_terms = TRUE, verbose = 0)
    expect_true("Fsp3*apol" %in% colnames(out))
    expect_equal(out$`Fsp3*apol`, X$Fsp3 * X$apol)
})

test_that("NA removal works", {
    M <- data.frame(NAME = c("A","B"), SMILES = c("CC","CCC"), RT = c(1,2))
    X <- data.frame(Fsp3 = c(1, NA), apol = c(3,4))
    d <- cbind(M, X)
    out <- preprocess_data(d, add_cds = FALSE, rm_na = TRUE, verbose = 0)
    expect_true("Fsp3" %notin% colnames(out))
    expect_true("apol" %in% colnames(out))
})

test_that("near-zero variance removal works", {
    M <- data.frame(NAME = c("A","B","C"), SMILES = c("CC","CCC","CCCC"), RT = c(1,2,3))
    X <- data.frame(Fsp3 = c(1,1,1), apol = c(3,4,5))
    d <- cbind(M, X)
    out <- preprocess_data(d, add_cds = FALSE, rm_near_zero_var = TRUE, verbose = 0)
    expect_true("Fsp3" %notin% colnames(out))
    expect_true("apol" %in% colnames(out))
})

test_that("adding interaction terms requires >=2 predictors", {
    M <- data.frame(NAME = c("A","B"), SMILES = c("CC","CCC"), RT = c(1,2))
    X <- data.frame(Fsp3 = c(1,2))
    d <- cbind(M, X)
    out <- preprocess_data(d, add_cds = FALSE, interaction_terms = TRUE, verbose = 0)
    expect_true(!any(grepl("*", colnames(out), fixed = TRUE)))
})

test_that("function removes unsupported columns", {
    d <- data.frame(NAME = "A", SMILES = "CC", RT = 1, OTHER = 5)
    out <- preprocess_data(d, add_cds = FALSE, verbose = 0)
    expect_true(!"OTHER" %in% colnames(out))
})

test_that("polynomial terms are ignored when constructing interaction terms", {
    M <- data.frame(NAME = c("A","B"), SMILES = c("CC","CCC"), RT = c(1,2))
    X <- data.frame(Fsp3 = c(1,2), apol = c(3,4), nAtom = c(5,6))
    d <- cbind(M, X)
    out <- preprocess_data(
        data = d, add_cds = FALSE, degree_polynomial = 2,
        interaction_terms = TRUE, verbose = 0
    )
    expect_identical(colnames(out), c(
        "NAME", "SMILES", "RT",
        "Fsp3", "apol", "nAtom",
        "Fsp3^2", "apol^2", "nAtom^2",
        "Fsp3:apol", "Fsp3:nAtom", "apol:nAtom"
    ))
})
