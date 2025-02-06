context("gsplit and rsplit")



wld2 <- wlddev
oldClass(wld2) <- NULL
vlabels(wld2) <- NULL
f <- wld2$iso3c
ind <- 1:1000
fss <- f[ind]
fl <- wld2[c("region", "income")]
flss <- ss(fl, ind)

test_that("gsplit / rsplit work like split", {

  for(i in seq_col(wld2)) {
    expect_equal(gsplit(wld2[[i]], f, TRUE), split(wld2[[i]], f))
    expect_equal(gsplit(wld2[[i]], f, FALSE), `names<-`(split(wld2[[i]], f), NULL))
    expect_equal(gsplit(wld2[[i]][ind], fss, TRUE), split(wld2[[i]][ind], fss))
    expect_equal(rsplit(wld2[[i]][ind], fss), split(wld2[[i]][ind], fss, drop = TRUE))
    # factor list
    expect_true(all_obj_equal(gsplit(wld2[[i]], fl, TRUE),
                              rsplit(wld2[[i]], fl, flatten = TRUE),
                              unlist(rsplit(wld2[[i]], fl), recursive = FALSE),
                              split(wld2[[i]], fl, drop = TRUE, lex.order = TRUE)))

    expect_true(all_obj_equal(gsplit(wld2[[i]][ind], flss, TRUE),
                              rsplit(wld2[[i]][ind], flss, flatten = TRUE),
                              unlist(rsplit(wld2[[i]][ind], flss), recursive = FALSE),
                              split(wld2[[i]][ind], flss, drop = TRUE, lex.order = TRUE)))
  }
})

test_that("rsplit matrix method works as intended", {
  m = qM(nv(GGDC10S))
  fl = lapply(GGDC10S[c("Country", "Variable")], qF, sort = FALSE)
  expect_equal(lapply(rsplit(m, GGDC10S$Country), unattrib), split(m, GGDC10S$Country))
  expect_equal(lapply(rsplit(m, itn(fl), flatten = TRUE), unattrib), split(m, itn(fl)))

  expect_equal(rsplit(m, fl, flatten = TRUE), unlist(rsplit(m, fl), FALSE))

  expect_true(all(vapply(rsplit(m, c(fl, GGDC10S["Year"]), flatten = TRUE), is.matrix, TRUE)))
  expect_true(!any(vapply(rsplit(m, c(fl, GGDC10S["Year"]), flatten = TRUE, drop.dim = TRUE), is.matrix, TRUE)))
})

test_that("rsplit data frame method works as intended", {

  expect_equal(rsplit(mtcars, mtcars$cyl), split(mtcars, mtcars$cyl))
  expect_equal(rsplit(mtcars, mpg ~ cyl), split(mtcars$mpg, mtcars$cyl))
  expect_equal(rsplit(mtcars, ~ cyl, nest = TRUE)$data, rsplit(mtcars, ~ cyl))
  expect_equal(rsplit(mtcars, mpg ~ cyl, simplify = FALSE), split(mtcars["mpg"], mtcars$cyl))

  expect_true(all_obj_equal(rsplit(mtcars, mtcars[.c(cyl, vs, am)], flatten = TRUE),
               rsplit(mtcars, ~ cyl + vs + am, flatten = TRUE, keep.by = TRUE),
               unlist(unlist(rsplit(mtcars, mtcars[.c(cyl, vs, am)]), FALSE), FALSE),
               unlist(unlist(rsplit(mtcars, ~ cyl + vs + am, keep.by = TRUE), FALSE), FALSE),
               split(mtcars, mtcars[.c(cyl, vs, am)], drop = TRUE, lex.order = TRUE)))

  expect_true(all_obj_equal(rsplit(mtcars, ~ cyl + vs + am, flatten = TRUE),
                            unlist(unlist(rsplit(mtcars, ~ cyl + vs + am), FALSE), FALSE),
                            split(mtcars[names(mtcars) %!in% .c(cyl, vs, am)],
                                  mtcars[.c(cyl, vs, am)], drop = TRUE, lex.order = TRUE)))

})




