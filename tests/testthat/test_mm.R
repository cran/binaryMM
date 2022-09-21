test_that("mm returns an object of class MMLong", {

  N      = 100
  nclust = sample( seq(10,10), N, replace=TRUE)
  id     = rep(seq(N), nclust)
  Xe     = rep(rbinom(N,size=1,prob=.5), nclust) # binary exposure
  time   = unlist( sapply( as.list(nclust), function(ZZ) seq(ZZ)-1 ) )
  data   = data.frame(id, time, Xe)
  data   = data[order(data$id, data$time),]
  newdat = GenBinaryY(mean.formula=~time*Xe, t.formula=~1,
                      beta=c(-2.5, 0.25, 0.25, 0.1),
                      gamma=1, id=id, data=data, q=20, Yname = "binY")

  fit <- mm(binY~time*Xe, t.formula=~1, data=newdat, id=id, step.max=4, verbose=FALSE)

  expect_equal(class(fit), "MMLong")

})

test_that("Inputs of mm are correct", {

  N      = 100
  nclust = sample( seq(10,10), N, replace=TRUE)
  id     = rep(seq(N), nclust)
  Xe     = rep(rbinom(N,size=1,prob=.5), nclust) # binary exposure
  time   = unlist( sapply( as.list(nclust), function(ZZ) seq(ZZ)-1 ) )
  data   = data.frame(id, time, Xe)
  data   = data[order(data$id, data$time),]
  newdat = GenBinaryY(mean.formula=~time*Xe, t.formula=~1,
                      beta=c(-2.5, 0.25, 0.25, 0.1),
                      gamma=1, id=id, data=data, q=20, Yname = "binY")

  expect_error(mm(binY~time*Xe, data=newdat,
                    id=id, step.max=4, verbose=FALSE))

})
