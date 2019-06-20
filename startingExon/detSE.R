
testfactor <- factor(rep(c("A","B"), c(3,2)))
print(testfactor)

# test of histogram
set.seed(1234)
dat <- data.frame(cond = factor(rep(c("A","B"), each=200)),
                  rating = c(rnorm(200),rnorm(200, mean=.8)))
head(dat)
library(ggplot2)
ggplot(dat, aes(x=rating, fill=cond)) +
  geom_histogram(binwidth=.5, position="dodge")

