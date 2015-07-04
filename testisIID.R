rm (list = ls())
graphics.off()

source ("./isIID.R")
set.seed (329588)
N <- 200
x1 <- 2 + rt (n=N, df=5)

ans <- isIID (x=x1)
print ("Testing x1")
if (ans$reject) {
	print (sprintf ("The test statistic is = %.4f, the %.0f%% critical value is %.4f, reject", 
                  ans$statistic, (100*ans$CL), ans$crit.val))
} else {
	print (sprintf ("The test statistic is = %.4f, the %.0f%% critical value is %.4f, can NOT reject",
		ans$statistic, (100*ans$CL), ans$crit.val)) 
}
print (sprintf ("p-value %.4f", ans$p.val))
print (Box.test (x1, lag=1, type="Ljung"))

cat ("\nTesting x2")
x2 <- x1
x2[1:(N/2)] <- rep(2,(N/2))
ans <- isIID (x=x2)
if (ans$reject) {
  print (sprintf ("The test statistic is = %.4f, the %.0f%% critical value is %.4f, reject", 
           ans$statistic, (100*ans$CL), ans$crit.val))
} else {
  print (sprintf ("The test statistic is = %.4f, the %.0f%% critical value is %.4f, can NOT reject",
           ans$statistic, (100*ans$CL), ans$crit.val) )
}
print (sprintf ("p-value %.4f", ans$p.val))
print (Box.test (x2, lag=1, type="Ljung"))

library (signal)

c <- 0.1
cat ("\nIntroduce dependency by filtering\n")
tmp <- filter (filt=c(c), a=c(1,c-1), x=x1)
x3 <- as.vector (tmp)

print ("Test x3")
ans <- isIID (x=x3)
if (ans$reject) {
  print (sprintf ("The test statistic is = %.4f, the %.0f%% critical value is %.4f, reject", 
           ans$statistic, (100*ans$CL), ans$crit.val))
} else {
  print (sprintf ("The test statistic is = %.4f, the %.0f%% critical value is %.4f, can NOT reject",
           ans$statistic, (100*ans$CL), ans$crit.val) )
}
print (sprintf ("p-value %.4f", ans$p.val))
print (Box.test (x3, lag=1, type="Ljung"))

c <- 0.2
cat ("\nIntroduce dependency by filtering\n")
tmp <- filter (filt=c(c), a=c(1,c-1), x=x1)
x4 <- as.vector (tmp)

print ("Test x4")
ans <- isIID (x=x4)
if (ans$reject) {
  print (sprintf ("The test statistic is = %.4f, the %.0f%% critical value is %.4f, reject", 
           ans$statistic, (100*ans$CL), ans$crit.val))
} else {
  print (sprintf ("The test statistic is = %.4f, the %.0f%% critical value is %.4f, can NOT reject",
           ans$statistic, (100*ans$CL), ans$crit.val) )
}
print (sprintf ("p-value %.4f", ans$p.val))
print (Box.test (x4, lag=1, type="Ljung"))

