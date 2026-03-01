\name{Peirce}
\alias{Peirce}
\docType{data}
\title{C.S. Peirce's Auditory Response Data}
\description{
Data from sequence experiments conducted by C.S. Pierce in 1872 to
determine the distribution of response times to an auditory stimulus.
}
\usage{data(Peirce)}
\format{
A \code{link{list}} of 24 objects each representing one day of the
experiment.  Each element of the list consists of three components:
the date the measurements were made, an \code{\link{x}} component
recording the response time in milliseconds, and an associated  \code{\link{y}}
component recording a count of the number of times that the response
was recorded to be equal to be equal to the  corresponding \code{\link{x}} entry.
There  are roughly 500 observations (counts) on each of the 24 days.
}
\details{
A detailed description of the experiment can be found in Peirce (1873).  A
young man of about 18 with no prior experience was  employed to respond to
a signal  ``consisting of a sharp sound like a rap, the answer being made
upon a telegraph-operator's key nicely adjusted.''  The response times,
made with the aid of a Hipp cronoscope were recorded to the nearest 
millisecond.  The data was analyzed by Peirce who concluded that after
the first day, when the the observer was entirely inexperienced, the
curves representing the densities of the response times ``differed very
little from that derived from the theory of least squares,'' i.e. from
the Gaussian density.

The data was subsequently analysed by Samama, in a diploma thesis supervised
by Maurice Frechet, who reported briefly the findings in Frechet (1924),
and by Wilson and Hilferty (1929).  In both instances the reanalysis showed
that Laplace's first law of error, the double exponential distribution, was
a better representation for the data than was the Gaussian law.

The data is provide in two formats:  The first in a "raw" form as 24 text
files as scanned from the reprinted Peirce source, the second as an R
dataset \code{\link{Peirce.rda}} containing the list.  See the examples
section below for details on how the text files were converted to the list,
and a provisional  attempt to reproduce part of the Wilson and Hilferty
analysis.  An open question regarding the dataset is:  How did Wilson
and Hilferty compute standard deviations for the median as they appear
in their table?  The standard textbook suggestion of Yule (1917) yields
far too small a bandwidth.  The methods employed in the example section
below, which rely on relatively recent proposals, are somewhat closer, 
but still deviate somewhat from the results reported by Wilson and Hilferty.
}
\source{
\textsc{Peirce, C.~S.}  (1873): ``On the Theory of Errors of Observation,''
  \emph{Report of the Superintendent of the U.S. Coast Survey}, pp. 200--224,
  Reprinted in {\it The New Elements of Mathematics}, (1976) collected papers
  of C.S. Peirce, ed. by C. Eisele, Humanities Press: Atlantic Highlands, N.J.,
  vol. 3, part 1, 639--676.
}
\references{
\textsc{Fr\'echet, M.}  (1924): ``Sur la loi des erreurs d'observation,''
  \emph{Matematichiskii Sbornik}, 32, 5--8.
\textsc{Wilson, E.~B.,  {\small and} M.~M. Hilferty}  (1929): ``Note on C.S.
  Peirce{'}s Experimental Discussion of the Law of Errors,'' \emph{Proceedings
  of the National Academy of Sciences of the U.S.A.}, 15, 120--125.
\textsc{Yule, G.~U.}  (1917): \emph{An Introduction to the Theory of
  Statistics}. Charles Griffen: London, 4 edn.
}
\examples{
data(day)

# Read Peirce's tables

#Usage:  Peirce <- ReadPeirceData()
ReadPeirceData <- function(){
  ReadPeirceFile <- function(file){
	d <- scan(file)
	n <- length(d)
	m <- (n %/% 18) + 1
	d[(n+1):(18*m)] <- NA
	D <- matrix(d,2)
	D <- t(matrix(D,18))
	D <- array(D,c(m,2,9))
	D <- aperm(D,c(1,3,2))
	x <- c(D[,,1])
	y <- c(D[,,2])
	x <- x[!is.na(x)]
	y <- y[!is.na(y)]
	for(i in 2:length(x))
		if(x[i] < x[i-1]) x[i] = x[i-1] + 1
	s <- y==0
	list(x = x[!s], y = y[!s])
	}
  files <- list.files(pattern = "txt")
  ndays <- length(files)
  day   <- list(NULL)
  dates <- c(1,5,6,8,9,10,15,16,17,18,19,20,22,23,24,25,26,27,29,30,31,1,2,3)
  dates <- paste(c(rep("July",21),rep("August",3)),dates,"1872")

  for(i in 1:ndays){
	print(files[i])
	day[[i]] <- ReadPeirceFile(files[i])
	day[[i]]$date <- dates[i]
	}
  day
  }
# Make table like Wilson and Hilferty

load("Peirce.rda")
set.seed(10) #Dither the counts
tab <- matrix(0,24,11)
for(i in 1:24){
	y <- rep(day[[i]]$x, day[[i]]$y) + runif(sum(day[[i]]$y), -.5, .5)
	f1 <- summary(rq(y~1),se="iid")$coef[1:2]
	n <- length(y)
	f0 <- 1/(2 * sum(abs(y-f1[1])/n)) #Laplace proposal
	f0 <- (1/(2 * f0))/ sqrt(n)
	f2 <- summary(lm(y~1))$coef[1:2]
	outm <- sum(y < (f1[1] - 3.1 * sqrt(n) * f2[2]))
	outp <- sum(y > (f1[1] + 3.1 * sqrt(n) * f2[2]))
	outt <- outm + outp
	inm <- y > (f1[1] - 0.25 * sqrt(n) * f2[2])
	inp <- y < (f1[1] + 0.25 * sqrt(n) * f2[2])
	int <- sum(inm * inp)
	Eint <- round(n * (pnorm(.25) - pnorm(-.25)))
	excess <- round(100*(int - Eint)/Eint)
	tab[i,] <- c(f1, f0, f2, outm, outp, outt,int,Eint,excess)
	cnames <- c("med","sdmed1","sdmed0","mean","sdmean","below","above","outliers",
		"inliers","Einliers","ExcessIns")
	dimnames(tab) <- list(paste("Day",1:24),cnames)
	}
\keyword{datasets}
