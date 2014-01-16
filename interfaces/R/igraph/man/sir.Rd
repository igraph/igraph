\name{sir}
\alias{sir}
\concept{SIR model}
\concept{Dynamics on graph}
\title{SIR model on graphs}
\description{Run simulations for an SIR (susceptible-infected-recovered)
  model, on a graph}
\usage{
sir(graph, beta, gamma, no.sim=100)
}
\arguments{
  \item{graph}{The graph to run the model on. If directed, then edge
    directions are ignored and a warning is given.}
  \item{beta}{Non-negative scalar. The rate of infection of an
    individual that is susceptible and has a single infected
    neighbor. The infection rate of a susceptible individual with n
    infected neighbors is n times beta. Formally
    this is the rate parameter of an exponential distribution.}
  \item{gamma}{Positive scalar. The rate of recovery of an infected
    individual. Formally, this is the rate parameter of an exponential
    distribution.}
  \item{no.sim}{Integer scalar, the number simulation runs to perform.}
}
\details{
  The SIR model is a simple model from epidemiology. The individuals
  of the population might be in three states: susceptible, infected
  and recovered. Recovered people are assumed to be immune to the 
  disease. Susceptibles become infected with a rate that depends on
  their number of infected neigbors. Infected people become recovered
  with a constant rate. 
}
\value{
  The results are returned in an object of class \sQuote{\code{sir}},
  which is a list, with one element for each simulation. Each
  simulation is itself a list with the following elements. They are all
  numeric vectors, with equal length:
  \itemize{
    \item{times}{The times of the events.}
    \item{NS}{The number of susceptibles in the population, over time.}
    \item{NI}{The number of infected individuals in the population, over
      time.}
    \item{NR}{The number of recovered individuals in the population,
      over time.}
  }
}
\references{
  Bailey, Norman T. J. (1975). The mathematical theory of infectious
  diseases and its applications (2nd ed.). London: Griffin.
}
\author{
  Gabor Csardi \email{csardi.gabor@gmail.com}. Eric Kolaczyk
  (\url{http://math.bu.edu/people/kolaczyk/}) wrote the initial
  version in R.
}
\seealso{\code{\link{plot.sir}} to conveniently plot the results}
\examples{
g <- erdos.renyi.game(100, 100, type="gnm")
sm <- sir(g, beta=5, gamma=1)
plot(sm)
}
\keyword{graphs}
