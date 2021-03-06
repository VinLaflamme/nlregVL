% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nlregVL.R
\name{graFromFittedDT}
\alias{graFromFittedDT}
\title{Construct graphs from data tables containing psyFitted objects}
\usage{
graFromFittedDT(dataT, byC, psyFCol, Xlab = "Duration (seconds)",
  Ylab = "Proportion of Long Responses", Annotate = T)
}
\arguments{
\item{dataT}{The data.table containing the psyFitted objects}

\item{byC}{Name of the columns to contain the factors. There should be no duplicate rows.Can be a character vector or can be a single string with each name separated by commas.}

\item{psyFCol}{Name of the column containing the psyFitted objects}

\item{Xlab}{X axis label}

\item{Ylab}{y axis label}

\item{Annotate}{If TRUE, graphs will be annotated.}
}
\value{
A data table with the byC columns and a column containing ggplots
}
\description{
Construct graphs from data tables containing psyFitted objects
}
\examples{
data("exampleData",package="nlregVL")
result<-nlregVL::fitLSDT(dataF = exampleData,cD = "Duration",cNbS = "Nb1",cNbL = "Nb3",cBy = "Sub")
result2<-graFromFittedDT(dataT = result,byC = "Sub",psyFCol = "psyF")

}

