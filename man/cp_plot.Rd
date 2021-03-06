% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{cp_plot}
\alias{cp_plot}
\title{Cluster participation plots}
\usage{
cp_plot(Panmatrix, show_cluster, plot_type, use_log)
}
\arguments{
\item{Panmatrix}{Panmatrix produced by make_panmatrix functions}

\item{show_cluster}{Integer value  optional parameter allows to user to "zoom",
ignoring clusters that have organism participation below it}

\item{plot_type}{Can be either line or bar}

\item{use_log}{Logical:Scale axis containing clusters with log}
}
\description{
This function produces a basic plot of the pangenome.
   This can be either a point plot or a bar plot.
}
\examples{
\dontrun{cp_plot(Panmatrix)}

}
