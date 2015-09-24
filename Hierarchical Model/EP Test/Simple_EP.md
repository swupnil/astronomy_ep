---
title: "Simple EP"
author: "Swupnil Sahai"
date: "July 27, 2015"
output: pdf_document
---

\section{Uknown Mean}
\subsection{Formulation}
We start by fitting one of the simplest possible models that isn't hierarchical in nature:

$$ \phi \sim N(\mu_0, \sigma_0^2) $$
$$ x_i \sim N(\phi) $$

We then generated 1000 data points with $\phi = 200$. We gave EP and HMC priors of $\mu_0 = 0$ and $\sigma_0 = 20$.

\subsection{Results}
EP and HMC converge to similar posterior centers of $\phi$, but the variances are quite different.



