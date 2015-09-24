Expectation Propagation with Astronomy Data
====================

This project explores data collected by the GALEX ultraviolet space telescope in an effort to model the relationship between infrared and ultrioviolet radiation in various parts of the all sky map. We attempt to implement EP (Expectation Propagation) in fitting a hierarchical model with nine global parameters. EP is an algorithm that approximates the posterior distribution with a multivariate normal and, most crucially, allows inference to be conducted simulatenously on different partitions of the data set. This makes posterior inference highly parallelizable and has massive implications for scalable Bayesian data analysis.
