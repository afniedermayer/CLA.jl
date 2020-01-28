# README #

This repository contains a reference implementation of the Critical Line Algorithm as described in [Niedermayer and Niedermayer (Applying Markowitz's Critical Line Algorithm)](http://www.springer.com/business+%26+management/finance/book/978-0-387-77438-1), 2010
Handbook of Portfolio Construction: Contemporary Applications of
Markowitz Techniques", edited by John Guerard, London: Springer, 2010, see also [the working paper version](http://www.vwl.unibe.ch/papers/dp/dp0701.pdf).

The reference implementation is intended to be as close as possible to the pseudo-code in the paper, even if this comes at the cost of performance.
For a fast implementation, see the [repository with the Fortran version](https://bitbucket.org/afniedermayer/fast_critical_line_algorithm) (which also
contains Python bindings for the Fortran code).

The IPython notebooks "example_call.ipynb" and "example_call_gen.ipynb" contain sample code that calls the Julia functions. The notebooks also compare
the results with the Python implementation (which just calls Fortran). To be able to call the Python code from the Julia notebooks, you need to get the [code](https://bitbucket.org/afniedermayer/fast_critical_line_algorithm) and set the path in the notebook to point to it.

*Note:* This code was written in 2014 using [Julia 0.3](https://julialang.org/downloads/oldreleases/). It most likely won't work on newer version of Julia. Feel free to make a pull request if you manage to make the code run on Julia 1.0+.