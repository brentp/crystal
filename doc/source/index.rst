.. crystal documentation master file, created by
   sphinx-quickstart on Mon Sep  8 11:12:44 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

crystal documentation
=====================

`crystal` is a series python modules to evaluate modeling strategies for correlated data.

Practically, `crystal` is a python module to assign significance to clusters of
correlated data. The most likely use-case is for DNA methylation data, but
it can be use for any correlated data.

New users should check-out the introductory IPython notebook `here <http://nbviewer.ipython.org/github/brentp/crystal/blob/master/notebooks/Introduction.ipynb>`_.

Developers interested in evaluating their own method should view `this notebook <http://nbviewer.ipython.org/github/brentp/crystal/blob/master/notebooks/crystal-methods-evaluation.ipynb>`_.

API Docs
========

.. currentmodule:: crystal

.. toctree::
   :maxdepth: 2

   crystal
   plot
   evaluate
   regions

Introduction
============

crystal models clusters of correlated data.

quick-start
-----------

All functions to model clusters take:

1. formula: a string model formula in R (or `patsy <http://patsy.readthedocs.org/en/latest/>`_ syntax). 
2. a methylation cluster from `aclust <https://github.com/brentp/aclust>`_
3. a pandas datafame containing the covariates of interest. The samples in this
   dataframe should be in the **same order** as they are in the methylation cluster.
4. a coefficient of interest from the model.

So a single call would look like:

.. code-block:: python

    import pandas as pd
    import aclust
    import crystal

    covs_df = pd.read_csv('covariates.csv')
    cluster_iter = aclust.mclust(feature_gen(), max_dist=100)
    crystal.zscore_cluster("methylation ~ disease + age + gender",
                            next(cluster_iter),
                            covs_df,
                            "disease")

Which will return a dictionary containing the `p` value, the `coef`
and the `t` statistic for 'disease' status. However it is more likely
that users will want a simple way to model every cluster:

.. code-block:: python 

    for cluster in crystal.model_clusters(cluster_iter,
                                          covs_df,
                                          formula,
                                          "disease",
                                          crystal.zscore_cluster,
                                          n_cpu=10):
        print cluster

Example
=======

Here we run through an actual example using data in the repository. It would be simple
to modify this script to run on any data.

.. plot:: _static/ex-pipeline.py
    :include-source: true

Note that we have this looks like a site where males generally have a lower level of 
methylation than females. The text output looks like this::

    chrom	start	end	coef	n_sites	
    chrX	2733163	2745940	0.0002336	-0.514	3
    chrX	2746332	2746420	0.006272	0.257	2
    chrX	2746696	2746697	0.5873	-0.077	1
    chrX	2747135	2747136	0.5996	0.054	1
    chrX	2800456	2800457	0.6346	0.041	1
    chrX	2822260	2822261	0.03331	0.168	1
    chrX	2825269	2825270	0.0005556	-0.361	1
    chrX	2825361	2825362	0.05747	-0.168	1
    chrX	2825595	2825596	0.8458	-0.039	1
    chrX	2826828	2835913	0.001088	-0.254	2
    chrX	2836026	2836084	0.1921	0.148	2
    chrX	2836113	2836114	0.8866	-0.024	1
    chrX	2836299	2836300	0.001849	-0.288	1
    chrX	2843195	2843196	0.01105	-0.173	1
    chrX	2844680	2844681	0.07179	-0.272	1
    chrX	2846195	2846196	0.5473	-0.075	1
    chrX	2846892	2846893	0.13	-0.101	1
    chrX	2847353	2847503	0.008624	-0.437	3
    chrX	2847509	2847510	0.81	0.039	1
    chrX	2847548	2852845	0.0003622	-0.490	6

Development
===========

`crystal` is developed at `https://github.com/brentp/crystal <https://github.com/brentp/crystal/>`_.

It includes example data in the `crystal/tests` directory.

It has decent test coverage and tests can be run with `./test.sh` from the root of the project
directory.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

