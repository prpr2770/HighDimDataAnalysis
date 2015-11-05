Information-theoretic Metric Learning (ITML-1.1) README File
=================================================

1. Quick-start

To test the code, we have provided a demo Matlab script Test.m, which runs metric learning with cross-validation and k-nearest neighbors over the example data set Iris.  Simply type 'Test' in the Matlab command prompt.  If everything is working properly, the accuracy of the nearest neighbor classifier using the learned metric should be around 95%.

2. Using the code

For fully supervised cases, metric learning is generally invoked after creating constraints from the labels---similarity constraints are generated for pairs of points in the same class, and dissimilarity constraints are generated for pairs of points in different classes.  The function MetricLearnAutotuneKnn.m creates random constraints from a labeled data set, runs metric learning while varying the gamma parameter, and returns a Mahalanobis matrix A.  When a nearest neighbor classifier is desired, the function CrossValidateKNN.m can be used to perform classification with cross-validation.  

For example, to generate cross-validated KNN accuracy on a data set X with labels y using standard Euclidean distance, you can use the following command:

acc=CrossValidateKNN(y,X,@(y,X) EuclideanDistance(X),num_folds,knn_neighbor_size);

The EuclideanDistance.m function is provided.

In many simple cases, the Test.m script can easily be modified to handle other labeled data sets.  In other cases, such as when there are user-defined constraints, the ItmlAlg.m program can be invoked directly.  As input, this function requires the constraints, the data matrix, the initial Mahalanobis matrix, and a set of parameters.  See ItmlAlg.m for information on the format for these inputs.

3. Citation

Finally, please acknowledge the use of our code with a citation:

Jason V. Davis, Brian Kulis, Prateek Jain, Suvrit Sra, and Inderjit S. Dhillon.  "Information-theoretic Metric Learning."  Proc. 24th International Conference on Machine Learning (ICML), 2007.


This software is currently being actively maintained; please check back often for updates. When using this code, please cite ITML and the relevant paper.

@inproceedings{davis-et-al-icml-2007,
  author    = {Jason V. Davis and
               Brian Kulis and
               Prateek Jain and
               Suvrit Sra and
               Inderjit S. Dhillon},
  title     = {Information-theoretic metric learning},
  booktitle = {ICML},
  year      = {2007},
  pages     = {209-216},
  address   = {Corvalis, Oregon, USA}
  month     = {June}
}

@manual{itml,
  title        = {Information Theoretic Metric Learning},
  author       = {Jason V. Davis and Brian Kulis and Prateek Jain and Suvrit Sra and Inderjit S. Dhillon},
  organization = {UT, Austin},
  address      = {http://www.cs.utexas.edu/users/pjain/itml/},
}