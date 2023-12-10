This section contains Python scripts that can be used to calculate the circle number Pi. Once completed, one will find here Python scripts that are based directly on Archimedes' ideas.

It should be noted that Python works with a limited number of decimal places by default. This is also valid for square roots from the standard module 'math'. To avoid this problem, the standard Python module 'decimal' can be used. Take also a look at the docstrings in the Python scripts with respect to the aforementioned comment.

Based on Archimedes' approach, further algorithms were developed to calculate both the lower limit and the upper limit as well as Pi as the mean value of the lower and upper limits.

As soon as I have a little more time, I will explain in more detail how the corresponding equations can be derived from the geometric relationships. 

<img src="\images/archimedes_figure3.png" alt="Figure 3"><br/>
Figure 3: Archimedes, Measurement of a Circle, geometrical model for one edge of the outer regular polygon
<br/>

<img src="\images/archimedes_figure4.png" alt="Figure 3"><br/>
Figure 4: Archimedes, Measurement of a Circle, geometrical model for one edge of the inner regular polygon
<br/>

The Indian Aryabhatta, who is quoted in Bhaskara I, apparently knew the relationship to the iterative calculation of the side lengths of a regular polygon using recurrence. This relationship can be found in the known sources under the name Archimedes algorithm.

The also so-called Archimedes algorithm, which can be found in a lot of resources, is called Archimedean iterative algorithm, Archimedean mean iteration or Pfaff-Borchardt-Schwab algorithm as well as Archimedes recurrence formula.

The latter method is similar to the Gauss-Legendre algorithm. In this method two numbers are repeatedly replaced by their arithmetic and geometric mean in order to approximate their arithmetic-geometric mean.

Taking a look at Figure 3 and Figure 4 shows, that it is in principle not possible, to derive a recurrence relationship from two different geometrical approaches.

In the context of the Archimedes algorithm or Archimedes method a distinction must be made between recursive and recurrence or iterative algorithm. The wording recurrence or iterative is correct and recursive is not.
