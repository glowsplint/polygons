# Polygons

### The Problem

Consider an n-sided polygon where n is even. We draw all the lines that connect the vertices, and mark out the intersection points. This creates a graph with the vertices and intersection points as the nodes of the graph, and the edges are the line segments between the vertices and the intersection points. We define a start point and an end point as the opposite vertices of the polygon. We also define a "decreasing path" as a path along the edges of the graph that links the start point to the end point, where the distance to the end point is strictly decreasing while moving along the path. Find the total number of decreasing paths available.

### The Solution

A solution is provided in `polygons.ipynb`. We think that the growth rate in _n_ (where _n_ is the number of sides of the polygon) is _exp(n^3)_.
