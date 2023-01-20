# Summary
This code simulates the cell differentiation model introduced in "Oscillatory Protein Expression Dynamics Endows Stem Cells with Robust Differentiation Potential" (*Nature* 470.7334 (2011), pp. 329–334) by Narito Suzuki, Chikara Furusawa and Kunihiko Kaneko. Stem cells are modelled as having gene and protein expression levels whose change is governed by two coupled differential equations and the inter-cellular interaction through a medium. I implemented this as part of the lecture "Universal Biology", held by Professor Furusawa, at the University of Tokyo in the winter semester of 2022/23.

# Code Structure
The code carrying out the simulation is contained in ``simulation.cpp``, which links to the header files ``functions.h`` (containing the cell struct and functions implementing time evolution and mitosos) and ``infrastructure.h`` (containing output routines and operator overloads). The gene and protein expression timelines obtained through simulation are saved in text files. ``evaluation.ipynb`` creates plots and tests for differentiated cell types based on these results.
