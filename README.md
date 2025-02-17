# Point Dragger and Geometry Visualizer

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [How to Run](#how-to-run)
- [Algorithms and Time Complexity](#algorithms-and-time-complexity)
- [Code Structure](#code-structure)

## Overview

This project is a Python Tkinter GUI application that allows you to:

- **Load** or **generate random points** on a canvas.
- **Drag points** interactively to change their positions.
- **Generate convex hull layers** for the set of points using the Graham scan algorithm.
- **Generate iterative minimum enclosing circles** (using a randomized MINIDISC algorithm) that enclose the remaining points.
- **Export the visualization as SVG files** in several formats:
  - **points.svg:** Contains only the input points (drawn as small black circles).
  - **hull.svg:** Overlays convex hull layers (drawn in alternating colors such as red and blue) with points colored according to the hull they belong to.
  - **cir.svg:** Overlays iterative minimum enclosing circles (drawn in alternating colors, e.g., neon-green and orange) with points colored according to the circle they are part of.
  - **all.svg:** Combines hulls, circles, and points (drawn in that order) into a single file.

![Point Dragger and Geometry Visualizer](./assets/point-dragger-gui.png)

## Features

- **Point Colors**  
  Each point is assigned a random color when loaded or generated. When saving to SVG, points are recolored to match the hull or circle layer they belong to.

- **Interactive Dragging**  
  Click and drag any point to reposition it on the canvas.

- **Convex Hull Generation**  
  The application computes convex hull layers (outer boundaries) using the Graham scan algorithm. Each hull is drawn in a random color.

- **Minimum Enclosing Circles**  
  The application iteratively computes the minimum enclosing circle (using a randomized MINIDISC algorithm) for all points, draws the circle, removes the enclosed points, and repeats for the remaining points.

## How to Run

1. **Requirements**
   - Python 3.x installed.
   - No external libraries are needed (only the standard library).

2. **Running the Application**
   - Open a terminal or command prompt.
   - Navigate to the project directory.
   - Run the command:
     ```
     python gui.py
     ```
3. **Using the Application**
   - **Load Points**  
     Click **Load Points from File** and select a text file formatted as follows:
     - The first line contains the number of points (n).
     - The following n lines each contain two integers (x and y coordinates) separated by a space.
   - **Generate Random Points**  
     Enter the desired number of points in the **Number of Points:** field, then click **Generate Random Points.**  
     Points are generated within the central portion of the canvas (with margins) to ensure clarity.
   - **Resize Canvas**  
     Change the canvas width and height in the controls and click **Resize** to scale the visualization.
   - **Generate Convex Hulls**  
     Click **Generate Convex Hulls** to compute and display convex hull layers.
   - **Generate Enclosing Circles**  
     Click **Generate Enclosing Circles** to iteratively compute and display minimum enclosing circles.
   - **Save SVG Files**  
     Click **Save SVG Files** to export the current visualization into four SVG files: `points.svg`, `hull.svg`, `cir.svg`, and `all.svg`.

## Algorithms and Time Complexity

### Graham Scan (Convex Hull)
- **Algorithm**  
  The Graham scan algorithm sorts the points based on their polar angle relative to the lowest point and then constructs the convex hull by iteratively checking the turn direction (using cross products).
- **Time Complexity**  
  O(n log n), dominated by the sorting step.

### MINIDISC (Minimum Enclosing Circle)
- **Algorithm**  
  The MINIDISC algorithm is a randomized incremental algorithm. It processes points in random order and updates the minimum enclosing circle when a point lies outside the current circle. Helper functions `MINIDISCWITHPOINT` and `MINIDISCWITH2POINTS` ensure that specific boundary conditions are met.
- **Time Complexity**  
  Expected O(n) for the randomized algorithm (though the worst-case can be higher).

## Code Structure

- **Geometry Class**  
  Contains all static methods for geometric computations, including:
  - Convex hull computation using the Graham scan.
  - Minimum enclosing circle (MINIDISC) computation with helper functions.
  - Utility functions to check if a point is inside or on the boundary of a circle.

- **PointDragger Class**  
  Handles the GUI logic for:
  - Point management (loading, generating, dragging, and resizing).
  - Visualizing convex hulls and minimum enclosing circles.
  - Exporting the visualization to SVG files using helper functions that write points, convex hulls, and circles to file.