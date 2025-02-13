# Point Dragger and Geometry Visualizer

This project is a Python Tkinter GUI application that allows you to:
- **Load** or **generate random points** on a canvas.
- **Drag points** interactively to change their positions.
- **Generate convex hull layers** for the set of points using the Graham scan algorithm.
- **Generate iterative minimum enclosing circles** (using the MINIDISC algorithm) that enclose the remaining points.

## Features

- **Point Colors:**  
  Each point is assigned a random color when loaded or generated, and this color is retained when drawing.

- **Interactive Dragging:**  
  Click and drag any point to reposition it on the canvas.

- **Convex Hull Generation:**  
  The application computes convex hull layers (outer boundaries) using the Graham scan algorithm. Each hull is drawn in a random color.

- **Minimum Enclosing Circles:**  
  The application iteratively computes the minimum enclosing circle (using a randomized MINIDISC algorithm) for all points, draws the circle, removes the enclosed points, and repeats for the remaining points.

## How to Run

1. **Requirements:**
   - Python 3.x installed.
   - No external libraries are needed (only the standard library).

2. **Run the Application:**
   - Open a terminal or command prompt.
   - Navigate to the project directory.
   - Run the command:
     ```
     python gui.py
     ```
3. **Using the Application:**
   - **Load Points:** Click on "Load Points from File" and select a text file with the following format:
     - The first line contains the number of points (n).
     - The following n lines each contain two integers (x and y coordinates) separated by a space.
   - **Generate Random Points:** Enter the desired number of points in the "Number of Points:" field, then click "Generate Random Points."
   - **Resize Canvas:** Change the width and height values in the canvas size controls and click "Resize."
   - **Generate Convex Hulls:** Click "Generate Convex Hulls" to compute and display convex hull layers.
   - **Generate Enclosing Circles:** Click "Generate Enclosing Circles" to iteratively compute and draw minimum enclosing circles.

## Algorithms and Time Complexity

### Graham Scan (Convex Hull)
- **Algorithm:**  
  The Graham scan algorithm sorts the points based on their polar angle relative to the lowest point and then constructs the convex hull by iteratively checking the turn direction (using cross products).
- **Time Complexity:**  
  O(n log n), dominated by the sorting step.

### MINIDISC (Minimum Enclosing Circle)
- **Algorithm:**  
  The MINIDISC algorithm is a randomized incremental algorithm. It processes points in random order and updates the minimum enclosing circle when a point lies outside the current circle. Helper functions `MINIDISCWITHPOINT` and `MINIDISCWITH2POINTS` ensure that specific boundary conditions are met.
- **Expected Time Complexity:**  
  Expected O(n) for the randomized algorithm, though the worst-case can be higher.

## Code Structure

- **Geometry Class:**  
  Contains all static methods related to geometric computations (convex hull, minimum enclosing circle, helper functions).

- **PointDragger Class:**  
  Handles the GUI logic, including point management, dragging, canvas resizing, and invoking geometry algorithms for visualization.
