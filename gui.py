import tkinter as tk
from tkinter import filedialog
import random
import math

class Geometry:
    """Collection of static methods for geometric computations.

    This class includes methods to compute convex hulls using the Graham Scan
    algorithm and to compute the minimum enclosing circle (MINIDISC) using a
    randomized incremental approach.
    """

    @staticmethod
    def graham_scan(points):
        """Compute the convex hull of a set of points using the Graham scan algorithm.
        
        Args:
            points (list): List of (x, y) tuples.
            
        Returns:
            list: The convex hull as a list of (x, y) tuples.
            
        Time Complexity:
            O(n log n) due to sorting.
        """
        def cross(o, a, b):
            return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])
        sorted_points = sorted(points, key=lambda p: (p[1], p[0]))
        lower = []
        for p in sorted_points:
            while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
                lower.pop()
            lower.append(p)
        upper = []
        for p in reversed(sorted_points):
            while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
                upper.pop()
            upper.append(p)
        return lower[:-1] + upper[:-1]

    @staticmethod
    def disc_two_points(p1, p2):
        """Return the smallest circle that passes through two points.
        
        Args:
            p1, p2 (tuple): Points as (x, y) tuples.
            
        Returns:
            tuple: (cx, cy, r) where (cx,cy) is the center and r is the radius.
            
        Time Complexity:
            O(1)
        """
        cx = (p1[0] + p2[0]) / 2.0
        cy = (p1[1] + p2[1]) / 2.0
        r = math.dist(p1, p2) / 2.0
        return (cx, cy, r)

    @staticmethod
    def disc_three_points(p1, p2, p3):
        """Return the circle passing through three non-collinear points.
        
        Args:
            p1, p2, p3 (tuple): Points as (x, y) tuples.
            
        Returns:
            tuple or None: (cx, cy, r) if a circle exists, else None.
            
        Time Complexity:
            O(1)
        """
        ax, ay = p1
        bx, by = p2
        cx, cy = p3
        d = 2 * (ax*(by-cy) + bx*(cy-ay) + cx*(ay-by))
        if d == 0:
            return None
        ux = ((ax**2+ay**2)*(by-cy) + (bx**2+by**2)*(cy-ay) + (cx**2+cy**2)*(ay-by)) / d
        uy = ((ax**2+ay**2)*(cx-bx) + (bx**2+by**2)*(ax-cx) + (cx**2+cy**2)*(bx-ax)) / d
        r = math.dist((ux, uy), p1)
        return (ux, uy, r)

    @staticmethod
    def is_point_in_disc(p, D):
        """Check if point p is inside or on circle D.
        
        Args:
            p (tuple): (x, y) point.
            D (tuple): (cx, cy, r) circle.
            
        Returns:
            bool: True if p is in or on D, else False.
        """
        if D is None:
            return False
        return math.dist(p, (D[0], D[1])) <= D[2] + 1e-6

    @staticmethod
    def is_point_on_disc(p, D):
        """Check if point p is exactly on the boundary of circle D.
        
        Args:
            p (tuple): (x, y) point.
            D (tuple): (cx, cy, r) circle.
            
        Returns:
            bool: True if p is on the boundary, else False.
        """
        if D is None:
            return False
        return math.isclose(math.dist(p, (D[0], D[1])), D[2])

    @staticmethod
    def MINIDISC(P):
        """Compute the minimum enclosing circle for set P using a randomized algorithm.
        
        Args:
            P (list): List of (x, y) points.
            
        Returns:
            tuple: (cx, cy, r) representing the minimum enclosing circle.
            
        Expected Time Complexity:
            Expected O(n), though worst-case can be higher.
        """
        P_copy = P.copy()
        random.shuffle(P_copy)
        if not P_copy:
            return None
        if len(P_copy) == 1:
            return (P_copy[0][0], P_copy[0][1], 0)
        D = Geometry.disc_two_points(P_copy[0], P_copy[1])
        for i in range(2, len(P_copy)):
            p = P_copy[i]
            if not Geometry.is_point_in_disc(p, D):
                D = Geometry.MINIDISCWITHPOINT(P_copy[:i], p)
        return D

    @staticmethod
    def MINIDISCWITHPOINT(P, q):
        """Compute the minimum enclosing circle for set P with q on its boundary.
        
        Args:
            P (list): List of (x, y) points.
            q (tuple): A point (x, y) that must lie on the circle's boundary.
            
        Returns:
            tuple: (cx, cy, r) for the minimum enclosing circle with q on its boundary.
        """
        P_copy = P.copy()
        random.shuffle(P_copy)
        D = Geometry.disc_two_points(q, P_copy[0])
        for i in range(1, len(P_copy)):
            p = P_copy[i]
            if not Geometry.is_point_in_disc(p, D):
                D = Geometry.MINIDISCWITH2POINTS(P_copy[:i], q, p)
        return D

    @staticmethod
    def MINIDISCWITH2POINTS(P, q1, q2):
        """Compute the minimum enclosing circle for set P with q1 and q2 on its boundary.
        
        Args:
            P (list): List of (x, y) points.
            q1, q2 (tuple): Points that must lie on the circle's boundary.
            
        Returns:
            tuple: (cx, cy, r) for the minimum enclosing circle with q1 and q2 on its boundary.
        """
        D = Geometry.disc_two_points(q1, q2)
        for p in P:
            if not Geometry.is_point_in_disc(p, D):
                D = Geometry.disc_three_points(q1, q2, p)
        return D


class PointDragger:
    """GUI application to visualize point dragging, convex hulls, and minimum enclosing circles.
    
    Also provides functionality to save the current output as SVG files and a text file summary.
    """
    def __init__(self, root, width=600, height=800):
        self.root = root
        self.root.title("Point Dragger")
        self.canvas_width = width
        self.canvas_height = height
        
        # --- Controls ---
        self.controls_frame = tk.Frame(root)
        self.controls_frame.pack(side=tk.TOP, fill=tk.X, pady=5)
        
        # Random points controls (Uniform, Gaussian, Exponential)
        self.random_frame = tk.Frame(self.controls_frame)
        self.random_frame.pack(fill=tk.X, pady=2)
        self.random_label = tk.Label(self.random_frame, text="Number of Points:")
        self.random_label.pack(side=tk.LEFT, padx=5)
        self.random_entry = tk.Entry(self.random_frame, width=5)
        self.random_entry.pack(side=tk.LEFT)
        self.uniform_button = tk.Button(self.random_frame, text="Generate Uniform Points", command=self.generate_uniform_points)
        self.uniform_button.pack(side=tk.LEFT, padx=5)
        self.gauss_button = tk.Button(self.random_frame, text="Generate Gaussian Points", command=self.generate_gaussian_points)
        self.gauss_button.pack(side=tk.LEFT, padx=5)
        self.exp_button = tk.Button(self.random_frame, text="Generate Exponential Points", command=self.generate_exponential_points)
        self.exp_button.pack(side=tk.LEFT, padx=5)
        
        # Load points button
        self.file_frame = tk.Frame(self.controls_frame)
        self.file_frame.pack(fill=tk.X, pady=2)
        self.load_button = tk.Button(self.file_frame, text="Load Points from File", command=self.load_points)
        self.load_button.pack(side=tk.LEFT, padx=5)
        
        # Canvas size controls
        self.size_frame = tk.Frame(self.controls_frame)
        self.size_frame.pack(fill=tk.X, pady=2)
        self.width_label = tk.Label(self.size_frame, text="Canvas Width:")
        self.width_label.pack(side=tk.LEFT, padx=5)
        self.width_entry = tk.Entry(self.size_frame, width=5)
        self.width_entry.insert(0, str(self.canvas_width))
        self.width_entry.pack(side=tk.LEFT)
        self.height_label = tk.Label(self.size_frame, text="Canvas Height:")
        self.height_label.pack(side=tk.LEFT, padx=5)
        self.height_entry = tk.Entry(self.size_frame, width=5)
        self.height_entry.insert(0, str(self.canvas_height))
        self.height_entry.pack(side=tk.LEFT)
        self.resize_button = tk.Button(self.size_frame, text="Resize", command=self.resize_canvas)
        self.resize_button.pack(side=tk.LEFT, padx=5)
        
        # --- Canvas ---
        self.canvas = tk.Canvas(root, width=self.canvas_width, height=self.canvas_height, bg="white")
        self.canvas.pack()
        
        # --- Additional Function Buttons ---
        self.output_frame = tk.Frame(root)
        self.output_frame.pack(fill=tk.X, pady=5)
        self.hull_button = tk.Button(self.output_frame, text="Generate Convex Hulls", command=self.draw_convex_hulls)
        self.hull_button.pack(side=tk.LEFT, padx=5)
        self.circle_button = tk.Button(self.output_frame, text="Generate Enclosing Circles", command=self.draw_enclosing_circles)
        self.circle_button.pack(side=tk.LEFT, padx=5)
        self.save_button = tk.Button(self.output_frame, text="Save SVG Files & Summary", command=self.save_svg)
        self.save_button.pack(side=tk.LEFT, padx=5)
        
        # --- Data ---
        self.points = []     # List of (x, y) tuples
        self.colors = []     # List of colors for points
        self.oval_ids = []   # Canvas IDs for points
        self.drag_data = {"x": 0, "y": 0, "idx": None, "initial_point": (0, 0)}
        
        self.canvas.bind("<ButtonRelease-1>", self.end_drag)
    
    # --- Point Generation with Minimum Distance Rejection ---
    def _min_distance(self, num_points):
        """Calculate the minimum allowed distance between points.
        
        Uses the canvas size and number of points.
        """
        return min(self.canvas_width, self.canvas_height) / (num_points)

    def generate_uniform_points(self):
        """Generate a random set of points using uniform randomness with rejection to avoid close points."""
        try:
            num_points = int(self.random_entry.get())
            if num_points <= 0:
                return
        except ValueError:
            return
        pts = []
        cols = []
        dmin = self._min_distance(num_points)
        max_tries = num_points * 1000
        tries = 0
        while len(pts) < num_points and tries < max_tries:
            x = random.randint(int(self.canvas_width/4), int(3*self.canvas_width/4))
            y = random.randint(int(self.canvas_height/4), int(3*self.canvas_height/4))
            candidate = (x, y)
            if all(math.dist(candidate, p) >= dmin for p in pts):
                pts.append(candidate)
                cols.append("#{:06x}".format(random.randint(0, 0xFFFFFF)))
            tries += 1
        self.points = pts
        self.colors = cols
        self.draw_points()

    def generate_gaussian_points(self):
        """Generate a set of points using Gaussian randomness centered at the canvas center."""
        try:
            num_points = int(self.random_entry.get())
            if num_points <= 0:
                return
        except ValueError:
            return
        pts = []
        cols = []
        # Set sigma based on canvas dimensions
        sigma_x = self.canvas_width / 8
        sigma_y = self.canvas_height / 8
        dmin = self._min_distance(num_points)
        max_tries = num_points * 1000
        tries = 0
        while len(pts) < num_points and tries < max_tries:
            x = int(random.gauss(self.canvas_width/2, sigma_x))
            y = int(random.gauss(self.canvas_height/2, sigma_y))
            # Ensure point is within margin of the canvas
            if not (self.canvas_width/4 <= x <= 3*self.canvas_width/4 and self.canvas_height/4 <= y <= 3*self.canvas_height/4):
                tries += 1
                continue
            candidate = (x, y)
            if all(math.dist(candidate, p) >= dmin for p in pts):
                pts.append(candidate)
                cols.append("#{:06x}".format(random.randint(0, 0xFFFFFF)))
            tries += 1
        self.points = pts
        self.colors = cols
        self.draw_points()

    def generate_exponential_points(self):
        """Generate a set of points using exponential randomness centered at the canvas center.
        
        The exponential distribution is shifted to be centered.
        """
        try:
            num_points = int(self.random_entry.get())
            if num_points <= 0:
                return
        except ValueError:
            return
        pts = []
        cols = []
        scale = self.canvas_width / 8
        dmin = self._min_distance(num_points)
        max_tries = num_points * 1000
        tries = 0
        while len(pts) < num_points and tries < max_tries:
            # Generate exponential value.
            exp_x = random.expovariate(1/scale)
            exp_y = random.expovariate(1/scale)
            # Assign random sign to x and y
            exp_x = exp_x if random.random() < 0.5 else -exp_x
            exp_y = exp_y if random.random() < 0.5 else -exp_y
            x = int(self.canvas_width/2 + exp_x)
            y = int(self.canvas_height/2 + exp_y)
            # Check bounds
            if not (self.canvas_width/4 <= x <= 3*self.canvas_width/4 and self.canvas_height/4 <= y <= 3*self.canvas_height/4):
                tries += 1
                continue
            candidate = (x, y)
            if all(math.dist(candidate, p) >= dmin for p in pts):
                pts.append(candidate)
                cols.append("#{:06x}".format(random.randint(0, 0xFFFFFF)))
            tries += 1
        self.points = pts
        self.colors = cols
        self.draw_points()

    # --- Helper Functions for SVG Writing ---
    def _write_svg_header(self, f):
        f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        f.write(f'<svg width="{self.canvas_width}" height="{self.canvas_height}" xmlns="http://www.w3.org/2000/svg">\n')
    
    def _write_svg_footer(self, f):
        f.write('</svg>\n')
    
    def _write_points(self, f, points, point_color_map=None, default_color="black"):
        """Write points as SVG circle elements.
        
        Args:
            f: File object.
            points: List of (x, y) tuples.
            point_color_map: Optional dict mapping point -> color.
            default_color: Default color for points.
        """
        for p in points:
            color = point_color_map.get(p, default_color) if point_color_map else default_color
            f.write(f'<circle cx="{p[0]}" cy="{p[1]}" r="3" fill="{color}" stroke="{color}" stroke-width="2"/>\n')
    
    def _write_hulls(self, f, hull_layers, hull_colors):
        """Write convex hull layers as SVG polygon elements.
        
        Args:
            f: File object.
            hull_layers: List of hull layers (each a list of (x, y) tuples).
            hull_colors: List of colors to use for hulls.
        """
        for i, hull in enumerate(hull_layers):
            color = hull_colors[i % len(hull_colors)]
            pts_str = " ".join(f"{p[0]},{p[1]}" for p in hull)
            f.write(f'<polygon points="{pts_str}" stroke="{color}" fill="none" stroke-width="2"/>\n')
    
    def _write_circles(self, f, circles, circle_point_map, circle_colors):
        """Write circles as SVG circle elements and their points in matching color.
        
        Args:
            f: File object.
            circles: List of (cx, cy, r) tuples.
            circle_point_map: List of lists of points assigned to each circle.
            circle_colors: List of colors for the circles.
        """
        for i, circle in enumerate(circles):
            color = circle_colors[i % len(circle_colors)]
            cx, cy, r = circle
            f.write(f'<circle cx="{int(cx)}" cy="{int(cy)}" r="{int(r)}" stroke="{color}" fill="none" stroke-width="2"/>\n')
            # Write the points belonging to this circle in the same color.
            pts = circle_point_map[i]
            for p in pts:
                f.write(f'<circle cx="{p[0]}" cy="{p[1]}" r="3" fill="{color}" stroke="{color}" stroke-width="2"/>\n')
    
    # --- Compute Layers ---
    def compute_hull_layers(self):
        """Compute convex hull layers (iterative removal) and return as a list of layers."""
        layers = []
        pts = self.points[:]
        while len(pts) >= 3:
            hull = Geometry.graham_scan(pts)
            if len(hull) < 3:
                break
            layers.append(hull)
            pts = [p for p in pts if p not in hull]
        if len(pts) == 2:
            layers.append(pts)
        return layers
    
    def compute_circle_layers(self):
        """Compute iterative minimum enclosing circles and return circles and their assigned points.
        
        Returns:
            tuple: (circles, circle_point_map) where circles is a list of (cx, cy, r) and
                   circle_point_map is a list of lists, each containing the points enclosed by that circle.
        """
        circles = []
        circle_point_map = []
        pts_remaining = self.points[:]
        while pts_remaining:
            circle = Geometry.MINIDISC(pts_remaining)
            if circle is None:
                break
            assigned = [p for p in pts_remaining if Geometry.is_point_on_disc(p, circle)]
            circles.append(circle)
            circle_point_map.append(assigned)
            pts_remaining = [p for p in pts_remaining if p not in assigned]
        return circles, circle_point_map

    def save_out_txt(self):
        """Save a text file 'out.txt' with hull layers and circle layers.
        
        The file will list:
          - #hull layers: Each hull layer and the points in that layer (outermost to innermost)
          - #circles: Each circle and the points on its boundary (biggest to smallest)
        """
        hull_layers = self.compute_hull_layers()
        circles, _ = self.compute_circle_layers()
        with open("out.txt", "w") as f:
            f.write("#hull layers\n")
            for i, layer in enumerate(hull_layers, start=1):
                f.write(f"Layer {i} ({len(layer)} points): " + " ".join(f"({p[0]},{p[1]})" for p in layer) + "\n")
            f.write("\n#circles\n")
            for i, circle in enumerate(circles, start=1):
                cx, cy, r = circle
                f.write(f"Circle {i} (center=({int(cx)},{int(cy)}), r={int(r)}):\n")
                # List points on the circle boundary
                pts = [p for p in self.points if Geometry.is_point_on_disc(p, circle)]
                f.write(" " + " ".join(f"({p[0]},{p[1]})" for p in pts) + "\n")
    
    # --- SVG Export Functions ---
    def save_points_svg(self):
        """Save an SVG file with the input points drawn in black."""
        with open("points.svg", "w") as f:
            self._write_svg_header(f)
            self._write_points(f, self.points, default_color="black")
            self._write_svg_footer(f)
    
    def save_hull_svg(self):
        """Save an SVG file with convex hulls (alternating red and blue) 
        overlaid with points colored according to the hull they belong to."""
        hull_layers = self.compute_hull_layers()
        hull_colors = ["red", "blue"]
        point_color_map = {}
        for i, hull in enumerate(hull_layers):
            color = hull_colors[i % len(hull_colors)]
            for p in hull:
                point_color_map[p] = color
        with open("hull.svg", "w") as f:
            self._write_svg_header(f)
            self._write_hulls(f, hull_layers, hull_colors)
            self._write_points(f, self.points, point_color_map=point_color_map, default_color="black")
            self._write_svg_footer(f)
    
    def save_cir_svg(self):
        """Save an SVG file with iterative minimum enclosing circles (alternating #AAFF00 and orange)
        overlaid with points colored according to the circle they are a part of."""
        circles, circle_point_map = self.compute_circle_layers()
        circle_colors = ["#AAFF00", "orange"]
        with open("cir.svg", "w") as f:
            self._write_svg_header(f)
            self._write_circles(f, circles, circle_point_map, circle_colors)
            self._write_svg_footer(f)
    
    def save_all_svg(self):
        """Save an SVG file that overlays hulls, then circles, then points."""
        hull_layers = self.compute_hull_layers()
        hull_colors = ["red", "blue"]
        circles, circle_point_map = self.compute_circle_layers()
        circle_colors = ["#AAFF00", "orange"]
        point_color_map = {}
        for i, hull in enumerate(hull_layers):
            color = hull_colors[i % len(hull_colors)]
            for p in hull:
                point_color_map[p] = color
        for i, pts in enumerate(circle_point_map):
            color = circle_colors[i % len(circle_colors)]
            for p in pts:
                point_color_map[p] = color
        with open("all.svg", "w") as f:
            self._write_svg_header(f)
            self._write_hulls(f, hull_layers, hull_colors)
            self._write_circles(f, circles, circle_point_map, circle_colors)
            self._write_points(f, self.points, point_color_map=point_color_map, default_color="black")
            self._write_svg_footer(f)
    
    def save_svg(self):
        """Save all four SVG files: points.svg, hull.svg, cir.svg, and all.svg, and save a summary in out.txt."""
        self.save_points_svg()
        self.save_hull_svg()
        self.save_cir_svg()
        self.save_all_svg()
        self.save_out_txt()
        print("SVG files saved: points.svg, hull.svg, cir.svg, all.svg, and summary out.txt")
    
    # --- Dragging Functions ---
    def start_drag(self, event, idx):
        """Initiate dragging of a point."""
        self.drag_data = {
            "x": event.x,
            "y": event.y,
            "idx": idx,
            "initial_point": self.points[idx]
        }
    
    def drag(self, event, idx):
        """Update the position of a point while dragging."""
        if self.drag_data["idx"] != idx:
            return
        dx = event.x - self.drag_data["x"]
        dy = event.y - self.drag_data["y"]
        new_x = self.drag_data["initial_point"][0] + dx
        new_y = self.drag_data["initial_point"][1] + dy
        self.canvas.coords(self.oval_ids[idx], new_x - 5, new_y - 5, new_x + 5, new_y + 5)
    
    def end_drag(self, event):
        """Finalize dragging and update the point's position."""
        if self.drag_data["idx"] is None:
            return
        idx = self.drag_data["idx"]
        dx = event.x - self.drag_data["x"]
        dy = event.y - self.drag_data["y"]
        self.points[idx] = (
            self.drag_data["initial_point"][0] + dx,
            self.drag_data["initial_point"][1] + dy
        )
        self.draw_points()
        self.drag_data["idx"] = None
    
    # --- Point Loading, Generation, and Drawing ---
    def load_points(self):
        """Load points from a text file.
        
        File format:
            First line: number of points (n)
            Next n lines: x y (space-separated integers)
        """
        file_path = filedialog.askopenfilename(filetypes=[("Text Files", "*.txt")])
        if not file_path:
            return
        with open(file_path, "r") as file:
            lines = file.readlines()
            num_points = int(lines[0].strip())
            self.points = []
            self.colors = []
            for line in lines[1:num_points + 1]:
                x, y = map(int, line.strip().split())
                self.points.append((x, y))
                self.colors.append("#{:06x}".format(random.randint(0, 0xFFFFFF)))
        self.draw_points()
    
    def generate_random_points(self):
        """Alias for uniform generation."""
        self.generate_uniform_points()
    
    def resize_canvas(self):
        """Resize the canvas and scale all points accordingly."""
        try:
            new_width = int(self.width_entry.get())
            new_height = int(self.height_entry.get())
        except ValueError:
            return
        if new_width <= 0 or new_height <= 0:
            return
        scale_x = new_width / self.canvas_width
        scale_y = new_height / self.canvas_height
        self.canvas_width = new_width
        self.canvas_height = new_height
        self.canvas.config(width=self.canvas_width, height=self.canvas_height)
        self.points = [(int(x * scale_x), int(y * scale_y)) for (x, y) in self.points]
        self.draw_points()
    
    def draw_points(self):
        """Clear and redraw all points on the canvas."""
        self.canvas.delete("all")
        self.oval_ids = []
        for idx, ((x, y), color) in enumerate(zip(self.points, self.colors)):
            oval = self.canvas.create_oval(x - 5, y - 5, x + 5, y + 5,
                                           fill=color, outline=color, tags=f"point_{idx}")
            self.oval_ids.append(oval)
            self.canvas.tag_bind(oval, "<Button-1>", lambda e, idx=idx: self.start_drag(e, idx))
            self.canvas.tag_bind(oval, "<B1-Motion>", lambda e, idx=idx: self.drag(e, idx))
    
    # --- Convex Hull Functions ---
    def draw_convex_hulls(self):
        """Draw convex hulls for the set of points."""
        hull_layers = self.compute_hull_layers()
        for i, hull in enumerate(hull_layers):
            color = f"#{random.randint(0, 0xFFFFFF):06x}"
            self.canvas.create_polygon(hull, outline=color, fill="", width=2)
    
    # --- Minimum Enclosing Circle Functions ---
    def draw_enclosing_circles(self):
        """Draw minimum enclosing circles for all the points."""
        circles, _ = self.compute_circle_layers()
        for cx, cy, r in circles:
            color = f"#{random.randint(0, 0xFFFFFF):06x}"
            self.canvas.create_oval(cx - r, cy - r, cx + r, cy + r, outline=color, fill="", width=2)

if __name__ == "__main__":
    root = tk.Tk()
    app = PointDragger(root, width=600, height=800)
    root.mainloop()