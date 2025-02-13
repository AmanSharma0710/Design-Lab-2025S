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
    """GUI application to visualize point dragging, convex hulls, and minimum enclosing circles."""
    def __init__(self, root, width=600, height=800):
        self.root = root
        self.root.title("Point Dragger")
        self.canvas_width = width
        self.canvas_height = height
        
        # --- Controls ---
        self.controls_frame = tk.Frame(root)
        self.controls_frame.pack(side=tk.TOP, fill=tk.X, pady=5)
        
        # File and random points controls
        self.file_frame = tk.Frame(self.controls_frame)
        self.file_frame.pack(fill=tk.X, pady=2)
        self.load_button = tk.Button(self.file_frame, text="Load Points from File", command=self.load_points)
        self.load_button.pack(side=tk.LEFT, padx=5)
        self.random_label = tk.Label(self.file_frame, text="Number of Points:")
        self.random_label.pack(side=tk.LEFT, padx=5)
        self.random_entry = tk.Entry(self.file_frame, width=5)
        self.random_entry.pack(side=tk.LEFT)
        # Renamed button text:
        self.random_button = tk.Button(self.file_frame, text="Generate Random Points", command=self.generate_random_points)
        self.random_button.pack(side=tk.LEFT, padx=5)
        
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
        # Renamed button text:
        self.hull_button = tk.Button(root, text="Generate Convex Hulls", command=self.calculate_convex_hulls)
        self.hull_button.pack(pady=5)
        self.circle_button = tk.Button(root, text="Generate Enclosing Circles", command=self.generate_enclosing_circles)
        self.circle_button.pack(pady=5)
        
        # --- Data ---
        self.points = []     # List of (x, y) tuples
        self.colors = []     # List of colors for points
        self.oval_ids = []   # Canvas IDs for points
        self.drag_data = {"x": 0, "y": 0, "idx": None, "initial_point": (0, 0)}
        
        self.canvas.bind("<ButtonRelease-1>", self.end_drag)
    
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
        """Generate a random set of points within the canvas dimensions."""
        try:
            num_points = int(self.random_entry.get())
            if num_points <= 0:
                return
        except ValueError:
            return
        self.points = []
        self.colors = []
        for _ in range(num_points):
            x = random.randint(0, self.canvas_width)
            y = random.randint(0, self.canvas_height)
            self.points.append((x, y))
            self.colors.append("#{:06x}".format(random.randint(0, 0xFFFFFF)))
        self.draw_points()
    
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
    def calculate_convex_hulls(self):
        """Generate and draw convex hull layers using the Graham scan algorithm."""
        pts = self.points[:]
        while len(pts) > 2:
            hull = Geometry.graham_scan(pts)
            if len(hull) < 3:
                break
            polygon_color = "#{:06x}".format(random.randint(0, 0xFFFFFF))
            self.canvas.create_polygon(hull, outline=polygon_color, fill="", width=2)
            pts = [p for p in pts if p not in hull]
        if len(pts) == 2:
            p1, p2 = pts
            color = "#{:06x}".format(random.randint(0, 0xFFFFFF))
            self.canvas.create_line(p1[0], p1[1], p2[0], p2[1], fill=color, width=2)
    
    # --- Minimum Enclosing Circle Functions ---
    def generate_enclosing_circles(self):
        """Iteratively compute and draw minimum enclosing circles for remaining points."""
        pts_remaining = self.points[:]
        while pts_remaining:
            circle = Geometry.MINIDISC(pts_remaining)
            if circle is None:
                break
            cx, cy, r = circle
            circle_color = "#{:06x}".format(random.randint(0, 0xFFFFFF))
            self.canvas.create_oval(cx - r, cy - r, cx + r, cy + r, outline=circle_color, width=2)
            # Remove points that are on or inside the circle.
            pts_remaining = [p for p in pts_remaining if not Geometry.is_point_on_disc(p, circle)]

if __name__ == "__main__":
    root = tk.Tk()
    app = PointDragger(root, width=600, height=800)
    root.mainloop()