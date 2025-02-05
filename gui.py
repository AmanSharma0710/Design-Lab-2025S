import tkinter as tk
from tkinter import filedialog
import random

class PointDragger:
    def __init__(self, root, width=600, height=800):
        self.root = root
        self.root.title("Point Dragger")
        self.canvas_width = width
        self.canvas_height = height
        
        self.controls_frame = tk.Frame(root)
        self.controls_frame.pack(side=tk.TOP, fill=tk.X, pady=5)

        self.file_frame = tk.Frame(self.controls_frame)
        self.file_frame.pack(fill=tk.X, pady=2)
        
        self.load_button = tk.Button(self.file_frame, text="Load Points from File", command=self.load_points)
        self.load_button.pack(side=tk.LEFT, padx=5)

        self.random_label = tk.Label(self.file_frame, text="Number of Points:")
        self.random_label.pack(side=tk.LEFT, padx=5)

        self.random_entry = tk.Entry(self.file_frame, width=5)
        self.random_entry.pack(side=tk.LEFT)

        self.random_button = tk.Button(self.file_frame, text="Generate", command=self.generate_random_points)
        self.random_button.pack(side=tk.LEFT, padx=5)

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
        
        self.canvas = tk.Canvas(root, width=self.canvas_width, height=self.canvas_height, bg="white")
        self.canvas.pack()

        self.hull_button = tk.Button(root, text="Calculate Convex Hull", command=self.calculate_convex_hulls)
        self.hull_button.pack(pady=5)

        self.points = []
        self.colors = []
        self.oval_ids = []
        self.drag_data = {"x": 0, "y": 0, "idx": None, "initial_point": (0, 0)}

        self.canvas.bind("<ButtonRelease-1>", self.end_drag)
    
    def start_drag(self, event, idx):
        self.drag_data = {
            "x": event.x,
            "y": event.y,
            "idx": idx,
            "initial_point": self.points[idx]
        }
    
    def drag(self, event, idx):
        if self.drag_data["idx"] != idx:
            return
        dx = event.x - self.drag_data["x"]
        dy = event.y - self.drag_data["y"]
        new_x = self.drag_data["initial_point"][0] + dx
        new_y = self.drag_data["initial_point"][1] + dy
        self.canvas.coords(self.oval_ids[idx], new_x - 5, new_y - 5, new_x + 5, new_y + 5)
    
    def end_drag(self, event):
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

    def load_points(self):
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
        self.points = [(int(x * scale_x), int(y * scale_y)) for x, y in self.points]
        self.draw_points()

    def draw_points(self):
        self.canvas.delete("all")
        self.oval_ids = []
        for idx, ((x, y), color) in enumerate(zip(self.points, self.colors)):
            oval = self.canvas.create_oval(x - 5, y - 5, x + 5, y + 5, fill=color, outline="black", tags=f"point_{idx}")
            self.oval_ids.append(oval)
            self.canvas.tag_bind(oval, "<Button-1>", lambda e, idx=idx: self.start_drag(e, idx))
            self.canvas.tag_bind(oval, "<B1-Motion>", lambda e, idx=idx: self.drag(e, idx))

    def calculate_convex_hulls(self):
        points = self.points[:]
        while len(points) > 2:
            hull = self.graham_scan(points)
            if len(hull) < 3:
                break
            polygon_color = "#{:06x}".format(random.randint(0, 0xFFFFFF))
            self.canvas.create_polygon(hull, outline=polygon_color, fill="", width=2)
            points = [p for p in points if p not in hull]

    def graham_scan(self, points):
        points = sorted(points, key=lambda p: (p[1], p[0]))
        def cross_product(o, a, b):
            return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
        lower, upper = [], []
        for p in points:
            while len(lower) >= 2 and cross_product(lower[-2], lower[-1], p) <= 0:
                lower.pop()
            lower.append(p)
        for p in reversed(points):
            while len(upper) >= 2 and cross_product(upper[-2], upper[-1], p) <= 0:
                upper.pop()
            upper.append(p)
        return lower[:-1] + upper[:-1]

if __name__ == "__main__":
    root = tk.Tk()
    app = PointDragger(root, width=600, height=800)
    root.mainloop()