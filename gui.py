import tkinter as tk
from tkinter import filedialog
import random

class PointDragger:
    def __init__(self, root):
        self.root = root
        self.root.title("Point Dragger")
        self.canvas = tk.Canvas(root, width=400, height=400, bg="white")
        self.canvas.pack()

        self.points = []  # Stores (x, y) coordinates in canvas scale
        self.colors = []
        self.oval_ids = []  # Track all point oval IDs
        self.drag_data = {"x": 0, "y": 0, "idx": None, "initial_point": (0, 0)}

        # Bind mouse release to entire canvas
        self.canvas.bind("<ButtonRelease-1>", self.end_drag)

        self.load_button = tk.Button(root, text="Load Points from File", command=self.load_points)
        self.load_button.pack()

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
                self.points.append((x * 2, y * 2))  # Scale to canvas size
                self.colors.append("#{:06x}".format(random.randint(0, 0xFFFFFF)))

        self.draw_points()

    def draw_points(self):
        self.canvas.delete("all")
        self.oval_ids = []
        for idx, ((x, y), color) in enumerate(zip(self.points, self.colors)):
            oval = self.canvas.create_oval(
                x - 5, y - 5, x + 5, y + 5,
                fill=color, outline="black", tags=f"point_{idx}"
            )
            self.oval_ids.append(oval)
            self.canvas.tag_bind(oval, "<Button-1>", lambda e, idx=idx: self.start_drag(e, idx))
            self.canvas.tag_bind(oval, "<B1-Motion>", lambda e, idx=idx: self.drag(e, idx))

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

        # Calculate movement from start position
        dx = event.x - self.drag_data["x"]
        dy = event.y - self.drag_data["y"]
        
        # Update visual position without modifying data
        new_x = self.drag_data["initial_point"][0] + dx
        new_y = self.drag_data["initial_point"][1] + dy
        
        # Move the oval directly on canvas
        self.canvas.coords(
            self.oval_ids[idx],
            new_x - 5, new_y - 5,
            new_x + 5, new_y + 5
        )

    def end_drag(self, event):
        if self.drag_data["idx"] is None:
            return

        # Calculate final position
        idx = self.drag_data["idx"]
        dx = event.x - self.drag_data["x"]
        dy = event.y - self.drag_data["y"]
        
        # Update data structure with new position
        self.points[idx] = (
            self.drag_data["initial_point"][0] + dx,
            self.drag_data["initial_point"][1] + dy
        )
        
        # Redraw all points to update positions and bindings
        self.draw_points()
        self.drag_data["idx"] = None  # Reset drag state

if __name__ == "__main__":
    root = tk.Tk()
    app = PointDragger(root)
    root.mainloop()