import tkinter as tk
from tkinter import ttk, messagebox
from math_engine import BeamSolver
from visualization import plot_diagrams

class BeamApp:
    def __init__(self, root):
        self.root = root
        self.root.title("CE 257 Beam Solver (Three-Moment Equation)")
        self.root.geometry("850x650")
        
        self.spans = []
        self.loads = {} # dict of list of loads {span_idx: [{'type': 'uniform', 'w': w}, ...]}
        
        self.create_widgets()
        self.draw_grid()

    def create_widgets(self):
        # Top Frame for Controls
        control_frame = ttk.LabelFrame(self.root, text="Beam Setup", padding="10")
        control_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=5)
        
        # Add Span
        ttk.Label(control_frame, text="Span Length (m):").grid(row=0, column=0, padx=5, pady=5)
        self.span_length_var = tk.StringVar()
        ttk.Entry(control_frame, textvariable=self.span_length_var, width=10).grid(row=0, column=1, padx=5, pady=5)
        ttk.Button(control_frame, text="Add Span", command=self.add_span).grid(row=0, column=2, padx=5, pady=5)
        
        ttk.Separator(control_frame, orient=tk.HORIZONTAL).grid(row=1, column=0, columnspan=9, sticky="ew", pady=10)
        
        # Add Load
        ttk.Label(control_frame, text="Load Type:").grid(row=2, column=0, padx=5, pady=5)
        self.load_type_var = tk.StringVar(value="uniform")
        ttk.Combobox(control_frame, textvariable=self.load_type_var, values=["uniform", "point"], width=10, state="readonly").grid(row=2, column=1, padx=5, pady=5)
        
        ttk.Label(control_frame, text="Span Index (0 to n-1):").grid(row=2, column=2, padx=5, pady=5)
        self.load_span_idx_var = tk.StringVar()
        ttk.Entry(control_frame, textvariable=self.load_span_idx_var, width=5).grid(row=2, column=3, padx=5, pady=5)
        
        ttk.Label(control_frame, text="Magnitude (kN or kN/m):").grid(row=2, column=4, padx=5, pady=5)
        self.load_mag_var = tk.StringVar()
        ttk.Entry(control_frame, textvariable=self.load_mag_var, width=8).grid(row=2, column=5, padx=5, pady=5)
        
        # Add Load - Second Row (so it fits on screen)
        ttk.Label(control_frame, text="Dist 'a' (m) (point load only):").grid(row=3, column=4, padx=5, pady=5)
        self.load_dist_var = tk.StringVar()
        ttk.Entry(control_frame, textvariable=self.load_dist_var, width=8).grid(row=3, column=5, padx=5, pady=5)
        
        ttk.Button(control_frame, text="Add Load", command=self.add_load).grid(row=3, column=6, padx=5, pady=5)
        
        # Action Buttons
        action_frame = ttk.Frame(self.root, padding="10")
        action_frame.pack(side=tk.TOP, fill=tk.X, padx=10)
        
        ttk.Button(action_frame, text="Clear All", command=self.clear_all).pack(side=tk.LEFT, padx=5)
        ttk.Button(action_frame, text="Solve & Plot BMD/SFD", command=self.solve_and_plot).pack(side=tk.RIGHT, padx=5)
        
        # Canvas for Drawing Beam
        canvas_frame = ttk.LabelFrame(self.root, text="Beam Visualization Grid", padding="10")
        canvas_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        self.canvas = tk.Canvas(canvas_frame, bg="white")
        self.canvas.pack(fill=tk.BOTH, expand=True)
        
        # Info text
        self.info_text = tk.Text(self.root, height=6, width=80)
        self.info_text.pack(side=tk.BOTTOM, fill=tk.X, padx=10, pady=10)
        self.info_text.insert(tk.END, "Status: Ready\n")
        self.info_text.config(state=tk.DISABLED)

    def log(self, message):
        self.info_text.config(state=tk.NORMAL)
        self.info_text.insert(tk.END, message + "\n")
        self.info_text.see(tk.END)
        self.info_text.config(state=tk.DISABLED)

    def draw_grid(self):
        self.canvas.delete("all")
        self.canvas.update()
        w = self.canvas.winfo_width()
        h = self.canvas.winfo_height()
        
        if w == 1: # Canvas not fully rendered yet
            w, h = 800, 300
        
        # Draw grid lines
        for i in range(0, w, 20):
            self.canvas.create_line(i, 0, i, h, fill="#f0f0f0")
        for i in range(0, h, 20):
            self.canvas.create_line(0, i, w, i, fill="#f0f0f0")
            
        # Draw base beam line (visual representation only, scaled)
        if not self.spans:
            self.canvas.create_text(w/2, h/2, text="Enter a span length to begin drawing.", fill="gray", font=("Arial", 14))
            return
            
        total_length = sum(self.spans)
        if total_length == 0:
            return
            
        scale = (w - 100) / total_length
        start_x = 50
        y = h / 2
        
        self.canvas.create_line(start_x, y, start_x + total_length * scale, y, width=6, fill="#3498db")
        
        # Draw supports
        current_x = start_x
        self.draw_support(current_x, y)
        for span in self.spans:
            current_x += span * scale
            self.draw_support(current_x, y)
            
        # Draw loads
        for span_idx, load_list in self.loads.items():
            span_start_x = start_x + sum(self.spans[:span_idx]) * scale
            span_len_px = self.spans[span_idx] * scale
            
            for load in load_list:
                if load['type'] == 'uniform':
                    # Draw series of arrows
                    num_arrows = max(5, int(span_len_px / 20))
                    for i in range(num_arrows + 1):
                        lx = span_start_x + (i / num_arrows) * span_len_px
                        self.canvas.create_line(lx, y - 40, lx, y, arrow=tk.LAST, fill="#e74c3c", width=2)
                    self.canvas.create_text(span_start_x + span_len_px/2, y - 50, text=f"{load['w']} kN/m", fill="#c0392b", font=("Arial", 10, "bold"))
                elif load['type'] == 'point':
                    lx = span_start_x + load['a'] * scale
                    self.canvas.create_line(lx, y - 60, lx, y, arrow=tk.LAST, fill="#e74c3c", width=3)
                    self.canvas.create_text(lx, y - 70, text=f"{load['P']} kN", fill="#c0392b", font=("Arial", 10, "bold"))

    def draw_support(self, x, y):
        self.canvas.create_polygon(x, y, x - 10, y + 20, x + 10, y + 20, fill="#2ecc71", outline="black")
        self.canvas.create_line(x - 15, y + 20, x + 15, y + 20, width=3, fill="black")

    def add_span(self):
        try:
            length = float(self.span_length_var.get())
            if length <= 0:
                raise ValueError
            self.spans.append(length)
            self.log(f"Added span of length {length}m. Total spans: {len(self.spans)}")
            self.draw_grid()
            self.span_length_var.set("")
        except ValueError:
            messagebox.showerror("Error", "Please enter a valid positive number for span length.")

    def add_load(self):
        try:
            span_idx = int(self.load_span_idx_var.get())
            if span_idx < 0 or span_idx >= len(self.spans):
                messagebox.showerror("Error", f"Invalid span index. Must be between 0 and {len(self.spans)-1}.")
                return
                
            l_type = self.load_type_var.get()
            mag = float(self.load_mag_var.get())
            
            load = {'type': l_type}
            if l_type == 'uniform':
                load['w'] = mag
                self.log(f"Added uniform load {mag} kN/m to span {span_idx}.")
            else:
                dist_str = self.load_dist_var.get()
                if not dist_str:
                    messagebox.showerror("Error", "Please enter distance 'a' for the point load.")
                    return
                dist = float(dist_str)
                if dist < 0 or dist > self.spans[span_idx]:
                    messagebox.showerror("Error", "Distance 'a' must be within the span length.")
                    return
                load['P'] = mag
                load['a'] = dist
                self.log(f"Added point load {mag} kN at {dist}m to span {span_idx}.")
                
            if span_idx not in self.loads:
                self.loads[span_idx] = []
            self.loads[span_idx].append(load)
            
            self.draw_grid()
            self.load_mag_var.set("")
            self.load_dist_var.set("")
        except ValueError:
            messagebox.showerror("Error", "Please enter valid numbers for load parameters.")

    def clear_all(self):
        self.spans = []
        self.loads = {}
        self.log("Cleared all spans and loads.")
        self.draw_grid()

    def solve_and_plot(self):
        if not self.spans:
            messagebox.showwarning("Warning", "Please add at least one span before solving.")
            return
            
        try:
            self.log("Solving using Three-Moment Equation...")
            solver = BeamSolver(self.spans, self.loads)
            moments = solver.solve_moments()
            
            self.log("Moments at supports:")
            for i, m in enumerate(moments):
                self.log(f"  Support {i} (M{i}): {m:.2f} kNm")
                
            x, V, M = solver.compute_internal_forces()
            self.log("Plotting BMD and SFD...")
            plot_diagrams(x, V, M, self.spans)
            self.log("Successfully plotted diagrams.")
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred during solution:\n{str(e)}")

def start_gui():
    root = tk.Tk()
    app = BeamApp(root)
    # Give canvas time to size before drawing
    root.after(100, app.draw_grid)
    root.mainloop()

if __name__ == "__main__":
    start_gui()
