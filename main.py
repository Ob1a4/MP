import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.animation import FuncAnimation
import tkinter as tk
from tkinter import ttk, messagebox
from matplotlib.ticker import FormatStrFormatter
import math

class HeatEquationApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Решение уравнения теплопроводности")
        self.animation = None
        self.is_paused = False
        self.current_frame = 0
        self.u = None
        self.x = None
        self.simulation_complete = False
        self.setup_ui()
        
    def setup_ui(self):
        ttk.Label(self.root, text="Длина стержня, L:").grid(row=0, column=0, padx=5, pady=5, sticky="e")
        self.L_entry = ttk.Entry(self.root)
        self.L_entry.grid(row=0, column=1, padx=5, pady=5)
        self.L_entry.insert(0, "1.0")
        
        ttk.Label(self.root, text="Количество узлов, N:").grid(row=1, column=0, padx=5, pady=5, sticky="e")
        self.N_entry = ttk.Entry(self.root)
        self.N_entry.grid(row=1, column=1, padx=5, pady=5)
        self.N_entry.insert(0, "50")
        
        ttk.Label(self.root, text="Коэффициент температуропроводности, a:").grid(row=2, column=0, padx=5, pady=5, sticky="e")
        self.a_entry = ttk.Entry(self.root)
        self.a_entry.grid(row=2, column=1, padx=5, pady=5)
        self.a_entry.insert(0, "1")
        
        ttk.Label(self.root, text="Начальное условие, φ(x):").grid(row=3, column=0, padx=5, pady=5, sticky="e")
        self.init_entry = ttk.Entry(self.root)
        self.init_entry.grid(row=3, column=1, padx=5, pady=5)
        self.init_entry.insert(0, "50*sin(pi*x)")

        ttk.Label(self.root, text="Температура на левом конце, T_0:").grid(row=4, column=0, padx=5, pady=5, sticky="e")
        self.left_temp_entry = ttk.Entry(self.root)
        self.left_temp_entry.grid(row=4, column=1, padx=5, pady=5)
        self.left_temp_entry.insert(0, "0")

        ttk.Label(self.root, text="Температура на правом конце, T_L:").grid(row=5, column=0, padx=5, pady=5, sticky="e")
        self.right_temp_entry = ttk.Entry(self.root)
        self.right_temp_entry.grid(row=5, column=1, padx=5, pady=5)
        self.right_temp_entry.insert(0, "0")
        
        ttk.Label(self.root, text="Шаг по времени (с):").grid(row=6, column=0, padx=5, pady=5, sticky="e")
        self.dt_entry = ttk.Entry(self.root)
        self.dt_entry.grid(row=6, column=1, padx=5, pady=5)
        self.dt_entry.insert(0, "0.1")
        
        ttk.Label(self.root, text="Общее время, T (с):").grid(row=7, column=0, padx=5, pady=5, sticky="e")
        self.total_time_entry = ttk.Entry(self.root)
        self.total_time_entry.grid(row=7, column=1, padx=5, pady=5)
        self.total_time_entry.insert(0, "10")
        
        self.run_button = ttk.Button(self.root, text="Старт", command=self.start_simulation)
        self.run_button.grid(row=8, column=0, padx=5, pady=5)
        
        self.pause_button = ttk.Button(self.root, text="Пауза", state="disabled", command=self.toggle_pause)
        self.pause_button.grid(row=8, column=1, padx=5, pady=5)
        
        self.figure = plt.Figure(figsize=(8, 5), dpi=100)
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.root)
        self.canvas.get_tk_widget().grid(row=0, column=2, rowspan=9, padx=10, pady=10)
        
        self.crosshair_v = self.ax.axvline(color='blue', alpha=0.5, lw=1, visible=False)
        self.crosshair_h = self.ax.axhline(color='blue', alpha=0.5, lw=1, visible=False)
        self.coord_label = self.ax.text(0.02, 1.02, '', transform=self.ax.transAxes, 
                                      bbox=dict(facecolor='white', alpha=0.8))
        
        self.canvas.mpl_connect("motion_notify_event", self.on_mouse_move)
        self.canvas.mpl_connect("axes_leave_event", self.on_axes_leave)
        
    def start_simulation(self):
        try:
            if self.animation is not None:
                if self.animation.event_source is not None:
                    self.animation.event_source.stop()
                self.animation = None
            
            self.run_button.config(state="disabled")
            self.simulation_complete = False
            
            L = float(self.L_entry.get())
            N_total = int(self.N_entry.get())
            a = float(self.a_entry.get())
            dt = float(self.dt_entry.get())
            total_time = float(self.total_time_entry.get())
            left_temp = float(self.left_temp_entry.get())
            right_temp = float(self.right_temp_entry.get())
            
            self.left_temp = left_temp
            self.right_temp = right_temp
            
            init_str = self.init_entry.get()
            init_str = init_str.replace('^', '**')

            x = np.linspace(0, L, N_total)
            
            math_functions = {
                'sin': np.sin, 'cos': np.cos, 'tan': np.tan,
                'exp': np.exp, 'log': np.log, 'log10': np.log10,
                'sqrt': np.sqrt, 'pi': np.pi, 'e': np.e
            }
            math_functions.update({'x': x, 'L': L})
            
            try:
                initial_temp = float(init_str)
                u = np.full(N_total, initial_temp)
            except ValueError:
                u = eval(init_str, {'__builtins__': None}, math_functions)
                if isinstance(u, (int, float)):
                    u = np.full(N_total, u)
            
            u[0] = left_temp
            u[-1] = right_temp
            
            h = L / (N_total - 1)
            r = a * dt / h**2
            
            self.ax.clear()
            self.line, = self.ax.plot(x, u, 'r-', linewidth=1)
            self.points, = self.ax.plot(x, u, 'ro', markersize=4)
            self.ax.set_xlim(0, L)
            self.ax.set_ylim(min(0, np.min(u))-0.1, max(0, np.max(u))+0.1)
            self.ax.set_xlabel('Положение вдоль стержня (м)')
            self.ax.set_ylabel('Температура (°C)')
            self.ax.grid(True)
            self.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            self.time_text = self.ax.text(0.02, 0.95, '', transform=self.ax.transAxes)
            
            self.u = u
            self.x = x
            self.N = N_total - 2
            self.dt = dt
            self.total_time = total_time
            self.r = r
            self.all_frames = [u.copy()]
            self.current_frame = 0
            self.is_paused = False
            self.pause_button.config(text="Пауза", state="normal")

            if total_time <= 0:
                self.canvas.draw()
                self.run_button.config(state="normal")
                self.pause_button.config(state="disabled")
                self.simulation_complete = True
            else:
                self.start_animation()
            
        except Exception as e:
            messagebox.showerror("Ошибка", f"Произошла ошибка:\n{str(e)}")
            self.run_button.config(state="normal")
    
    def start_animation(self):
        frames = int(self.total_time / self.dt)
        self.animation = FuncAnimation(
            self.figure, 
            self.update_plot, 
            frames=frames,
            init_func=lambda: [self.line, self.points, self.time_text],
            interval=50, 
            blit=False,
            repeat=False,
            save_count=frames
        )
        self.canvas.draw()
    
    def toggle_pause(self):
        if self.animation:
            if self.is_paused:
                self.animation.event_source.start()
                self.is_paused = False
                self.pause_button.config(text="Пауза")
            else:
                self.animation.event_source.stop()
                self.is_paused = True
                self.pause_button.config(text="Продолжить")
    
    def update_plot(self, frame):
        if not self.is_paused:
            n = self.N
            alpha = np.zeros(n)
            beta = np.zeros(n)
            u_internal = self.u[1:-1].copy()

            a = np.full(n, -self.r)
            b = np.full(n, 1 + 2*self.r)
            c = np.full(n, -self.r)
            d = u_internal.copy()

            d[0] += self.r * self.left_temp
            d[-1] += self.r * self.right_temp

            alpha[0] = -c[0] / b[0]
            beta[0] = d[0] / b[0]
            for i in range(1, n):
                denominator = b[i] + a[i] * alpha[i-1]
                alpha[i] = -c[i] / denominator
                beta[i] = (d[i] - a[i] * beta[i-1]) / denominator

            u_internal[n-1] = beta[n-1]
            for i in range(n-2, -1, -1):
                u_internal[i] = alpha[i] * u_internal[i+1] + beta[i]

            self.u[1:-1] = u_internal

            self.current_frame = frame
            self.all_frames.append(self.u.copy())

        self.line.set_ydata(self.u)
        self.points.set_ydata(self.u)
        current_time = (self.current_frame+1)*self.dt
        self.time_text.set_text(f'Время: {current_time:.2f} с')

        if self.current_frame >= int(self.total_time / self.dt) - 1:
            self.run_button.config(state="normal")
            self.pause_button.config(state="disabled")
            self.simulation_complete = True

        return [self.line, self.points, self.time_text, self.crosshair_v, self.crosshair_h, self.coord_label]

    
    def on_mouse_move(self, event):
        if event.inaxes == self.ax and self.x is not None:
            self.crosshair_v.set_xdata([event.xdata, event.xdata])
            self.crosshair_h.set_ydata([event.ydata, event.ydata])
            self.crosshair_v.set_visible(True)
            self.crosshair_h.set_visible(True)
            
            idx = np.argmin(np.abs(self.x - event.xdata))
            x_val = self.x[idx]
            
            if self.simulation_complete and len(self.all_frames) > 0:
                y_val = self.all_frames[-1][idx]
            elif self.u is not None:
                y_val = self.u[idx]
            else:
                y_val = 0
            
            self.coord_label.set_text(f'x={x_val:.3f} м, T={y_val:.3f} °C')
            self.canvas.draw_idle()
    
    def on_axes_leave(self, event):
        self.crosshair_v.set_visible(False)
        self.crosshair_h.set_visible(False)
        self.coord_label.set_text('')
        self.canvas.draw_idle()

if __name__ == "__main__":
    root = tk.Tk()
    app = HeatEquationApp(root)
    root.mainloop()
