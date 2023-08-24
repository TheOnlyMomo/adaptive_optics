import tkinter as tk
from tkinter import ttk
import tkinter.messagebox as messagebox
import inspect
import importlib
from gaia_import_for_gui import *

class PageManager:
    def __init__(self, window):
        self.window = window
        self.min_mag_entry = None
        self.max_mag_entry = None
        self.telescope = None
        self.sci_camera = None
        self.wfs_camera = None
        self.exp_time_entry = None
        self.wfs = None
        self.save_plots_var = None
        self.log_yaxis_shwf_var = None
        self.save_plots_shwf_var = None
        self.shwf_exp_time_entry = None
        self.log_yaxis_var = None
        self.save_data_var = None
        self.l_poi_entry = None
        self.b_poi_entry = None
        self.num_searches_entry  = None

    def set_telescope(self, telescope_class):
        self.telescope = telescope_class

    def set_sci_camera(self, sci_camera_class):
        self.sci_camera = sci_camera_class

    def set_wfs_camera(self, wfs_camera_class):
        self.wfs_camera = wfs_camera_class

    def set_wfs(self, wfs_class):
        self.wfs = wfs_class


def get_classes(file_name):
    # Import the module
    module = importlib.import_module(file_name)

    # Create a dictionary of telescope names and their corresponding classes
    sci_camera_classes = [member for name, member in inspect.getmembers(module) if isinstance(member, module.sci_cameras)]
    wfs_camera_classes = [member for name, member in inspect.getmembers(module) if isinstance(member, module.wfs_cameras)]
    telescope_classes = {member for name, member  in inspect.getmembers(module) if isinstance(member, module.telescopes)}
    wfs_classes = [member for name, member in inspect.getmembers(module) if isinstance(member, module.shwf_sensor)]

    return sci_camera_classes,wfs_camera_classes,wfs_classes,telescope_classes




def open_next_page(self):
    
    # delete old widgets
    for widget in window.winfo_children():
        widget.destroy()

    sci_camera_classes,wfs_camera_classes,wfs_classes,telescope_classes = get_classes("gui_classes")

    selected_telescope = tk.StringVar()

    universal_frame = tk.Frame(window, bd=1, relief="solid")
    universal_frame.grid(sticky='nsew', padx=5, pady=5)

    tk.Label(universal_frame, text="Universal Parameters").grid(row=0, column=0, padx=10, pady=10, sticky='w')

    tk.Label(universal_frame, text="Min Mag").grid(row=1, column=0, padx=10, pady=5, sticky='w')
    self.min_mag_entry = tk.Entry(universal_frame)
    self.min_mag_entry.grid(row=1, column=1, padx=5, pady=5, sticky='ew')

    tk.Label(universal_frame, text="Max Mag").grid(row=2, column=0, padx=10, pady=5, sticky='w')
    self.max_mag_entry = tk.Entry(universal_frame)
    self.max_mag_entry.grid(row=2, column=1, padx=5, pady=5, sticky='ew')

    tk.Label(universal_frame, text="Telescope").grid(row=3, column=0, padx=10, pady=5, sticky='w')
    telescope_combobox = ttk.Combobox(universal_frame, textvariable=selected_telescope,state="readonly")
    telescope_combobox['values'] = [telescope.name for telescope in telescope_classes]
    telescope_combobox.grid(row=3, column=1, padx=5, pady=5, sticky='ew')

    def on_telescope_selected(event):
        selected_telescope_name = selected_telescope.get()

        for telescope_class in telescope_classes:
            if telescope_class.name == selected_telescope_name:
                 PageManager.set_telescope(self, telescope_class)

    def on_sci_camera_selected(event):
        selected_sci_camera_name = selected_sci_camera.get()

        for sci_camera_class in sci_camera_classes:
            if sci_camera_class.name == selected_sci_camera_name:
                 PageManager.set_sci_camera(self, sci_camera_class)

    def on_wfs_camera_selected(event):
        selected_wfs_camera_name = selected_wfs_camera.get()

        for wfs_camera_class in wfs_camera_classes:
            if wfs_camera_class.name == selected_wfs_camera_name:
                 PageManager.set_wfs_camera(self, wfs_camera_class)

    def on_wfs_selected(event):
        selected_wfs_name = selected_wfs.get()

        for wfs_class in wfs_classes:
            if wfs_class.name == selected_wfs_name:
                PageManager.set_wfs(self, wfs_class)

    telescope_combobox.bind("<<ComboboxSelected>>", on_telescope_selected)

    self.save_data_var = tk.IntVar()
    save_data_checkbutton = tk.Checkbutton(universal_frame, text="Save Star Data", variable= self.save_data_var)
    save_data_checkbutton.grid(row=4, column=0, padx=10, pady=5, sticky='w')

            
    # Adds new widgets
    row_index = 5
    for option in selected_options:
        if option == 'Plot Inset Graph':
            frame = tk.Frame(window, bd=1, relief="solid")
            frame.grid(sticky='nsew', padx=5, pady=5)

            tk.Label(frame, text=option).grid(row=row_index, column=0, padx=10, pady=10, sticky='w')
            row_index += 1

            tk.Label(frame, text="l Point of Interest").grid(row=row_index, column=0, padx=10, pady=5, sticky='w')
            self.l_poi_entry = tk.Entry(frame)
            self.l_poi_entry.grid(row=row_index, column=1, padx=5, pady=5, sticky='ew')
            row_index += 1

            tk.Label(frame, text="b Point of Interest").grid(row=row_index, column=0, padx=10, pady=5, sticky='w')
            self.b_poi_entry = tk.Entry(frame)
            self.b_poi_entry.grid(row=row_index, column=1, padx=5, pady=5, sticky='ew')
            row_index += 1

        elif option == 'Plot SNR for Science Camera' or option == 'Plot SNR for Science Camera + Noise':
            frame = tk.Frame(window, bd=1, relief="solid")
            frame.grid(sticky='nsew', padx=5, pady=5)

            tk.Label(frame, text=option).grid(row=row_index, column=0, padx=10, pady=10, sticky='w')
            row_index += 1

            selected_sci_camera = tk.StringVar()
            tk.Label(frame, text="Science Camera").grid(row=row_index, column=0, padx=10, pady=5, sticky='w')
            sci_cam_combobox = ttk.Combobox(frame, textvariable = selected_sci_camera, state="readonly")
            sci_cam_combobox['values'] = [sci_camera.name for sci_camera in sci_camera_classes]
            sci_cam_combobox.grid(row=row_index, column=1, padx=5, pady=5, sticky='ew')
            row_index += 1

            tk.Label(frame, text="Exposure Time").grid(row=row_index, column=0, padx=10, pady=5, sticky='w')
            self.exp_time_entry = tk.Entry(frame)
            self.exp_time_entry.grid(row=row_index, column=1, padx=5, pady=5, sticky='ew')
            row_index += 1

            self.log_yaxis_var = tk.IntVar()
            log_yaxis_checkbutton = tk.Checkbutton(frame, text="Log y axis", variable=self.log_yaxis_var)
            log_yaxis_checkbutton.grid(row=row_index, column=0, padx=10, pady=5, sticky='w')
            row_index += 1

            self.save_plots_var = tk.IntVar()
            save_plots_checkbutton = tk.Checkbutton(frame, text="Save plots", variable=self.save_plots_var)
            save_plots_checkbutton.grid(row=row_index, column=0, padx=10, pady=5, sticky='w')
            row_index += 1

            sci_cam_combobox.bind("<<ComboboxSelected>>", on_sci_camera_selected)

        elif option == 'Plot SHWF Sesnor Camera SNR' or option == 'Plot SHWF Sesnor Camera SNR + Noise':
            frame = tk.Frame(window, bd=1, relief="solid")
            frame.grid(sticky='nsew', padx=5, pady=5)

            tk.Label(frame, text=option).grid(row=row_index, column=0, padx=10, pady=10, sticky='w')
            row_index += 1

            selected_wfs = tk.StringVar()
            tk.Label(frame, text="Wavefront Sensor").grid(row=row_index, column=0, padx=10, pady=5, sticky='w')
            shwf_sen_combobox = ttk.Combobox(frame,textvariable = selected_wfs, state="readonly")
            shwf_sen_combobox['values'] = [wfs.name for wfs in wfs_classes]
            shwf_sen_combobox.grid(row=row_index, column=1, padx=5, pady=5, sticky='ew')
            row_index += 1

            selected_wfs_camera = tk.StringVar()
            tk.Label(frame, text="Wavefront Sensor Camera").grid(row=row_index, column=0, padx=10, pady=5, sticky='w')
            shwf_cam_combobox = ttk.Combobox(frame,textvariable = selected_wfs_camera, state="readonly")
            shwf_cam_combobox['values'] = [wfs_camera.name for wfs_camera in wfs_camera_classes]
            shwf_cam_combobox.grid(row=row_index, column=1, padx=5, pady=5, sticky='ew')
            row_index += 1

            tk.Label(frame, text="Exposure Time").grid(row=row_index, column=0, padx=10, pady=5, sticky='w')
            self.shwf_exp_time_entry = tk.Entry(frame)
            self.shwf_exp_time_entry.grid(row=row_index, column=1, padx=5, pady=5, sticky='ew')
            row_index += 1

            self.log_yaxis_shwf_var = tk.IntVar()
            log_yaxis_shwf_checkbutton = tk.Checkbutton(frame, text="Log y axis", variable=self.log_yaxis_shwf_var)
            log_yaxis_shwf_checkbutton.grid(row=row_index, column=0, padx=10, pady=5, sticky='w')
            row_index += 1

            self.save_plots_shwf_var = tk.IntVar()
            save_plots_shwf_checkbutton = tk.Checkbutton(frame, text="Save plots", variable=self.save_plots_shwf_var)
            save_plots_shwf_checkbutton.grid(row=row_index, column=0, padx=10, pady=5, sticky='w')
            row_index += 1

            shwf_cam_combobox.bind("<<ComboboxSelected>>", on_wfs_camera_selected)
            shwf_sen_combobox.bind("<<ComboboxSelected>>", on_wfs_selected)

        elif option == 'Quantitative Random Star Search':
            frame = tk.Frame(window, bd=1, relief="solid")
            frame.grid(sticky='nsew', padx=5, pady=5)

            tk.Label(frame, text=option).grid(row=row_index, column=0, padx=10, pady=10, sticky='w')
            row_index += 1

            tk.Label(frame, text="Number of Searches").grid(row=row_index, column=0, padx=10, pady=5, sticky='w')
            self.num_searches_entry = tk.Entry(frame)
            self.num_searches_entry.grid(row=row_index, column=1, padx=5, pady=5, sticky='ew')
            row_index += 1

    button_next_2 = tk.Button(window, text="Next", command=lambda: execute(PageManager))
    button_next_2.grid(row=row_index, column=0, padx=5, pady=5)

def execute(self):

    if self.save_data_var == 0:
        save_data = False
    else:
        save_data = True

    star_data = fetch_data(self.min_mag_entry.get(), self.max_mag_entry.get(), save_data)

    for option in selected_options:
        if option == 'Plot Inset Graph':
            plot_inset_graph(star_data, self.telescope, float(self.l_poi_entry.get()), float(self.b_poi_entry.get()))

        elif option == 'Plot Mignitude Histogram':
            plot_mag_hist(star_data)

        elif option == 'Plot SNR for Science Camera':

            plot_snr_standard(star_data, float(self.exp_time_entry.get()), self.telescope,self.sci_camera, log_yscale = True)

        elif option == 'Plot SNR for Science Camera + Noise':

            plot_snr_standard(star_data, float(self.exp_time_entry.get()), self.telescope,self.sci_camera, log_yscale = True)
            plot_standard_noise_contributions(star_data, float(self.exp_time_entry.get()), self.telescope,self.sci_camera)


        elif option == 'Plot SHWF Sesnor Camera SNR':
            plot_shwf_snr_standard(star_data, float(self.shwf_exp_time_entry.get()), self.telescope,self.wfs_camera,self.wfs, log_yscale = True)



        elif option == 'Plot SHWF Sesnor Camera SNR + Noise':
            plot_shwf_snr_standard(star_data, float(self.shwf_exp_time_entry.get()), self.telescope,self.wfs_camera,self.wfs, log_yscale = True)
            plot_shwf_noise_contributions(star_data, float(self.shwf_exp_time_entry.get()), self.telescope,self.wfs_camera,self.wfs)

        elif option == 'Quantitative Random Star Search':
            mag_threshold= calc_mag(float(self.shwf_exp_time_entry.get()), self.telescope ,min_photon_standard(float(self.shwf_exp_time_entry.get()),self.telescope,self.wfs_camera))
            search_stars(star_data, self.telescope, mag_threshold, num_searches = int(self.num_searches_entry.get()))
            pass

    pass



def toggle_sub_checkbox(parent_var, sub_checkbutton):
    if parent_var.get() == 1:
        sub_checkbutton.configure(state="normal")
    else:
        sub_checkbutton.configure(state="disabled")

def check_options():
    global selected_options

    selected_options = []
    
    if var1.get() == 1:
        selected_options.append('Plot Inset Graph')

    if var2.get() == 1:
        selected_options.append('Plot Magnitude Histogram')
        
    if var3.get() == 1 and var4.get() == 1:
        selected_options.append('Plot SNR for Science Camera + Noise')
    if var3.get() == 1 and var4.get() != 1:
        selected_options.append('Plot SNR for Science Camera')

    if var5.get() == 1 and var6.get() == 1:
        selected_options.append('Plot SHWF Sesnor Camera SNR + Noise')
    if var5.get() == 1 and var6.get() != 1:
        selected_options.append('Plot SHWF Sesnor Camera SNR')

    if var7.get() == 1:
        selected_options.append('Quantitative Random Star Search')

    if len(selected_options) == 0:
        messagebox.showerror("Error", "No checkboxes are selected. Please select at least one.")
    else:
        open_next_page(PageManager)


window = tk.Tk()
window.title("My Program")

var1 = tk.IntVar()
check1 = tk.Checkbutton(window, text='Plot Inset Graph', variable=var1, anchor="w", justify="left")
check1.pack(fill='x')

var2 = tk.IntVar()
check2 = tk.Checkbutton(window, text='Plot Magnitude Histogram', variable=var2, anchor="w", justify="left")
check2.pack(fill='x')

var3 = tk.IntVar()
check3 = tk.Checkbutton(window, text='Plot SNR for Science Camera', variable=var3, anchor="w", justify="left", command=lambda: toggle_sub_checkbox(var3, check4))
check3.pack(fill='x')

var4 = tk.IntVar()
check4 = tk.Checkbutton(window, text='Plot Noise Contributions', variable=var4, anchor="w", justify="left", state="disabled")
check4.pack(fill='x', padx=20)

var5 = tk.IntVar()
check5 = tk.Checkbutton(window, text='Plot SHWF Sesnor Camera SNR', variable=var5, anchor="w", justify="left", command=lambda: toggle_sub_checkbox(var5, check6))
check5.pack(fill='x')

var6 = tk.IntVar()
check6 = tk.Checkbutton(window, text='Plot Noise Contributions', variable=var6, anchor="w", justify="left", state="disabled")
check6.pack(fill='x', padx=20)

var7 = tk.IntVar()
check7 = tk.Checkbutton(window, text='Quantitative Random Star Search', variable=var7, anchor="w", justify="left")
check7.pack(fill='x')

button_next = tk.Button(window, text="Next", command=check_options)
button_next.pack()

window.mainloop()


