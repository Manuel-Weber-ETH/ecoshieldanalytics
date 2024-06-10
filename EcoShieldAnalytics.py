import tkinter as tk
from tkinter import filedialog, ttk
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon, LineString, MultiLineString, GeometryCollection, Point
import gpxpy
from fastkml import kml
import csv
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from PIL import Image, ImageTk
import sys
import os
from datetime import datetime
from mpl_toolkits.axes_grid1 import make_axes_locatable

class GridApp:
    def __init__(self, root):
        self.root = root
        self.setup_ui()
        self.study_area = None
        self.gps_tracks = None
        self.grid = None
        self.dynamic_rows = []
        self.hq_entries = []

        # Initialize variables for negative and positive event indicators
        self.mean_coverage_var20 = tk.StringVar(value="0.00")
        self.baseline_coverage_var20 = tk.DoubleVar(value=0.0)
        self.mean_score_var20 = tk.StringVar(value="0.00%")

        self.mean_coverage_var30 = tk.StringVar(value="0.00")
        self.baseline_coverage_var30 = tk.DoubleVar(value=0.0)
        self.mean_score_var30 = tk.StringVar(value="0.00%")

    
    def setup_ui(self):
        self.root.title("EcoShield Analytics")

        # Main layout division: Control Frame (Left) & Visualization Frame (Right)
        control_frame = ttk.Frame(self.root)
        control_frame.grid(row=0, column=0, sticky="nsew", padx=5, pady=5)

        visualization_frame = ttk.Frame(self.root)
        visualization_frame.grid(row=0, column=1, sticky="nsew", padx=5, pady=5)

        # Configuring the root window's columns and rows for proper resizing
        self.root.columnconfigure(0, weight=1)
        self.root.columnconfigure(1, weight=3)  # Visualization frame takes more space
        self.root.rowconfigure(0, weight=1)

        # Top Section for Shapefile and Grid Resolution within Control Frame
        top_frame = ttk.Frame(control_frame)
        top_frame.grid(row=0, column=0, sticky="ew", padx=10, pady=5, columnspan=5)

        ttk.Label(top_frame, text="Protected area shapefile (.shp):").grid(row=0, column=0, sticky="w")
        ttk.Button(top_frame, text="Browse", command=self.load_shp).grid(row=0, column=1, padx=5)
        
        ttk.Label(top_frame, text="Grid cell resolution (meters):").grid(row=0, column=2, sticky="w")
        self.grid_resolution_var = tk.StringVar()
        ttk.Entry(top_frame, textvariable=self.grid_resolution_var, width=10).grid(row=0, column=3, padx=5)

        # Input for HQ Coordinates
        self.hq_frame = ttk.Frame(top_frame)
        self.hq_frame.grid(row=2, column=0, columnspan=7, sticky="ew", padx=5, pady=5)

        ttk.Label(self.hq_frame, text="Post Name").grid(row=0, column=0, sticky="w")
        ttk.Label(self.hq_frame, text="Longitude").grid(row=0, column=1, sticky="w")
        ttk.Label(self.hq_frame, text="Latitude").grid(row=0, column=2, sticky="w")
        ttk.Label(self.hq_frame, text="Coefficient").grid(row=0, column=3, sticky="w")

        self.hq_name_var = tk.StringVar()
        self.hq_longitude_var = tk.DoubleVar()
        self.hq_latitude_var = tk.DoubleVar()
        self.hq_coefficient_var = tk.DoubleVar(value=1)

        # Button to add more HQs
        ttk.Button(self.hq_frame, text="Add Post", command=self.add_hq).grid(row=0, column=4, padx=5)

        # Button to import HQs from csv
        ttk.Button(self.hq_frame, text="Import posts from .csv", command=self.import_posts_from_csv).grid(row=0, column=5, padx=5)

        # Button to import HQs from csv
        ttk.Button(self.hq_frame, text="Export posts to .csv", command=self.export_posts_to_csv).grid(row=0, column=6, padx=5)

        # Indicator Table Setup within Control Frame
        self.indicator_frame = ttk.Frame(control_frame)
        self.indicator_frame.grid(row=2, column=0, sticky="nsew", padx=10, pady=5, columnspan=5)

        titles = ["Indicator Type", "Indicator Name", "Target value", "Value", "Score", "Coefficient"]
        column_widths = [15, 10, 8, 8, 8, 8]  # Example widths in characters

        for i, (title, width) in enumerate(zip(titles, column_widths)):
            label = ttk.Label(self.indicator_frame, text=title)
            label.grid(row=0, column=i, padx=5, pady=5)
            self.indicator_frame.columnconfigure(i, minsize=width*10)  # Width in pixels (rough estimate)

        # Button to add a new row
        ttk.Button(self.indicator_frame, text="Add Indicator", command=self.add_row).grid(row=0, column=6, padx=5, pady=5)

        # Button to import indicators from csv
        ttk.Button(self.indicator_frame, text="Import indicators from .csv", command=self.import_indicators_from_csv).grid(row=0, column=7, padx=5, pady=5)

        control_frame.columnconfigure(0, weight=1)

        # Add progress_text to display progress updates
        progress_frame = ttk.Frame(control_frame)
        progress_frame.grid(row=3, column=0, sticky="nsew", padx=10, pady=5, columnspan=5)
        control_frame.rowconfigure(3, weight=1)  # Make the progress frame expandable

        self.progress_text = tk.Text(progress_frame, height=5, width=40)
        self.progress_text.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # Visualization Setup on the Right Side
        self.fig, self.ax = plt.subplots(figsize=(6, 4))
        self.canvas = FigureCanvasTkAgg(self.fig, master=visualization_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # Button for exporting the table to CSV
        ttk.Button(top_frame, text="Export indicators to .csv", command=self.export_indicators_to_csv).grid(row=0, column=4, sticky="w")

        # Button to export grid
        ttk.Button(top_frame, text="Save grid as .shp", command=self.save_grid).grid(row=0, column=5, sticky="w")

        # Button to save current plot
        ttk.Button(top_frame, text="Save plot as .png", command=self.save_plot).grid(row=0, column=6, sticky="w")

        # Label and entry for the overall score
        self.weighted_mean_var = tk.StringVar(value="0.00%")
        self.weighted_mean_entry = ttk.Entry(visualization_frame, textvariable=self.weighted_mean_var, state='readonly', width=10)
        self.weighted_mean_entry.pack(side=tk.BOTTOM, padx=5, pady=2)
        ttk.Label(visualization_frame, text="Overall Control Score:").pack(side=tk.BOTTOM, padx=5, pady=2)


    def save_plot(self):
        # Open a file dialog to select where to save the file
        file_path = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[("PNG files", "*.png"), ("All files", "*.*")],
            title="Save plot as"
        )
        
        # Check if a file path was provided (i.e., the user didn't cancel the dialog)
        if file_path:
            # Save the current figure to the specified file path
            self.fig.savefig(file_path, format='png')
            self.update_progress(f"Plot saved as {file_path}")


    def export_indicators_to_csv(self):
        today_date = datetime.today().strftime('%Y-%m-%d')
        indicators_file_name = filedialog.asksaveasfilename(
            initialfile=f"ecoshieldanalytics_{today_date}_indicators.csv",
            defaultextension=".csv",
            filetypes=[("CSV (Comma delimited)", "*.csv")]
        )
        
        if not indicators_file_name:  # User cancelled the dialog
            return

        # Export indicators to CSV
        with open(indicators_file_name, mode='w', newline='') as indicators_file:
            csv_writer = csv.writer(indicators_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            csv_writer.writerow(["Indicator Name", "Target value", "Value", "Score [%]", "Coefficient", "Overall Control Score"])
            
            for vars, _ in self.dynamic_rows:
                csv_writer.writerow([
                    vars["name_var"].get(), 
                    vars["baseline_value_var"].get(), 
                    vars["value_var"].get(), 
                    vars["score_var"].get().rstrip('%'),  # Remove '%' for consistency
                    vars["coefficient_var"].get(), 
                    self.weighted_mean_var.get().rstrip('%')  # Include overall score in each row for simplicity
                ])

        # Optionally, update the user about the export
        self.update_progress("Posts and indicators exported to CSV successfully.")

    def export_posts_to_csv(self):
        today_date = datetime.today().strftime('%Y-%m-%d')
        posts_file_name = filedialog.asksaveasfilename(
            initialfile=f"ecoshieldanalytics_{today_date}_posts.csv",
            defaultextension=".csv",
            filetypes=[("CSV (Comma delimited)", "*.csv")]
        )
    
        if not posts_file_name:  # User cancelled the dialog
            return

        # Export posts to CSV
        with open(posts_file_name, mode='w', newline='') as posts_file:
            csv_writer = csv.writer(posts_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            csv_writer.writerow(["Post Name", "Longitude", "Latitude", "Coefficient"])
            
            for _, hq_name_entry, hq_longitude_entry, hq_latitude_entry, hq_coefficient_entry, _ in self.hq_entries:
                csv_writer.writerow([
                    hq_name_entry.get(), 
                    hq_longitude_entry.get(), 
                    hq_latitude_entry.get(), 
                    hq_coefficient_entry.get()
                ])

        # Optionally, update the user about the export
        self.update_progress("Posts exported to CSV successfully.")


    def import_posts_from_csv(self):
        file_path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if not file_path:  # User cancelled the dialog
            return

        try:
            df = pd.read_csv(file_path)
            required_columns = ["Post Name", "Longitude", "Latitude", "Coefficient"]
            for col in required_columns:
                if col not in df.columns:
                    df[col] = ""  # Add missing columns with empty values

            # Remove any additional columns
            df = df[required_columns]

            # Clear existing HQ entries
            for _, hq_name_entry, hq_longitude_entry, hq_latitude_entry, hq_coefficient_entry, remove_button in self.hq_entries:
                hq_name_entry.grid_forget()
                hq_longitude_entry.grid_forget()
                hq_latitude_entry.grid_forget()
                hq_coefficient_entry.grid_forget()
                remove_button.grid_forget()
            self.hq_entries.clear()

            # Add new rows from CSV
            for _, row in df.iterrows():
                self.add_hq()
                self.hq_entries[-1][1].insert(0, row["Post Name"])
                self.hq_entries[-1][2].insert(0, row["Longitude"])
                self.hq_entries[-1][3].insert(0, row["Latitude"])
                self.hq_entries[-1][4].insert(0, row["Coefficient"])
            self.update_progress("Posts imported successfully from CSV.")
        except Exception as e:
            self.update_progress(f"Error importing posts from CSV: {e}")

    def import_indicators_from_csv(self):
        file_path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if not file_path:  # User cancelled the dialog
            return

        try:
            df = pd.read_csv(file_path)
            required_columns = ["Indicator Name", "Target value", "Value", "Score [%]", "Coefficient"]
            for col in required_columns:
                if col not in df.columns:
                    df[col] = ""  # Add missing columns with empty values

            # Remove any additional columns
            df = df[required_columns]

            # Clear existing indicator rows
            for _, widgets in self.dynamic_rows:
                for widget in widgets.values():
                    widget.destroy()
            self.dynamic_rows.clear()

            # Add new rows from CSV
            for _, row in df.iterrows():
                self.add_row()
                last_row_vars, last_row_widgets = self.dynamic_rows[-1]
                last_row_vars["indicator_type_var"].set("4. Spatially implicit indicator")  # Set the indicator type
                last_row_vars["name_var"].set(row["Indicator Name"])
                last_row_vars["baseline_value_var"].set(row["Target value"])
                last_row_vars["value_var"].set(row["Value"])
                last_row_vars["score_var"].set(row["Score [%]"])
                last_row_vars["coefficient_var"].set(row["Coefficient"])

                # Ensure that the correct buttons are displayed based on the indicator type
                self.update_row_functionality(last_row_vars, len(self.dynamic_rows))

            self.update_progress("Indicators imported successfully from CSV.")
        except Exception as e:
            self.update_progress(f"Error importing indicators from CSV: {e}")










    def add_hq(self):
        row = len(self.hq_entries) + 2  # Start after the initial row

        hq_longitude_var = tk.DoubleVar(value=0.0)
        hq_latitude_var = tk.DoubleVar(value=0.0)
        hq_coefficient_var = tk.DoubleVar(value=1)
        hq_name_var = tk.StringVar()

        hq_name_entry = ttk.Entry(self.hq_frame, textvariable=hq_name_var, width=10)
        hq_name_entry.grid(row=row, column=0, padx=5)

        hq_longitude_entry = ttk.Entry(self.hq_frame, textvariable=hq_longitude_var, width=10)
        hq_longitude_entry.grid(row=row, column=1, padx=5)
        
        hq_latitude_entry = ttk.Entry(self.hq_frame, textvariable=hq_latitude_var, width=10)
        hq_latitude_entry.grid(row=row, column=2, padx=5)
        
        hq_coefficient_entry = ttk.Entry(self.hq_frame, textvariable=hq_coefficient_var, width=10)
        hq_coefficient_entry.grid(row=row, column=3, padx=5)
        
        remove_button = ttk.Button(self.hq_frame, text="Remove", command=lambda r=row: self.remove_hq(r))
        remove_button.grid(row=row, column=4, padx=5)

        self.hq_entries.append((row, hq_name_entry, hq_longitude_entry, hq_latitude_entry, hq_coefficient_entry, remove_button))

    def remove_hq(self, row):
        for i, entry in enumerate(self.hq_entries):
            if entry[0] == row:
                for widget in entry[1:]:
                    widget.grid_forget()
                self.hq_entries.pop(i)
                break

        # Update the grid positions of the remaining entries
        for new_row, entry in enumerate(self.hq_entries, start=2):
            entry[0] = new_row
            entry[1].grid(row=new_row, column=0, padx=5)
            entry[2].grid(row=new_row, column=1, padx=5)
            entry[3].grid(row=new_row, column=2, padx=5)
            entry[4].grid(row=new_row, column=3, padx=5)
            entry[5].grid(row=new_row, column=4, padx=5)







    def load_negative_incidents_csv(self):
        filepath = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if not filepath:
            return  # User cancelled the dialog
        try:
            incidents_df = pd.read_csv(filepath)
            incidents_gdf = gpd.GeoDataFrame(incidents_df, geometry=gpd.points_from_xy(incidents_df['Longitude'], incidents_df['Latitude']), crs="EPSG:4326")
            
            # Spatial join incidents with the grid to count incidents per cell
            self.grid = self.grid.to_crs("EPSG:4326")  # Ensure CRS matches for spatial join
            joined = gpd.sjoin(self.grid, incidents_gdf, how="left", predicate="contains")
            
            # Count incidents per grid cell
            incidents_count = joined.groupby(joined.index).size()
            self.grid['negative_event'] = incidents_count.reindex(self.grid.index, fill_value=0)

            # Perform calculations as specified
            mean_distance = self.grid['distance'].mean()
            self.grid['calculation'] = (self.grid['negative_event'] / self.grid['distance']) * mean_distance
            mean_calculation = self.grid.loc[self.grid['coverage'] > 0, 'calculation'].mean()

            # Round mean_calculation to two decimal places
            mean_calculation = round(mean_calculation, 2)

            # Update UI
            self.mean_coverage_var20.set(mean_calculation)

            # Fetch the baseline value from the appropriate dynamic row
            for vars, _ in self.dynamic_rows:
                if vars["indicator_type_var"].get() == "2. Negative events [/cell/d]":
                    baseline_value = vars["baseline_value_var"].get()
                    break
            else:
                baseline_value = 0  # Default value if no matching row is found

            self.update_progress(f"Baseline Value (Negative): {baseline_value}")
            
            percentage = 0  # Initialize percentage variable
            if baseline_value > 0:  # Ensure we don't divide by zero
                percentage = (1 - ((mean_calculation / baseline_value))) * 100
                if percentage > 100:
                    percentage = 100
                self.mean_score_var20.set(f"{percentage:.2f}%")
            else:
                self.mean_score_var20.set("Baseline = 0")

            # Update the corresponding dynamic row
            for vars, _ in self.dynamic_rows:
                if vars["indicator_type_var"].get() == "2. Negative events [/cell/d]":
                    vars["value_var"].set(f"{mean_calculation:.2f}")
                    vars["score_var"].set(f"{percentage:.2f}%")
                    break
                    
            if hasattr(self, 'ax'):
                self.negative_incidents_gdf = incidents_gdf
                self.plot_incidents()
        except Exception as e:
            self.update_progress(f"Failed to load or process CSV: {e}")



    def load_positive_incidents_csv(self):
        filepath = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if not filepath:
            return  # User cancelled the dialog
        try:
            incidents_df = pd.read_csv(filepath)
            incidents_gdf = gpd.GeoDataFrame(incidents_df, geometry=gpd.points_from_xy(incidents_df['Longitude'], incidents_df['Latitude']), crs="EPSG:4326")
            
            self.grid = self.grid.to_crs("EPSG:4326")  # Ensure CRS matches for spatial join
            joined = gpd.sjoin(self.grid, incidents_gdf, how="left", predicate="contains")
            incidents_count = joined.groupby(joined.index).size()
            self.grid['positive_event'] = incidents_count.reindex(self.grid.index, fill_value=0)

            # Perform calculations as specified for positive incidents
            mean_distance = self.grid['distance'].mean()
            self.grid['calculation'] = (self.grid['positive_event'] / mean_distance) * self.grid['distance']
            mean_calculation = self.grid.loc[self.grid['coverage'] > 0, 'calculation'].mean()

            # Round mean_calculation to two decimal places
            mean_calculation = round(mean_calculation, 2)

            # Update UI
            self.mean_coverage_var30.set(mean_calculation)

            # Fetch the baseline value from the appropriate dynamic row
            for vars, _ in self.dynamic_rows:
                if vars["indicator_type_var"].get() == "3. Positive events [/cell*d]":
                    baseline_value = vars["baseline_value_var"].get()
                    break
            else:
                baseline_value = 0  # Default value if no matching row is found

            self.update_progress(f"Baseline Value (Positive): {baseline_value}")
            
            percentage = 0  # Initialize percentage variable
            # Calculate percentage of baseline coverage
            if baseline_value > 0:  # Ensure we don't divide by zero
                percentage = (mean_calculation / baseline_value) * 100
                if percentage > 100:
                    percentage = 100
                self.mean_score_var30.set(f"{percentage:.2f}%")
            else:
                self.mean_score_var30.set("Baseline = 0")

            # Update the corresponding dynamic row
            for vars, _ in self.dynamic_rows:
                if vars["indicator_type_var"].get() == "3. Positive events [/cell*d]":
                    vars["value_var"].set(f"{mean_calculation:.2f}")
                    vars["score_var"].set(f"{percentage:.2f}%")
                    break

            if hasattr(self, 'ax'):
                self.positive_incidents_gdf = incidents_gdf
                self.plot_incidents()
        except Exception as e:
            self.update_progress(f"Failed to load or process CSV: {e}")


    def remove_specific_row(self, row_index):
        # Adjust row_index to match the internal index
        internal_index = row_index - 1

        # Find the row to remove by its index in the dynamic_rows list
        if 0 <= internal_index < len(self.dynamic_rows):
            _, widgets = self.dynamic_rows.pop(internal_index)
            for widget in widgets.values():
                widget.grid_forget()  # Use grid_forget instead of destroy

            # Re-grid remaining rows to fill the gap
            for new_index, (vars, new_widgets) in enumerate(self.dynamic_rows):
                for i, (key, widget) in enumerate(new_widgets.items()):
                    widget.grid(row=new_index + 1, column=i, padx=5, pady=5)

                # Ensure that the correct buttons are displayed based on the indicator type
                self.update_row_functionality(vars, new_index + 1)

            # Update the commands for the remove buttons to keep row indexes correct
            for index, (_, new_widgets) in enumerate(self.dynamic_rows, start=1):
                new_widgets["remove_button"].config(command=lambda idx=index: self.remove_specific_row(idx))

            # Update the baseline values for negative and positive events
            self.update_baseline_values()


    def update_baseline_values(self):
        for vars, _ in self.dynamic_rows:
            if vars["indicator_type_var"].get() == "2. Negative events [/cell/d]":
                self.baseline_coverage_var20.set(vars["baseline_value_var"].get())
            elif vars["indicator_type_var"].get() == "3. Positive events [/cell*d]":
                self.baseline_coverage_var30.set(vars["baseline_value_var"].get())




    def update_weighted_mean_score(self):
        total_weighted_score = 0
        total_weight = 0

        # Process dynamically added rows
        for vars, widgets in self.dynamic_rows:
            try:
                coefficient_dynamic = vars["coefficient_var"].get()
                score_str_dynamic = vars["score_var"].get().rstrip('%')
                if score_str_dynamic and score_str_dynamic != "N/A" and score_str_dynamic != "" and coefficient_dynamic:
                    score_dynamic = float(score_str_dynamic)
                    total_weighted_score += score_dynamic * coefficient_dynamic
                    total_weight += coefficient_dynamic
            except ValueError as e:
                print(f"Error processing dynamic row scores: {e}")

        # Calculate the overall score
        if total_weight > 0:
            overall_score = total_weighted_score / total_weight
            self.weighted_mean_var.set(f"{overall_score:.2f}%")
        else:
            self.weighted_mean_var.set("N/A")



    

    def visualize_grid(self):
        self.ax.clear()

        # Remove any existing colorbars from the figure
        fig = self.ax.get_figure()
        for ax in fig.get_axes()[1:]:  # Remove all axes except the main one
            ax.remove()

        # Plot the grid with a different color gradient and thinner lines
        if hasattr(self, 'grid'):
            divider = make_axes_locatable(self.ax)
            cax = divider.append_axes("right", size="5%", pad=0.1)

            self.grid.plot(column='track_length', ax=self.ax, legend=False, cmap='binary', edgecolor='black', linewidth=0.1)
            cbar = fig.colorbar(plt.cm.ScalarMappable(cmap='binary', norm=plt.Normalize(vmin=self.grid['track_length'].min(), vmax=self.grid['track_length'].max())), cax=cax)
            cbar.set_label("Patrol length within cell [m]")
            self.ax.set_title("Patrol coverage")

        # Add HQs as large black points with their names as labels
        hq_points = []
        if self.hq_longitude_var.get() != 0.0 or self.hq_latitude_var.get() != 0.0:
            hq_points.append((float(self.hq_longitude_var.get()), float(self.hq_latitude_var.get()), self.hq_name_var.get()))

        for _, hq_name_entry, hq_longitude_entry, hq_latitude_entry, hq_coefficient_entry, remove_button in self.hq_entries:
            hq_points.append((float(hq_longitude_entry.get()), float(hq_latitude_entry.get()), hq_name_entry.get()))

        for lon, lat, name in hq_points:
            self.ax.scatter(lon, lat, color='black', s=10, label='HQ', zorder=5)
            self.ax.text(lon, lat, name, color='black', fontsize=12, ha='right')

        # Plot the actual tracks as very fine black lines
        if hasattr(self, 'gps_tracks'):
            self.gps_tracks.plot(ax=self.ax, color='black', linewidth=0.5, alpha=0.7)

        self.canvas.draw()

    def plot_incidents(self):
        # Clear the previous plot
        self.ax.clear()

        # Remove any existing colorbars from the figure
        fig = self.ax.get_figure()
        for ax in fig.get_axes()[1:]:  # Remove all axes except the main one
            ax.remove()

        # Plot the grid with patrol lengths
        if hasattr(self, 'grid'):
            divider = make_axes_locatable(self.ax)
            cax = divider.append_axes("right", size="5%", pad=0.1)

            self.grid.plot(column='track_length', ax=self.ax, legend=False,
                        cmap='binary', edgecolor='black', linewidth=0.1)
            cbar = fig.colorbar(plt.cm.ScalarMappable(cmap='binary', norm=plt.Normalize(vmin=self.grid['track_length'].min(), vmax=self.grid['track_length'].max())), cax=cax)
            cbar.set_label("Patrol length within cell [m]")
            self.ax.set_title("Patrol coverage")

        # Plot HQs with names as labels
        hq_points = []
        if self.hq_longitude_var.get() != 0.0 or self.hq_latitude_var.get() != 0.0:
            hq_points.append((float(self.hq_longitude_var.get()), float(self.hq_latitude_var.get()), self.hq_name_var.get()))

        for _, hq_name_entry, hq_longitude_entry, hq_latitude_entry, hq_coefficient_entry, remove_button in self.hq_entries:
            hq_points.append((float(hq_longitude_entry.get()), float(hq_latitude_entry.get()), hq_name_entry.get()))

        for lon, lat, name in hq_points:
            self.ax.scatter(lon, lat, color='black', s=10, label='HQ', zorder=5)
            self.ax.text(lon, lat, name, color='black', fontsize=12, ha='right')

        # Plot the actual tracks as very fine black lines
        if hasattr(self, 'gps_tracks'):
            self.gps_tracks.plot(ax=self.ax, color='black', linewidth=0.5, alpha=0.7)

        # Plot negative incidents
        if hasattr(self, 'negative_incidents_gdf'):
            negative_lons = self.negative_incidents_gdf.geometry.x
            negative_lats = self.negative_incidents_gdf.geometry.y
            self.ax.scatter(negative_lons, negative_lats, color='red', label='Negative Incidents', alpha=0.6, edgecolor='none')

        # Plot positive incidents
        if hasattr(self, 'positive_incidents_gdf'):
            positive_lons = self.positive_incidents_gdf.geometry.x
            positive_lats = self.positive_incidents_gdf.geometry.y
            self.ax.scatter(positive_lons, positive_lats, color='green', label='Positive Incidents', alpha=0.6, edgecolor='none')

        self.canvas.draw()



    def export_table_to_csv(self):
        # Define the CSV file name
        csv_file_name = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV (Comma delimited)", "*.csv")])
        if not csv_file_name:  # User cancelled the dialog
            return

        with open(csv_file_name, mode='w', newline='') as csv_file:
            csv_writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            
            # Write the header
            csv_writer.writerow(["Indicator Name", "Baseline Value", "Value", "Score [%]", "Coefficient", "Overall Control Score"])
            
            # Write rows for the dynamic indicators
            for vars, _ in self.dynamic_rows:
                csv_writer.writerow([vars["name_var"].get(), 
                    vars["baseline_value_var"].get(), 
                    vars["value_var"].get(), 
                    vars["score_var"].get().rstrip('%'),  # Remove '%' for consistency
                    vars["coefficient_var"].get(), 
                    self.weighted_mean_var.get().rstrip('%')])  # Include overall score in each row for simplicity


        # Optionally, update the user about the export
        self.update_progress("Table exported to CSV successfully.")


    # Add a method to update the progress text
    def update_progress(self, message):
        self.progress_text.insert(tk.END, message + "\n")
        self.progress_text.see(tk.END)

    def calculate_dynamic_score(self, value_var, baseline_var, score_var):
        try:
            value = value_var.get()
            baseline = baseline_var.get()
            if baseline > 0:  # Prevent division by zero
                score = (value / baseline) * 100
                if score > 100:
                    score = 100
                score_var.set(f"{score:.2f}%")
            else:
                score_var.set("N/A")
        except Exception as e:
            score_var.set("Error")
            print(f"Error calculating dynamic score: {e}")
        self.update_weighted_mean_score()


    
    def update_row_functionality(self, vars, row_index):
        indicator_type = vars["indicator_type_var"].get()
        widgets = self.dynamic_rows[row_index - 1][1]

        # Remove existing file button if it exists
        if "file_button" in widgets and widgets["file_button"]:
            widgets["file_button"].destroy()
            widgets["file_button"] = None

        if indicator_type == "2. Negative events [/cell/d]":
            widgets["file_button"] = ttk.Button(self.indicator_frame, text="Add Negative Events (.csv)", command=lambda: self.load_negative_incidents_csv_and_update_row(vars))
            widgets["file_button"].grid(row=row_index, column=7, padx=5, pady=5, sticky="ew")
            widgets["value_entry"].configure(state='readonly')  # Ensure value entry is read-only
        elif indicator_type == "3. Positive events [/cell*d]":
            widgets["file_button"] = ttk.Button(self.indicator_frame, text="Add Positive Events (.csv)", command=lambda: self.load_positive_incidents_csv_and_update_row(vars))
            widgets["file_button"].grid(row=row_index, column=7, padx=5, pady=5, sticky="ew")
            widgets["value_entry"].configure(state='readonly')  # Ensure value entry is read-only
        elif indicator_type == "1. Patrol [m/cell]":
            widgets["file_button"] = ttk.Button(self.indicator_frame, text="Add Patrol Tracks (.kml)", command=self.load_gps)
            widgets["file_button"].grid(row=row_index, column=7, padx=5, pady=5, sticky="ew")
            widgets["value_entry"].configure(state='readonly')  # Ensure value entry is read-only
        elif indicator_type == "4. Spatially implicit indicator":
            if "file_button" in widgets and widgets["file_button"]:
                widgets["file_button"].destroy()
                widgets["file_button"] = None
            widgets["value_entry"].configure(state='normal')  # Allow user to input value

        # Update the dynamic rows structure with potential changes
        self.dynamic_rows[row_index - 1] = (vars, widgets)

    def load_negative_incidents_csv_and_update_row(self, vars):
        self.load_negative_incidents_csv()
        vars["value_var"].set(self.mean_coverage_var20.get())
        vars["score_var"].set(self.mean_score_var20.get())

    def load_positive_incidents_csv_and_update_row(self, vars):
        self.load_positive_incidents_csv()
        vars["value_var"].set(self.mean_coverage_var30.get())
        vars["score_var"].set(self.mean_score_var30.get())





    def add_row(self):
        row_index = len(self.dynamic_rows) + 1  # starting index after existing rows
        vars = {
            "indicator_type_var": tk.StringVar(),
            "name_var": tk.StringVar(),
            "baseline_value_var": tk.DoubleVar(),
            "value_var": tk.DoubleVar(),
            "score_var": tk.StringVar(),
            "coefficient_var": tk.IntVar()
        }

        indicator_types = ['1. Patrol [m/cell]', '2. Negative events [/cell/d]', '3. Positive events [/cell*d]', '4. Spatially implicit indicator']
        widgets = {
            "type_menu": ttk.OptionMenu(self.indicator_frame, vars["indicator_type_var"], indicator_types[0], *indicator_types, command=lambda event: self.update_row_functionality(vars, row_index)),
            "name_entry": ttk.Entry(self.indicator_frame, textvariable=vars["name_var"], width=10),
            "baseline_entry": ttk.Entry(self.indicator_frame, textvariable=vars["baseline_value_var"], width=10),
            "value_entry": ttk.Entry(self.indicator_frame, textvariable=vars["value_var"], width=10, state='readonly'),  # Set as read-only initially
            "score_entry": ttk.Entry(self.indicator_frame, textvariable=vars["score_var"], state='readonly', width=10),  # Set as read-only
            "coefficient_entry": ttk.Entry(self.indicator_frame, textvariable=vars["coefficient_var"], width=10),
            "remove_button": ttk.Button(self.indicator_frame, text="Remove", command=lambda idx=row_index: self.remove_specific_row(idx)),
            "file_button": None  # Initialize the file_button key
        }

        # Trace callbacks to automatically compute score and update the overall score
        vars["value_var"].trace_add("write", lambda *args: self.calculate_dynamic_score(vars["value_var"], vars["baseline_value_var"], vars["score_var"]))
        vars["baseline_value_var"].trace_add("write", lambda *args: self.calculate_dynamic_score(vars["value_var"], vars["baseline_value_var"], vars["score_var"]))
        vars["coefficient_var"].trace_add("write", lambda *args: self.update_weighted_mean_score())

        # Grid the widgets, skip file_button as it is None initially
        for i, (key, widget) in enumerate(widgets.items()):
            if key != "file_button":
                if key != "type_menu":
                    widget.grid(row=row_index, column=i, padx=5, pady=5, sticky="ew")
                else:
                    widget.grid(row=row_index, column=0, padx=5, pady=5, sticky="ew")

        widgets["remove_button"].grid(row=row_index, column=6, padx=5, pady=5, sticky="ew")
        self.dynamic_rows.append((vars, widgets))

        # Initially update functionality to set the correct button commands
        self.update_row_functionality(vars, row_index)
        self.update_weighted_mean_score()  # Update overall score when a new row is added

        self.update_weighted_mean_score()  # Update overall score when a new row is added

        # Set the baseline values for negative and positive events
        if vars["indicator_type_var"].get() == "2. Negative events [/cell/d]":
            self.baseline_coverage_var20.set(vars["baseline_value_var"].get())
        elif vars["indicator_type_var"].get() == "3. Positive events [/cell*d]":
            self.baseline_coverage_var30.set(vars["baseline_value_var"].get())



    def load_shp(self):
        file_path = filedialog.askopenfilename(filetypes=[("Shapefiles", "*.shp")])
        if file_path:
            self.study_area = gpd.read_file(file_path)
            if self.study_area.crs != 'epsg:4326':
                self.study_area = self.study_area.to_crs('epsg:4326')
            self.update_progress("Shapefile loaded and projected to EPSG:4326")
        pass

    def load_gps(self):
        file_path = filedialog.askopenfilename(filetypes=[("GPS files", "*.kml *.gpx")])
        if file_path:
            if file_path.endswith('.kml'):
                self.gps_tracks = self.parse_kml(file_path)
            elif file_path.endswith('.gpx'):
                self.gps_tracks = self.parse_gpx(file_path)
            
            self.update_progress("GPS tracks loaded")

            # Check if both study area and GPS tracks are loaded before processing
            if hasattr(self, 'study_area') and self.gps_tracks is not None:
                self.process_data()
            else:
                self.update_progress("Please load a study area shapefile before processing GPS tracks.")
        pass

    def parse_kml(self, file_path):
        def extract_geometries(feature, lines):
            if hasattr(feature, 'features'):  # It's a Folder or Document
                for inner_feature in feature.features():
                    extract_geometries(inner_feature, lines)
            elif hasattr(feature, 'geometry'):  # It's a Placemark
                geom = feature.geometry
                if hasattr(geom, 'geoms'):  # MultiGeometry
                    for sub_geom in geom.geoms:
                        lines.append(LineString(sub_geom.coords))
                else:
                    lines.append(LineString(geom.coords))

        with open(file_path, 'rt', encoding="utf-8") as myfile:
            doc = myfile.read()
        k = kml.KML()
        k.from_string(doc.encode('utf-8'))
        lines = []
        for feature in k.features():
            extract_geometries(feature, lines)
        return gpd.GeoDataFrame(geometry=lines, crs='epsg:4326')

    def parse_gpx(self, file_path): ## this function is not functional
        gpx_file = open(file_path, 'r')
        gpx = gpxpy.parse(gpx_file)
        lines = []
        for track in gpx.tracks:
            for segment in track.segments:
                points = segment.points
                lines.append(LineString([(point.longitude, point.latitude) for point in points]))
        return gpd.GeoDataFrame(geometry=lines, crs='epsg:4326')

    def generate_grid(self, cell_size):
        # Determine the UTM zone for the study area's centroid for more accurate distance calculations
        centroid = self.study_area.unary_union.centroid
        utm_zone = int((centroid.x + 180) / 6) + 1
        utm_crs = f'+proj=utm +zone={utm_zone} +ellps=WGS84 +datum=WGS84 +units=m +no_defs'

        # Reproject study area to UTM
        study_area_utm = self.study_area.to_crs(utm_crs)

        # Get bounds of the study area
        minx, miny, maxx, maxy = study_area_utm.total_bounds

        # Generate grid cells
        rows = int((maxy - miny) / cell_size)
        cols = int((maxx - minx) / cell_size)
        cells = []
        for i in range(cols):
            for j in range(rows):
                cells.append(Polygon([
                    (minx + i * cell_size, miny + j * cell_size),
                    (minx + (i + 1) * cell_size, miny + j * cell_size),
                    (minx + (i + 1) * cell_size, miny + (j + 1) * cell_size),
                    (minx + i * cell_size, miny + (j + 1) * cell_size)
                ]))

        # Create a GeoDataFrame
        grid = gpd.GeoDataFrame(geometry=cells, crs=utm_crs)

        # Reproject grid back to EPSG:4326
        self.grid = grid.to_crs('epsg:4326')

        # Collect HQ points
        hq_points = []
        if (self.hq_longitude_var.get() != 0.0 or self.hq_latitude_var.get() != 0.0):
            hq_points.append((float(self.hq_longitude_var.get()), float(self.hq_latitude_var.get()), float(self.hq_coefficient_var.get())))
        
        for entry in self.hq_entries:
            hq_longitude_var = entry[2]
            hq_latitude_var = entry[3]
            hq_coefficient_var = entry[4]
            hq_points.append((float(hq_longitude_var.get()), float(hq_latitude_var.get()), float(hq_coefficient_var.get())))

        def calculate_weighted_distance(row, hq_points):
            distances = []
            for lon, lat, coeff in hq_points:
                hq_point = Point(lon, lat)
                distance = row.geometry.centroid.distance(hq_point) / coeff
                distances.append(distance)
            return min(distances)

        self.grid['distance'] = self.grid.apply(lambda row: calculate_weighted_distance(row, hq_points), axis=1)







    def calculate_lengths(self):
        # Determine the UTM zone for the study area's centroid for more accurate distance calculations
        centroid = self.study_area.unary_union.centroid
        utm_zone = int((centroid.x + 180) / 6) + 1
        utm_crs = f'+proj=utm +zone={utm_zone} +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
        
        # Project grid and GPS tracks to UTM CRS for accurate length calculation
        grid_utm = self.grid.to_crs(utm_crs)
        gps_tracks_utm = self.gps_tracks.to_crs(utm_crs)

        # Initialize a column for lengths in the original grid dataframe
        self.grid['track_length'] = 0.0

        for index, cell in grid_utm.iterrows():
            total_length = 0
            cell_geom = cell.geometry
            for _, track in gps_tracks_utm.iterrows():
                track_geom = track.geometry
                if cell_geom.intersects(track_geom):
                    intersection = cell_geom.intersection(track_geom)
                    # Calculate length in meters
                    if isinstance(intersection, LineString):
                        total_length += intersection.length
                    elif isinstance(intersection, (MultiLineString, GeometryCollection)):
                        total_length += sum(line.length for line in intersection.geoms if isinstance(line, LineString))
            # Update the original grid dataframe with calculated lengths in meters
            self.grid.at[index, 'track_length'] = total_length

    def save_grid(self):
        output_path = filedialog.asksaveasfilename(defaultextension=".shp", filetypes=[("Shapefiles", "*.shp")])
        if output_path:
            self.grid.to_file(output_path)
            print("Grid saved to shapefile.")

    def process_data(self):
        self.update_progress("Starting process...")
        try:
            cell_size = int(self.grid_resolution_var.get())
        except ValueError:
            self.update_progress("Invalid grid resolution.")
            return

        if hasattr(self, 'study_area') and hasattr(self, 'gps_tracks'):
            self.generate_grid(cell_size)
            if not hasattr(self, 'grid') or self.grid.empty:
                self.update_progress("Grid generation failed or is empty.")
                return
            
            # Ensure the grid is clipped to the study area
            self.grid = gpd.clip(self.grid, self.study_area)
            self.calculate_lengths()

            # Calculate the mean patrol length in meters per cell
            mean_patrol_length = self.grid['track_length'].mean()

            # Assuming there is a baseline coverage value in the dynamically added rows
            for vars, widgets in self.dynamic_rows:
                if vars["indicator_type_var"].get() == "1. Patrol [m/cell]":
                    coverage_distance = vars["baseline_value_var"].get()
                    if coverage_distance is not None and coverage_distance > 0:
                        # Adjust 'track_length' by capping it at 'coverage_distance' (baseline value) before calculating 'coverage'
                        self.grid['capped_track_length'] = self.grid['track_length'].apply(lambda x: min(x, coverage_distance))

                        # Calculate 'coverage' as the ratio of capped track length to baseline, ensuring it's capped at 100%
                        self.grid['coverage'] = self.grid['capped_track_length'] / coverage_distance

                        # Calculate the mean of the 'coverage' values for the coverage score
                        mean_coverage = self.grid['coverage'].mean() * 100  # Convert to percentage for the score

                        # Prevent it from going beyond 100
                        if mean_coverage > 100:
                            mean_coverage = 100

                        # Update the UI: mean_patrol_length for the "Value" field, mean_coverage for the "Score" field
                        vars["value_var"].set(f"{mean_patrol_length:.2f}")
                        vars["score_var"].set(f"{mean_coverage:.2f}%")
                        self.update_progress(f"Mean patrol length in each grid cell: {mean_patrol_length:.2f} meters/cell")
                        self.update_progress(f"Mean coverage score: {mean_coverage:.2f}%")
                    else:
                        self.update_progress("Baseline coverage value not found or invalid.")
                    break
            
            self.visualize_grid()
        else:
            self.update_progress("Please load both a study area shapefile and a GPS track file.")




if __name__ == "__main__":
    root = tk.Tk()
    app = GridApp(root)
    root.mainloop()