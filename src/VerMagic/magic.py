import tkinter as tk
from tkinter import filedialog
import pandas as pd
from os import path
from plotnine import *
from numpy import nan

def open_file():
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename()
    return(file_path)

def new_paths():
    global newpath, newxlsx
    newpath = filedialog.askdirectory(title = 'Select an empty folder to store all files and plots.', initialdir = input_path)    
    newxlsx = path.join(newpath,path.splitext(path.basename(input_path))[0] + '_MAGIC.xlsx')

def save_excel(data: pd.DataFrame):
    data.to_excel(newxlsx)

def save_plot(plot, name: str):
    plot.save(path = newpath, filename = name + '.png', format = 'png')

def read_data():
    data = pd.read_excel(input_path)
    data = data.drop(columns=["flowcell","batch_name"]).drop_duplicates()
    data = data[~data["sample_barcode"].isna()]
    metrics_of_global_interest = ['NES','fetal_fraction', 'NCV_X', 'NCV_Y', 'number_of_cnv_events', 'non_excluded_sites']
    metrics_of_chr_interest = ['region_classification', 'region_llr_trisomy', 'region_llr_monosomy', 'region_t_stat_long_reads']
    data1 = data[data["metric_name"].isin(metrics_of_global_interest)]
    data2 = data[data["metric_name"].isin(metrics_of_chr_interest)]
    data1_pivot = data1.drop(columns='region').pivot(index='sample_barcode', columns="metric_name")
    data1_pivot.columns = data1_pivot.columns.droplevel()
    data2["variable"] = data2['metric_name'] + "-" + data2["region"]
    data2_pivot = data2.drop(columns=['metric_name','region']).pivot(index='sample_barcode', columns="variable")
    data2_pivot.columns = data2_pivot.columns.droplevel()
    data2_pivot = data2_pivot[data2_pivot.columns[(~data2_pivot.columns.str.startswith("region_classification")) | (data2_pivot.columns.str.endswith("X"))]]
    data_merged = data1_pivot.merge(data2_pivot, left_index = True, right_index = True)
    return(data, data_merged)

def plot_bars(data: pd.DataFrame):
    plot_data = data[data["metric_name"] == "non_excluded_sites"] 
    plot_data['metric_value'] = plot_data['metric_value'].astype(float)
    plot = ggplot(plot_data, aes(x="sample_barcode", y="metric_value")) +  geom_bar(stat = "identity") + theme(axis_text_x=element_text(rotation=90, hjust=1))
    save_plot(plot, name = "Bars")

def plot_ncvx_y(data: pd.DataFrame):
    plot_data = data[(data["metric_name"] == "NCV_X") | (data["metric_name"] == "NCV_Y")] 
    plot_data['metric_value'] = plot_data['metric_value'].astype(float)
    plot_data = plot_data.pivot_table(values=["metric_value"], columns="metric_name",index=["sample_barcode"])
    plot_data.columns = plot_data.columns.droplevel()
    plot = ggplot(plot_data, aes(x="NCV_X", y="NCV_Y")) +  geom_point() + theme(axis_text_x=element_text(rotation=90, hjust=1)) + theme_light()
    save_plot(plot, "NCVX_vs_NCVY")

def plot_ncvx_ff(data: pd.DataFrame):
    plot_data = data[(data["metric_name"] == "fetal_fraction") | (data["metric_name"] == "NCV_Y")] 
    plot_data['metric_value'] = plot_data['metric_value'].astype(float)
    plot_data = plot_data.pivot_table(values=["metric_value"], columns="metric_name",index=["sample_barcode"])
    plot_data.columns = plot_data.columns.droplevel()
    plot = ggplot(plot_data, aes(x="NCV_Y", y="fetal_fraction")) +  geom_point() + theme(axis_text_x=element_text(rotation=90, hjust=1)) + theme_light()
    save_plot(plot, "NCVX_vs_FetalFraction")

def plot_ncvx_ff_per_chr(data: pd.DataFrame):
    plot_data1 = data[data.metric_name == "fetal_fraction"]
    plot_data1.drop(columns="region", inplace= True)
    plot_data2 = data[(data.metric_name == "region_llr_trisomy") & ((data.region == "chr13") | (data.region == "chr18") | (data.region == "chr21"))]
    p1 = plot_data1.pivot_table(values=["metric_value"], columns=["metric_name"],index=["sample_barcode"]).droplevel(0, axis = 1)
    p2 = plot_data2.pivot_table(values=["metric_value"], columns=["metric_name"],index=["sample_barcode","region"]).reset_index().droplevel([1], axis=1)
    plot_data = p1.merge(p2, left_on = "sample_barcode", right_on = "sample_barcode", how = "outer")
    plot = ggplot(plot_data, aes(x="metric_value", y="fetal_fraction", colour = "region")) +  geom_point() + theme(axis_text_x=element_text(rotation=90, hjust=1)) + theme_light()
    save_plot(plot, "NCVX_vs_FetalFraction_per_chromosome")

if __name__ == "__main__":
    input_path = open_file()
    data, data_merged = read_data()
    new_paths()
    save_excel(data_merged)
    plot_bars(data)
    plot_ncvx_y(data)
    plot_ncvx_ff(data)
    plot_ncvx_ff_per_chr(data)