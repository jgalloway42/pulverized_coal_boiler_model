import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
from itertools import chain
import os
from datetime import datetime
import joblib
import logging

'''General Helpers'''
def search_columns(search_str,df):
    return list(filter(lambda x: search_str.upper() in x.upper(),[y for y in df.columns]))

def filter_list(search_str,searchable_list):
    return list(filter(lambda x: search_str.upper() in x.upper(),[y for y in searchable_list]))

def flatten_list(list_of_lists):
    return list(chain(*list_of_lists))

'''Plotting Functions'''
def plot_df(data, cols_list, save_to_path = None,
             figsize=(15,20), linestyle='none', marker=','):
    if save_to_path:
        if os.path.exists(save_to_path):
            plt.imshow(plt.imread(save_to_path), aspect='auto',interpolation='nearest')
            return f'Read From File: {save_to_path}'
        else:
            _ = data[cols_list].plot(subplots=True,linestyle=linestyle,
                                     marker=marker,figsize=figsize)
            plt.savefig(save_to_path)
            return f'Saved to File: {save_to_path}'
    else:
        _ = data[cols_list].plot(subplots=True,linestyle=linestyle,marker=',',figsize=(15,30))
        return 'graph not saved to file...'

def plotly_graph(df, cols_to_graph, left_legend = False, save_to_path = None):
    f = px.line(data_frame=df,x=df.index,y=cols_to_graph)
    if left_legend:
        f.update_layout(legend=dict(yanchor="top",y=0.99,xanchor="left",x=0.01))
    if save_to_path:
        print(f'...saving graph to {save_to_path}')
        f.write_html(save_to_path)
    else:
        f.show()
    
'''Logger and Saves'''
def save_joblib(object ,folder_path, file_name, add_timestamp = False):
    if add_timestamp:
        current_timestamp = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
        root, ext = os.path.splitext(file_name)
        file_name = root + '_' + current_timestamp + ext # ext retains the .xxxx including the .
    joblib.dump(object,os.path.join(folder_path,file_name))
    print(f'File Saved: {os.path.join(folder_path,file_name)}')

def get_logger(folder_path ,file_name, logging_level=logging.DEBUG, add_timestamp = True):
    if add_timestamp:
        current_timestamp = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
        root, ext = os.path.splitext(file_name)
        file_name = root + '_' + current_timestamp + ext # ext retains the .xxxx including the .
    
    logger = logging.getLogger('Jupyter Notebook')
    logger.setLevel(logging_level)
    handler = logging.FileHandler(os.path.join(folder_path,file_name))
    handler.setLevel(logging_level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    print(f'Log Started: {os.path.join(folder_path,file_name)}')
    return  logger


