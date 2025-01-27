# -*- coding: utf-8 -*-
"""
Created on Sat Jan 18 12:40:58 2025

@author: Nemat002
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import re
import subprocess
import os
import sys
import time

def save_dict_to_file(filename, data_dict):
    """
    Saves the keys and values of a dictionary to a file with keys in the first column,
    values in subsequent columns, and column indices at the top for existing columns only.
    
    If the file exists, appends the values as new columns in the same row.
    If the file does not exist, creates it and writes the keys and column indices.
    
    Parameters:
    - filename: str, name of the file to write to.
    - data_dict: dict, dictionary with keys and float values.
    """
    file_exists = os.path.exists(filename)
    column_width = 20
    
    if not file_exists:
        with open(filename, 'w') as file:
            # Write initial column indices
            file.write(' ' * column_width)  # Indent to align with keys
            
            for i in range(1, 2):  # First set of values is the first column
                file.write(f"{i:<{column_width}}")
            file.write('\n')
            
            # Write keys in the second row
            for key in data_dict.keys():
                file.write(f"{key:<{column_width}}\n")

    with open(filename, 'r+') as file:
        lines = file.readlines()
        file.seek(0)
        
        # Count existing columns by counting entries in the first line
        num_existing_columns = len(lines[0].split())
        new_index = num_existing_columns + 1
        
        # Add the new column index
        lines[0] = lines[0].rstrip('\n') + f"{new_index:<{column_width}}\n"
        
        # Append new values under their respective keys
        for i, (key, value) in enumerate(data_dict.items()):
            formatted_value = f"{value:.8f}" if isinstance(value, float) else str(value)
            lines[i + 1] = lines[i + 1].rstrip('\n') + f"{formatted_value:<{column_width}}\n"
        
        file.writelines(lines)

def read_file_to_dict(filename):
    """
    Reads a text file with keys in the first column, column indices on top, and values in subsequent columns.
    Returns a dictionary with keys from the first column and values from the last column.
    
    Parameters:
    - filename: str, name of the file to read.
    
    Returns:
    - dict: dictionary with keys from the first column and values from the last column.
    """
    result_dict = {}
    
    with open(filename, 'r') as file:
        lines = file.readlines()
        
        # Skip the first line with column indices
        for line in lines[1:]:
            columns = line.split()
            key = columns[0]
            value = float(columns[-1]) if '.' in columns[-1] else int(columns[-1])
            result_dict[key] = value
    
    return result_dict
        
def indices_finder(string, begin_char, end_char):
    
    start_index = -1 + len(begin_char)
    
    while 1:
        start_index += 1
        
        if string[start_index-len(begin_char):start_index] == begin_char:
            
            end_index = string.index(end_char, start_index)
            
            try:
                value = float(string[start_index:end_index])
                break
            except:
                continue
        
    
    return start_index, end_index

def read_file_to_dict_list(file_path):
    """
    Reads a file and converts it into a list of dictionaries.
    The first row of the file is used as the keys for the dictionaries.

    :param file_path: Path to the text file
    :return: List of dictionaries containing the file's data
    """
    dict_list = []
    total_quantities = []
    updating_quantities = []
    
    with open(file_path, 'r') as file:
        # Read the first line and split it into keys
        keys = file.readline().strip().split()
        
        # Iterate over the remaining lines
        for line in file:
            # Split the line into values
            values = line.strip().split()
            
            # Create a dictionary by zipping keys with values
            row_dict = dict(zip(keys, values))
            
            # Convert 'update_switch' and 'factor' to integers
            row_dict['update_switch'] = int(row_dict['update_switch'])
            row_dict['factor'] = int(row_dict['factor'])
            
            if row_dict['end_char']=='\\n':
                row_dict['end_char'] ='\n'
            
            if row_dict['nondim_quantity'] not in total_quantities:
                total_quantities.append(row_dict['nondim_quantity'])
                
            if row_dict['update_switch']==1:
                # Append the dictionary to the list
                dict_list.append(row_dict)
                if (row_dict['nondim_quantity'] not in updating_quantities):
                        updating_quantities.append(row_dict['nondim_quantity'])
                    
    return dict_list, total_quantities, updating_quantities

def read_SGD_params(file_path):
    """
    Reads a file and assigns values to separate variables based on the file content.
    
    :param file_path: Path to the configuration file
    :return: Tuple containing max_iter, l_rate, and stop_thresh
    """
    with open(file_path, 'r') as file:
        for line in file:
            # Split the line into name and value
            name, value = line.strip().split()
            
            # Assign values to variables based on the name
            if name == 'max_iter':
                max_iter = int(value)
            elif name == 'l_rate':
                l_rate = float(value)
            elif name == 'stop_thresh':
                stop_thresh = float(value)
            elif name == 'dh':
                dh = float(value)
    
    return max_iter, l_rate, stop_thresh, dh

def read_init_vals(list_of_modif_tasks, updating_quantities, file_path):
    
    vals_init = {key: 0.0 for key in updating_quantities}
    
    read_switches = {key: 0 for key in updating_quantities}
        
    for task_info in list_of_modif_tasks:

        with open(file_path, 'r') as file:
            for line in file:
                if (line[:len(task_info['line_word'])] in  task_info['line_word']) and  read_switches[task_info['nondim_quantity']] == 0:
                    
                    # start_index = line.index(task_info['begin_char']) + len(task_info['begin_char'])
                    # end_index = line.index(task_info['end_char'], start_index)
                    
                    start_index, end_index = indices_finder(line, task_info['begin_char'], task_info['end_char'])
                    
                    vals_init[task_info['nondim_quantity']] = task_info['factor']*float(line[start_index:end_index])
                    read_switches[task_info['nondim_quantity']] = 1

    return vals_init

def modify_func(task_info, file_path, new_val):
    
    modified_lines = []
    
    with open(file_path, 'r') as file:
        for line in file:
            if line[:len(task_info['line_word'])] == task_info['line_word']:
                
                # start_index = line.index(task_info['begin_char']) + len(task_info['begin_char'])
                # end_index = line.index(task_info['end_char'])
                
                start_index, end_index = indices_finder(line, task_info['begin_char'], task_info['end_char'])
                
                modified_line = line[:start_index]+str(task_info['factor']*new_val)+line[end_index:]
                modified_lines.append(modified_line)
                
            else:
                modified_lines.append(line)
    
    with open(file_path, 'w') as file:
        file.writelines(modified_lines)
       
    return 1


def indiv_type_optimize(optimization_type, dict_of_folders):
    
    os.chdir(dict_of_folders['main'][0])
    
    list_of_modif_tasks,\
        total_quantities,\
            updating_quantities\
                = read_file_to_dict_list("cost_opt/switches_addresses.txt")

    n_dim = len(updating_quantities)

    max_iter, l_rate, stop_thresh, dh = read_SGD_params("cost_opt/SGD_params.txt")

    vals_init = read_init_vals(list_of_modif_tasks, updating_quantities, "source/params.csv")
    cost_init = read_file_to_dict("cost_opt/cost_hist.txt")

    save_dict_to_file("cost_opt/quan_hist.txt", vals_init)
    
    os.chdir(mother_directory)
    
    new_vals = vals_init.copy()
    gradient = vals_init.copy()
    # for key in new_vals.keys():
    #     new_vals[key] = new_vals[key]*2

        
    plus_Dh_vals  = new_vals.copy()
    minus_Dh_vals = new_vals.copy()

    iteration_c = 0
    cost_var_frac = 1.0
    while((iteration_c<max_iter) and (cost_var_frac > stop_thresh) ):
        
        os.chdir(dict_of_folders['main'][0])
        
        cost_old = read_file_to_dict("cost_opt/cost_hist.txt")
        
        ############ finding two points (+dh, -dh) ##############
        while_cond = True
        while while_cond:
            # randomize_a_vector()
            rand_vector = np.random.randn(n_dim)
            # Step 2: Normalize the vector to make it unit length
            magnitude = np.linalg.norm(rand_vector)
            unit_vector = rand_vector / magnitude
            # unit_vector_inv = 1.0 / unit_vector
            
            u_vec_dict = dict()
            # u_vec_inv_dict = dict()
            for key_c in range(len(updating_quantities)):
                key = updating_quantities[key_c]
                u_vec_dict[key]     = unit_vector[key_c]
                # u_vec_inv_dict[key] = unit_vector_inv[key_c]
                
                plus_Dh_vals[key]  = new_vals[key] + dh * u_vec_dict[key]
                minus_Dh_vals[key] = new_vals[key] - dh * u_vec_dict[key]
            
            if (plus_Dh_vals['psi_w'] > plus_Dh_vals['psi_g0'] \
                and \
            	minus_Dh_vals['psi_w'] > minus_Dh_vals['psi_g0']):
                
                while_cond = 0
        ############ finding two points (+dh, -dh) ##############
        os.chdir(mother_directory)
        
        
        ############ Gradient Estimation #############
        
        # do_simulations_for_+Dh
        for modif_task_c in range(len(list_of_modif_tasks)):
            modif_task = list_of_modif_tasks[modif_task_c]
            modify_func(modif_task, dict_of_folders['plus_dh'][0] +"/"+"source/params.csv", plus_Dh_vals[modif_task['nondim_quantity']])

        os.chdir(dict_of_folders['plus_dh'][0])
        N_runs_plus_dh = np.loadtxt("N_runs.csv", delimiter=',', dtype=int)
        os.system('./clear_copy_init_run_pp.sh')
        os.system('./terminal_opener.sh')
        
        print("****************************************")
        print("Optimization type: "+str(optimization_type))
        print("Now doing : "+dict_of_folders['plus_dh'][0])
        print("iteration_c : " + str(iteration_c))
        # print("cost_var_frac : " + str(cost_var_frac))
        print("****************************************")
        
        switch_still_running = 1
        while switch_still_running: # THis makes sure that all runs are finished
            time.sleep(60)
            for run_c in range(N_runs_plus_dh):
                if os.path.exists("run_"+str(run_c+1)+"/pp_data/time.txt"):
                    switch_still_running = 0
                else:
                    switch_still_running = 1
                    break
                
        time.sleep(20)
        os.system('./pp_plus_cost.sh')
        cost_plus_Dh = read_file_to_dict("cost_opt/cost_hist.txt")
        # do_simulations_for_+Dh
        
        os.chdir(mother_directory) # get back to mother directory
        
        # do_simulations_for_-Dh
        for modif_task_c in range(len(list_of_modif_tasks)):
            modif_task = list_of_modif_tasks[modif_task_c]
            modify_func(modif_task, dict_of_folders['minus_dh'][0] +"/"+"source/params.csv", minus_Dh_vals[modif_task['nondim_quantity']])

        os.chdir(dict_of_folders['minus_dh'][0])
        N_runs_minus_dh = np.loadtxt("N_runs.csv", delimiter=',', dtype=int)
        os.system('./clear_copy_init_run_pp.sh')
        os.system('./terminal_opener.sh')
        
        print("****************************************")
        print("Optimization type: "+str(optimization_type))
        print("Now doing : "+dict_of_folders['minus_dh'][0])
        print("iteration_c : " + str(iteration_c))
        # print("cost_var_frac : " + str(cost_var_frac))
        print("****************************************")
        
        switch_still_running = 1
        while switch_still_running: # THis makes sure that all runs are finished
            time.sleep(60)
            for run_c in range(N_runs_minus_dh):
                if os.path.exists("run_"+str(run_c+1)+"/pp_data/time.txt"):
                    switch_still_running = 0
                else:
                    switch_still_running = 1
                    break
        time.sleep(20)
        os.system('./pp_plus_cost.sh')
        cost_minus_Dh = read_file_to_dict("cost_opt/cost_hist.txt")
        # do_simulations_for_-Dh
        
        os.chdir(mother_directory) # get back to mother directory
        os.chdir(dict_of_folders['main'][0])
        
        directional_grad = (cost_plus_Dh['tot'] - cost_minus_Dh['tot']) / (2*dh)
        
        for key in gradient.keys():
            # gradient[key] = directional_grad * unit_vector_inv[key]
            gradient[key] = directional_grad * u_vec_dict[key]
        save_dict_to_file("cost_opt/grad_hist.txt", gradient)
        ############ Gradient Estimation #############
        
        
        
        ############ Finding new_vals #############
        for key in new_vals.keys():
            new_vals[key] = new_vals[key] - l_rate * gradient[key]
        save_dict_to_file("cost_opt/quan_hist.txt", new_vals)
        ############ Finding new_vals #############
        
        ############ Do the simulation with new_vals #############
        for modif_task_c in range(len(list_of_modif_tasks)):
            modif_task = list_of_modif_tasks[modif_task_c]
            modify_func(modif_task, "source/params.csv", new_vals[modif_task['nondim_quantity']])
        
        N_runs_main = np.loadtxt("N_runs.csv", delimiter=',', dtype=int)
        os.system('./clear_copy_init_run_pp.sh')
        os.system('./terminal_opener.sh')
        
        print("****************************************")
        print("Optimization type: "+str(optimization_type))
        print("Now doing : "+dict_of_folders['main'][0])
        print("iteration_c : " + str(iteration_c))
        # print("cost_var_frac : " + str(cost_var_frac))
        print("****************************************")
        
        switch_still_running = 1
        while switch_still_running: # THis makes sure that all runs are finished
            time.sleep(60)
            for run_c in range(N_runs_main):
                if os.path.exists("run_"+str(run_c+1)+"/pp_data/time.txt"):
                    switch_still_running = 0
                else:
                    switch_still_running = 1
                    break
        time.sleep(20)
        os.system('./pp_plus_cost.sh')
        ############ Do the simulation with new_vals #############
        
        cost_new = read_file_to_dict("cost_opt/cost_hist.txt")
        
        iteration_c += 1
        cost_var_frac = np.abs(cost_new['tot']-cost_old['tot'])/cost_old['tot']
        
        print("****************************************")
        print("Optimization type: "+str(optimization_type))
        print("iteration_c : " + str(iteration_c)+" done!")
        print("cost_var_frac : " + str(cost_var_frac))
        print("****************************************")
        os.chdir(mother_directory)

    return

def overal_type_optimize(dict_of_folders):
    
    os.chdir(mother_directory)
    
    list_of_modif_tasks,\
        total_quantities,\
            updating_quantities\
                = read_file_to_dict_list("overal_cost_opt/switches_addresses.txt")
    
    n_dim = len(updating_quantities)
    max_iter, l_rate, stop_thresh, dh = read_SGD_params("overal_cost_opt/SGD_params.txt")
    
    vals_init = dict()
    cost_init = read_file_to_dict(dict_of_folders['main'][0]+"/cost_opt/cost_hist.txt")
    cost_init = {key: 0.0 for key in cost_init.keys()} # making zero first
    
    for folder_name in dict_of_folders['main']:
        
        os.chdir(folder_name)#############
        
        list_of_modif_tasks_part,\
            total_quantities_part,\
                updating_quantities_part\
                    = read_file_to_dict_list("cost_opt/switches_addresses.txt")
        
        vals_init_part = read_init_vals(list_of_modif_tasks_part, updating_quantities_part, "source/params.csv")
        vals_init.update(vals_init_part)
        
        cost_init_part = read_file_to_dict("cost_opt/cost_hist.txt")
        cost_init = {key: cost_init[key] + cost_init_part[key] for key in cost_init.keys()}
        
        vals_init_part.clear()
        cost_init_part.clear()
        list_of_modif_tasks_part.clear()
        total_quantities_part.clear()
        updating_quantities_part.clear()
        os.chdir(mother_directory)########
    

    

    # vals_init = read_init_vals(list_of_modif_tasks, updating_quantities, "source/params.csv")
    # cost_init = read_file_to_dict("cost_opt/cost_hist.txt")

    save_dict_to_file("overal_cost_opt/quan_hist.txt", vals_init)
    save_dict_to_file("overal_cost_opt/cost_hist.txt", cost_init)
    
    new_vals = vals_init.copy()
    gradient = vals_init.copy()
    # for key in new_vals.keys():
    #     new_vals[key] = new_vals[key]*2

    
    plus_Dh_vals  = new_vals.copy()
    minus_Dh_vals = new_vals.copy()

    iteration_c = 0
    cost_var_frac = 1.0
    while((iteration_c<max_iter) and (cost_var_frac > stop_thresh) ):
        
        os.chdir(mother_directory)
        
        cost_old = read_file_to_dict("overal_cost_opt/cost_hist.txt")

        ############ finding two points (+dh, -dh) ##############
        while_cond = True
        while while_cond:
            # randomize_a_vector()
            rand_vector = np.random.randn(n_dim)
            # Step 2: Normalize the vector to make it unit length
            magnitude = np.linalg.norm(rand_vector)
            unit_vector = rand_vector / magnitude
            # unit_vector_inv = 1.0 / unit_vector
            
            u_vec_dict = dict()
            # u_vec_inv_dict = dict()
            for key_c in range(len(updating_quantities)):
                key = updating_quantities[key_c]
                u_vec_dict[key]     = unit_vector[key_c]
                # u_vec_inv_dict[key] = unit_vector_inv[key_c]
                
                plus_Dh_vals[key]  = new_vals[key] + dh * u_vec_dict[key]
                minus_Dh_vals[key] = new_vals[key] - dh * u_vec_dict[key]
            
            if (plus_Dh_vals['psi_w'] > plus_Dh_vals['psi_g0'] \
                and \
            	minus_Dh_vals['psi_w'] > minus_Dh_vals['psi_g0']):
                
                while_cond = 0
        ############ finding two points (+dh, -dh) ##############
        
        # os.chdir(dict_of_folders['main'][0])
        # os.chdir(mother_directory)
        
        
        ############ Gradient Estimation #############
        
        # do_simulations_for_+Dh
        for folder_name in dict_of_folders['plus_dh']:
            
            os.chdir(folder_name)#############
            
            for modif_task_c in range(len(list_of_modif_tasks)):
                modif_task = list_of_modif_tasks[modif_task_c]
                modify_func(modif_task, "source/params.csv", plus_Dh_vals[modif_task['nondim_quantity']])
    
            N_runs_plus_dh = np.loadtxt("N_runs.csv", delimiter=',', dtype=int)
            os.system('./clear_copy_init_run_pp.sh')
            os.system('./terminal_opener.sh')
            
            print("****************************************")
            print("Optimization type: "+str(optimization_type))
            print("Now doing : "+folder_name)
            print("iteration_c : " + str(iteration_c))
            # print("cost_var_frac : " + str(cost_var_frac))
            print("****************************************")
            
            switch_still_running = 1
            while switch_still_running: # THis makes sure that all runs are finished
                time.sleep(60)
                for run_c in range(N_runs_plus_dh):
                    try:
                        np.loadtxt("run_"+str(run_c)+"/pp_data/time.txt", delimiter=',')
                    except:
                        switch_still_running = 1
                        break
                    switch_still_running = 0
            time.sleep(20)
            os.system('./pp_plus_cost.sh')
            cost_plus_Dh = read_file_to_dict("cost_opt/cost_hist.txt")
        # do_simulations_for_+Dh
        
        os.chdir(mother_directory) # get back to mother directory
        
        # do_simulations_for_-Dh
        for folder_name in dict_of_folders['minus_dh']:
            
            os.chdir(folder_name)#############
            
            for modif_task_c in range(len(list_of_modif_tasks)):
                modif_task = list_of_modif_tasks[modif_task_c]
                modify_func(modif_task, "source/params.csv", minus_Dh_vals[modif_task['nondim_quantity']])
    
            N_runs_minus_dh = np.loadtxt("N_runs.csv", delimiter=',', dtype=int)
            os.system('./clear_copy_init_run_pp.sh')
            os.system('./terminal_opener.sh')
            
            print("****************************************")
            print("Optimization type: "+str(optimization_type))
            print("Now doing : "+folder_name)
            print("iteration_c : " + str(iteration_c))
            # print("cost_var_frac : " + str(cost_var_frac))
            print("****************************************")
            
            switch_still_running = 1
            while switch_still_running: # THis makes sure that all runs are finished
                time.sleep(60)
                for run_c in range(N_runs_minus_dh):
                    try:
                        np.loadtxt("run_"+str(run_c)+"/pp_data/time.txt", delimiter=',')
                    except:
                        switch_still_running = 1
                        break
                    switch_still_running = 0
            time.sleep(20)
            os.system('./pp_plus_cost.sh')
            cost_minus_Dh = read_file_to_dict("cost_opt/cost_hist.txt")
        # do_simulations_for_-Dh
        
        os.chdir(mother_directory) # get back to mother directory
        
        directional_grad = (cost_plus_Dh['tot'] - cost_minus_Dh['tot']) / (2*dh)
        
        for key in gradient.keys():
            # gradient[key] = directional_grad * unit_vector_inv[key]
            gradient[key] = directional_grad * u_vec_dict[key]
        save_dict_to_file("overal_cost_opt/grad_hist.txt", gradient)
        ############ Gradient Estimation #############
        
        
        ############ Finding new_vals #############
        for key in new_vals.keys():
            new_vals[key] = new_vals[key] - l_rate * gradient[key]
        save_dict_to_file("overal_cost_opt/quan_hist.txt", new_vals)
        ############ Finding new_vals #############
        
        ############ Do the simulation with new_vals #############
        cost_new = {key: 0.0 for key in cost_old.keys()}
        
        for folder_name in dict_of_folders['main']:
            
            os.chdir(folder_name)#############
            
            for modif_task_c in range(len(list_of_modif_tasks)):
                modif_task = list_of_modif_tasks[modif_task_c]
                modify_func(modif_task, "source/params.csv", new_vals[modif_task['nondim_quantity']])
            
            N_runs_main = np.loadtxt("N_runs.csv", delimiter=',', dtype=int)
            os.system('./clear_copy_init_run_pp.sh')
            os.system('./terminal_opener.sh')
            
            print("****************************************")
            print("Optimization type: "+str(optimization_type))
            print("Now doing : "+folder_name)
            print("iteration_c : " + str(iteration_c))
            # print("cost_var_frac : " + str(cost_var_frac))
            print("****************************************")
            
            switch_still_running = 1
            while switch_still_running: # THis makes sure that all runs are finished
                time.sleep(60)
                for run_c in range(N_runs_main):
                    try:
                        np.loadtxt("run_"+str(run_c)+"/pp_data/time.txt", delimiter=',')
                    except:
                        switch_still_running = 1
                        break
                    switch_still_running = 0
            time.sleep(20)
            os.system('./pp_plus_cost.sh')
            
            cost_new_part = read_file_to_dict("cost_opt/cost_hist.txt")
            cost_new = {key: cost_new[key] + cost_new_part[key] for key in cost_new.keys()}
            
            os.chdir(mother_directory)#############
        ############ Do the simulation with new_vals #############
        
        save_dict_to_file("overal_cost_opt/cost_hist.txt", cost_new)
        
        iteration_c += 1
        cost_var_frac = np.abs(cost_new['tot']-cost_old['tot'])/cost_old['tot']
        
        
        print("****************************************")
        print("Optimization type: "+str(optimization_type))
        print("iteration_c : " + str(iteration_c)+ " done!")
        print("cost_var_frac : " + str(cost_var_frac))
        print("****************************************")
        os.chdir(mother_directory)
    
    return

optimization_type = sys.argv[1]

mother_directory = os.getcwd()
os.chdir(mother_directory)

dict_of_folders = dict()

# Conditional logic to print based on the input string
if optimization_type == "WT":
    dict_of_folders['main'] = ['main/WT']
    dict_of_folders['plus_dh'] = ['plus_dh/WT']
    dict_of_folders['minus_dh'] = ['minus_dh/WT']
    
    indiv_type_optimize(optimization_type, dict_of_folders)
    
elif optimization_type == "C":
    dict_of_folders['main'] = ['main/C']
    dict_of_folders['plus_dh'] = ['plus_dh/C']
    dict_of_folders['minus_dh'] = ['minus_dh/C']
    
    indiv_type_optimize(optimization_type, dict_of_folders)
    
elif optimization_type == "mix":
    dict_of_folders['main'] = ['main/mix']
    dict_of_folders['plus_dh'] = ['plus_dh/mix']
    dict_of_folders['minus_dh'] = ['minus_dh/mix']
    
    indiv_type_optimize(optimization_type, dict_of_folders)

elif optimization_type == "overal":
    dict_of_folders['main'] = ['main/WT', 'main/C', 'main/mix']
    dict_of_folders['plus_dh'] = ['plus_dh/WT', 'plus_dh/C', 'plus_dh/mix']
    dict_of_folders['minus_dh'] = ['minus_dh/WT', 'minus_dh/C', 'minus_dh/mix']
    
    overal_type_optimize(dict_of_folders)


    

