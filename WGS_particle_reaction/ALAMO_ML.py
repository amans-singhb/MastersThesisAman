from alamopy import almain as alamo
import numpy as np
import pandas as pd
import time

# Load data
raw_data = pd.read_csv('WGS_particle_reaction\\imp\\test_data_init_imp_no_t_10.csv', header=None, delimiter=',')

input_labels = ["T", "C_1", "C_2", "C_3", "C_4", "C_5", "C_c_1", "C_c_2", "C_c_3", "C_c_4", "C_c_5"]
output_labels = ["dC_c_1_dt", "dC_c_2_dt", "dC_c_3_dt", "dC_c_4_dt", "dC_c_5/dt"]

l_input = len(input_labels)


# print(raw_data)
print("Data shape: ", raw_data.shape)

data = raw_data.to_numpy()
inputs = data[:, 0:l_input]
outputs = data[:, l_input:]

dCc1_data = outputs[:, 0]
dCc2_data = outputs[:, 1]
dCc3_data = outputs[:, 2]
dCc4_data = outputs[:, 3]
dCc5_data = outputs[:, 4]

print("Data dcc1: ", dCc1_data)
print("Data dcc2: ", dCc2_data)
print("Data dcc3: ", dCc3_data)
print("Data dcc4: ", dCc4_data)
print("Data dcc5: ", dCc5_data)


print("Inputs shape: ", inputs.shape)
print("Outputs shape: ", outputs.shape)
print("Inputs: ", inputs)
print("Outputs: ", outputs)

print(data[0, :])
print(inputs[0, :])
print(outputs[0, :])

# ALAMO
start = time.time()

#opts = alamo.doalamo(inputs, outputs, xlabels = input_labels, zlabels = output_labels, linfcns = 1, expfcns = 1, logfcns = 1, sinfcns = 1, cosfns = 1, monomialpower = [1,2,3,4], multi2power = [1,2,3,4], ratiopower =[1,2],  keep_alm_file=True, keep_lst_file=True, print_alm_output=True)

opts_dcc1 = alamo.doalamo(inputs, dCc1_data, xlabels = input_labels, zlabels = ["d_C_c_1/dt"], linfcns = 1, expfcns = 1, logfcns = 1, sinfcns = 1, cosfns = 1, monomialpower = [1,2,3,4], multi2power = [1,2,3,4], ratiopower =[1,2],  keep_alm_file=True, keep_lst_file=True, print_alm_output=False)
opts_dcc2 = alamo.doalamo(inputs, dCc2_data, xlabels = input_labels, zlabels = ["d_C_c_2/dt"], linfcns = 1, expfcns = 1, logfcns = 1, sinfcns = 1, cosfns = 1, monomialpower = [1,2,3,4], multi2power = [1,2,3,4], ratiopower =[1,2],  keep_alm_file=True, keep_lst_file=True, print_alm_output=False)
opts_dcc3 = alamo.doalamo(inputs, dCc3_data, xlabels = input_labels, zlabels = ["d_C_c_3/dt"], linfcns = 1, expfcns = 1, logfcns = 1, sinfcns = 1, cosfns = 1, monomialpower = [1,2,3,4], multi2power = [1,2,3,4], ratiopower =[1,2],  keep_alm_file=True, keep_lst_file=True, print_alm_output=False)
opts_dcc4 = alamo.doalamo(inputs, dCc4_data, xlabels = input_labels, zlabels = ["d_C_c_4/dt"], linfcns = 1, expfcns = 1, logfcns = 1, sinfcns = 1, cosfns = 1, monomialpower = [1,2,3,4], multi2power = [1,2,3,4], ratiopower =[1,2],  keep_alm_file=True, keep_lst_file=True, print_alm_output=False)
opts_dcc5 = alamo.doalamo(inputs, dCc5_data, xlabels = input_labels, zlabels = ["d_C_c_5/dt"], linfcns = 1, expfcns = 1, logfcns = 1, sinfcns = 1, cosfns = 1, monomialpower = [1,2,3,4], multi2power = [1,2,3,4], ratiopower =[1,2],  keep_alm_file=True, keep_lst_file=True, print_alm_output=False)


end = time.time()

# opts["out"]["z"]
#print(opts["out"]["dC_c_1_dt"]["model_str"])

print("Time: ", end - start)