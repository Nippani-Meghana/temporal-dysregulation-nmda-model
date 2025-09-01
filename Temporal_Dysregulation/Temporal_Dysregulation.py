#MODELING TEMPORAL DYSREGULATION IN SCHIZOPHRENIA USING THE MORRIS-LECAR NEURON MODEL

import numpy as np
import matplotlib.pyplot as plt
np.random.seed(0)

# Parameters
C = 1.0
V_ca = 1
V_k = -0.7
V_L = -0.5
g_k = 2.0
g_ca = 1.0
g_L = 0.5
rate_const = 1.0/3.0
v0 = -0.5
v1 = -0.01
v2 = 0.15
v3 = 0.1
v4 = 0.145
I_app = 0.3
I_ca = 0.0
I_k = 0.0
I_L = 0.0 
V_th = -0.40
k = 0.06
w0 = 0.5 * (1 + np.tanh((v0-v3)/v4))
s0 = 1/(1 + np.exp(-(v0 - V_th)/k))
tau_nmda = 100
tau_gaba = 10.0
E_exc = 0.0
E_inh = -0.8
V_th_nmda = -0.4
dt = 0.01
V_exc = 0.0
V_inh = -0.8
V_E = 0.0
V_I = 0.0
w_e = w0
w_i = w0


#1. MODEL DEFINITION

def morris_lecar_model(v_e,w_e,s_E,v_i,w_i,s_I,params,W_IE = 1.0):
    g_nmda_local = params['g_nmda']
    g_gaba_local = params['g_gaba']
    W_EI = params['W_EI']
    I_app_e = params['I_app_e']
    I_app_i = params['I_app_i']
    V_th = params['V_th_gaba']
    #Excitatory stuff
    m_inf_e = 0.5 * (1 + np.tanh((v_e - v1)/v2))
    w_inf_e = 0.5 * (1 + np.tanh((v_e-v3)/v4))
    tau_w_e = 1/(np.cosh((v_e-v3)/(2.0*v4)))
    I_ca_e = g_ca * m_inf_e * (v_e - V_ca)
    I_k_e = g_k * w_e * (v_e - V_k)
    I_L_e = g_L * (v_e - V_L)

    #Inhibitory Stuff
    m_inf_i = 0.5 * (1 + np.tanh((v_i - v1)/v2))
    w_inf_i = 0.5 * (1 + np.tanh((v_i-v3)/v4))
    tau_w_i = 1/(np.cosh((v_i-v3)/(2.0*v4)))
    I_ca_i = g_ca * m_inf_i * (v_i - V_ca)
    I_k_i = g_k * w_i * (v_i - V_k)
    I_L_i = g_L * (v_i - V_L)

    s_inf_e = 1/(1 + np.exp(-(v_e - V_th_nmda)/k))
    s_e_next = s_E + dt * ((s_inf_e - s_E) / tau_nmda)

    s_inf_i = 1/(1 + np.exp(-(v_i - V_th)/k))   
    s_i_next = s_I + dt * ((s_inf_i - s_I) / tau_gaba)

    S_NMDA_on_I = W_IE * s_e_next   
    S_GABA_on_E = W_EI * s_i_next

    I_nmda_on_I = g_nmda_local* (v_i - E_exc) * S_NMDA_on_I
    I_gaba_on_E = g_gaba_local * (v_e - E_inh) * S_GABA_on_E
    I_syn_E = I_gaba_on_E
    I_syn_I = I_nmda_on_I

    dv_e = (I_app_e - (I_ca_e + I_k_e + I_L_e + I_syn_E)) / C
    dv_i = (I_app_i - (I_ca_i + I_k_i + I_L_i + I_syn_I)) / C
    dw_e = rate_const * (w_inf_e - w_e) / tau_w_e
    dw_i = rate_const * (w_inf_i - w_i) / tau_w_i

    v_e_next = v_e + dt*dv_e
    v_i_next = v_i + dt*dv_i
    w_e_next = w_e + dt*dw_e
    w_i_next = w_i + dt*dw_i

    return v_e_next,w_e_next,s_e_next,v_i_next,w_i_next,s_i_next

#2. RUNNING THE MODEL

def call_method(parameters,noise_level):
    time = 300
    dt = 0.01
    s_E = s0
    s_I = s0
    v_i = v0
    w_i = w0
    v_e = v0
    w_e = w0
    voltages_exc = []
    voltages_exc.append(v0)
    voltages_inh = []
    spike_timings_e = []
    spike_timings_i = []
    times = []
    refractory_period = 8.0
    last_spike_e = -np.inf
    last_spike_i = -np.inf
    for t in np.arange(0,time,dt):
        v_e += noise_level * np.random.randn()
        v_i += noise_level * np.random.randn()
        v_e,w_e,s_E,v_i,w_i,s_I = morris_lecar_model(v_e,w_e,s_E,v_i,w_i,s_I, parameters)
        voltages_exc.append(v_e)
        voltages_inh.append(v_i)
        times.append(t)
        

        if v_e>=V_th and (t - last_spike_e)>refractory_period:
            spike_timings_e.append(t)
            last_spike_e = t

        if v_i>=V_th and (t - last_spike_i)>refractory_period:
            spike_timings_i.append(t)
            last_spike_i = t

    return voltages_exc, times, spike_timings_e, voltages_inh,spike_timings_i


#3. PLOTTING THE MODEL

def plot_spikes(ve,ste,vi,sti):
    print("Excitatory Neuron Spikes at (ms):", np.round(ste, 2))
    print("The number of spikes of Excitatory Neuron is : ",len(ste))
    print("Inhibitory Neuron Spikes at (ms):", np.round(sti, 2))
    print("The number of spikes of Inhibitory Neuron is : ",len(sti))
    t_e = np.arange(len(ve)) * dt
    t_i = np.arange(len(vi)) * dt
    plt.xlabel("Time (ms)")
    plt.plot(t_e, ve, label = "Excitatory Neuron", color = 'crimson')
    plt.plot(t_i,vi, label = "Inhibitory Neuron", color = 'midnightblue')
    plt.ylabel("Membrane Potential (a.u)")
    plt.title("Morris Lecar Model")
    plt.legend()
    plt.show()

#4. NORMAL SIMULATION 
normal_parameters = {
    'g_nmda': 0.18,
    'g_gaba': 0.5,
    'W_EI': 0.3,
    'I_app_e': 0.3,
    'I_app_i': 0.09,
    'V_th_gaba': -0.4
}


#5. SCHIZOPHRENIA SIMULATION
Schizophrenia_parameters = {
    'g_nmda': 0.08,      
    'g_gaba': 0.45,       
    'W_EI': 0.5,        
    'I_app_e': 0.2,     
    'I_app_i': 0.03,     
    'V_th_gaba': -0.4   
}


#6. RUNNING THE NORMAL SIMULATION
print("Normal Simulation Results :")
noise = 0.001
V_E,times,spike_timings_e,V_I,spike_timings_i = call_method(normal_parameters,noise)
plot_spikes(V_E,spike_timings_e,V_I,spike_timings_i)

#7. RUNNING THE SCHIZOPHRENIA SIMUATION
print("Schizophrenia Simulation Results :")
noise = 0.005
V_E_Schizo,times_Schizo,spike_timings_e_Schizo,V_I_Schizo,spike_timings_i_Schizo = call_method(Schizophrenia_parameters,noise)
plot_spikes(V_E_Schizo,spike_timings_e_Schizo,V_I_Schizo,spike_timings_i_Schizo)


#8. COMPARISON 
print("Comparison of Normal and Schizophrenia Simulations ")
normal_isi_e = np.diff(spike_timings_e)
normal_isi_i = np.diff(spike_timings_i)
schizo_isi_e = np.diff(spike_timings_e_Schizo)
schizo_isi_i = np.diff(spike_timings_i_Schizo)
bins = np.linspace(7.5, 15, 20)  
plt.figure(figsize=(12,6))
plt.subplot(1,2,1)
n1, bin1, p1 = plt.hist(normal_isi_e,bins =bins,label = "Normal Excitatory Neuron", alpha = 0.5, color = 'deeppink')
n2, bin2, p2 = plt.hist(schizo_isi_e,bins =bins,label = "Schizophrenia Excitatory Neuron", alpha = 0.5,color = 'darkcyan')
ax = plt.gca()
bg = np.array(ax.get_facecolor()[:3])
c1_rgba = p1[0].get_facecolor()
c2_rgba = p2[0].get_facecolor()
c1, a1 = np.array(c1_rgba[:3]), c1_rgba[3]
c2, a2 = np.array(c2_rgba[:3]), c2_rgba[3]
blend = c2*a2 + (c1*a1 + bg*(1-a1))*(1-a2)
plt.fill_between([], [], [], color=blend, label='Overlap')
plt.xlabel("Inter-Spike Interval (ms)")
plt.ylabel("Frequency")
plt.title("Excitatory Neuron ISI Distribution")
plt.legend()
plt.subplot(1,2,2)
plt.hist(normal_isi_i, bins=bins, label="Normal I", color='deeppink', alpha=0.5)
plt.hist(schizo_isi_i, bins=bins, label="Schizo I", color='darkcyan', alpha=0.5)
plt.fill_between([], [], [], color=blend, label='Overlap')
plt.xlabel("Inter-Spike Interval (ms)")
plt.ylabel("Frequency")
plt.title("Inhibitory Neuron ISI Distribution")
plt.legend()
plt.show()

#Calculations
normal_cv_e = np.std(normal_isi_e)/np.mean(normal_isi_e)
normal_cv_i = np.std(normal_isi_i)/np.mean(normal_isi_i)
schizo_cv_e = np.std(schizo_isi_e)/np.mean(schizo_isi_e)
schizo_cv_i = np.std(schizo_isi_i)/np.mean(schizo_isi_i)
normal_fr_e = len(spike_timings_e)/(times[-1]-times[0])
normal_fr_i = len(spike_timings_i)/(times[-1]-times[0])
schizo_fr_e = len(spike_timings_e_Schizo)/(times[-1]-times[0])
schizo_fr_i = len(spike_timings_i_Schizo)/(times[-1]-times[0])
normal_ei_ratio = normal_fr_e/normal_fr_i
schizo_ei_ratio = schizo_fr_e/schizo_fr_i

#CV Plotting
cv_values = [normal_cv_e, normal_cv_i, schizo_cv_e, schizo_cv_i]
labels = ["E Normal", "I Normal", "E Schizo", "I Schizo"]
plt.figure(figsize=(6,4))
plt.bar(labels, cv_values, color=['deeppink','deeppink','darkcyan','darkcyan'], alpha=0.7)
plt.ylabel("CV(ISI)")
plt.title("Comparison of CV(ISI)")
plt.show()

#Firing Rate Plotting
fr_values = [normal_fr_e, normal_fr_i, schizo_fr_e, schizo_fr_i]
labels = ["E Normal", "I Normal", "E Schizo", "I Schizo"]
plt.figure(figsize=(6,4))
plt.bar(labels, fr_values, color=['deeppink','deeppink','darkcyan','darkcyan'], alpha=0.7)
plt.ylabel("Firing Rate (Hz)")
plt.title("Comparison of Firing Rates")
plt.show()

#E/I Ratio Plotting
ei_values = [normal_ei_ratio, schizo_ei_ratio]
labels = ["Normal", "Schizo"]
plt.figure(figsize=(4,4))
plt.bar(labels, ei_values, color=['deeppink','darkcyan'], alpha=0.7)
plt.ylabel("E/I Ratio")
plt.title("Excitatory/Inhibitory Ratio")
plt.show()
