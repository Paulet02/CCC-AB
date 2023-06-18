from collections import Counter
import math
import numpy as np
import matplotlib.pyplot as plt
import qsimov as qj


# Count the total results of the iterations
def count_results(results, integer=False):
    bool_to_str = lambda x: str(int(x))
    bool_to_str_list = lambda x: ''.join(list(map(bool_to_str, x))[::-1])
    bool_to_int_list = lambda x: int(''.join(list(map(bool_to_str, x))[::-1]), 2)
    map_results = lambda x: list(map(bool_to_str_list, x))
    map_results_int = lambda x: list(map(bool_to_int_list, x))

    if integer:
        return dict(Counter(map_results_int(results)))
    return dict(Counter(map_results(results)))


# Graficar la distribución de probabilidad de variables
def plot_state_prob(state, list_qubits):
    mapping = {}
    probs = {}
    list_qubits.sort()
    n_qubits = int(np.log2(len(state)))
    for i in range(len(list_qubits)):
        mapping[i] = list_qubits[i]
    for i in range(2**n_qubits):
        bin_number = format(i, '0' + str(n_qubits) + 'b')
        bin_number = bin_number[::-1]
        bin_state = ''
        for j in range(len(list_qubits)):
            bin_state = bin_state + bin_number[mapping[j]]
        bin_state = bin_state[::-1]
        if int(bin_state,2) not in probs.keys():
            probs[int(bin_state,2)] = state[i]**2
        else:
            probs[int(bin_state,2)] = probs[int(bin_state,2)] + state[i]**2
    
    x = list(probs.keys())
    x.sort()
    x_ticks = []
    x_ticks_labels = []
    y = []
    for i in range(len(x)):
        x_ticks.append(x[i])
        x_ticks_labels.append(format(x[i], '0' + str(len(list_qubits)) + 'b'))
        y.append(probs[x[i]])
    fig, ax = plt.subplots()
    ax.bar(x,y,align='center')
    ax.set_xlabel('State')
    ax.set_ylabel('Probability')
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_ticks_labels, rotation=65)
    plt.show()
    
    
# Imprimir la información de un estado cuántico
def print_amplitudes(state):
    print("Información del estado")
    print("---------------------------------------")
    n_qubits=int(np.log2(len(state)))
    for j in range(len(state)):
        bin_number=[0 if (len("{0:b}".format(j))+i)-(n_qubits)<0 else int("{0:b}".format(j)[len("{0:b}".format(j))+i-(n_qubits)]) for i in range(n_qubits)]
        print(str(j)+"   "+str(bin_number)+"     "+str(state[j]))
        
    
# Tabla de verdad de una función clásica representada con un operador cuántico
def truth_table_classical_function(qGate, n_qubits, n_qubits_in, n_ancilla=0):#El resultado se vuelca en el menos significativo
    print("Truth table")
    print("------------------------")
    results = {}
    for j in range(2**(n_qubits_in)):
        bin_number=[0 if (len("{0:b}".format(j))+i)-(n_qubits_in)<0 else int("{0:b}".format(j)[len("{0:b}".format(j))+i-(n_qubits_in)]) for i in range(n_qubits_in)]
        bin_number.reverse()
        O=([['X', [], []] if i else None for i in bin_number])
        q_operator = qj.QGate(n_qubits, 0, "")
        q_operator.add_line(*[None for i in range(n_qubits-n_qubits_in)], *O)
        q_operator.add_line(qGate)
        bin_number.reverse()
        registry = qj.QRegistry(n_qubits)
        registry.apply_gate(q_operator)
        results[int(''.join([str(s) for s in bin_number]),2)] = str(bin_number)+"    "+str(registry.measure([1 if i < (n_qubits - n_qubits_in - n_ancilla) else 0 for i in range(n_qubits)], remove=False))
    for i in range(2**(n_qubits_in)):
        print(results[i])
    print("------------------------")
    
    
    
# Gráfica de la probabilidad de medir un estado solución
def plot_prob_sol(n_qubits, k_sol, max_iterations):
    y = []
    x = []
    for i in range(max_iterations):
        x.append(i)
        y.append(prob_sol(n_qubits, k_sol, i))
        
    plt.xticks(ticks=x)
    plt.bar(x, y)
    plt.show()
    
    
    
    
# Cálculo de la probabilidad de medir un estado solución
def alpha(n_qubits, k_sol, iterations):
    if iterations == 0:
        return 1 / math.sqrt(2**n_qubits)
    else:
        return beta(n_qubits, k_sol, iterations - 1) * (2 - (2 * k_sol) / 2**n_qubits) + alpha(n_qubits, k_sol, iterations - 1) * (1 - (2 * k_sol) / 2**n_qubits)

def beta(n_qubits, k_sol, iterations):
    if iterations == 0:
        return 1 / math.sqrt(2**n_qubits)
    else:
        return beta(n_qubits, k_sol, iterations - 1) * (1 - (2 * k_sol) / 2**n_qubits) - alpha(n_qubits, k_sol, iterations - 1) * ((2 * k_sol) / 2**n_qubits)

    
def prob_sol(n_qubits, k_sol, iterations):
    prob = alpha(n_qubits, k_sol, iterations)
    return prob**2