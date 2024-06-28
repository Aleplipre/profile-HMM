# Este script contiene funciones para el estudio de secuencias biológicas con técnicas basadas en modelos ocultos de Markov.
from Bio.motifs import Motif # Consenso.
# Documentación de biopython: https://biopython.org/wiki/Documentation
import numpy as np
import math


# Función para la obtención de un profile-HMM basado en las observaciones.
def HMM(alignment, alphabet):
    # Selección de alfabeto:
    if alphabet == 'DNA' or 'ADN':
        alphabet = 'ACGT'
    
    else:
        alphabet = 'ACDEFGHIKLMNPQRSTVWY'
    
    # Secuencia de consenso:
    # Por cada columna de alineamiento, devuelve un símbolo si aparece en la mitad o más de las celdas.
    # Si ningún símbolo supera el umbral del 50%, entonces devuelve un símbolo '-' si hay más huecos que símbolos y, en
    # caso contrario, devuelve no incluído en el alfabeto 'N' o 'X'.
    motif = Motif(alphabet + '-', alignment)  # Se incluye '-' en el alfabeto.
    consensus = motif.counts.calculate_consensus(identity=0.5)
    
    # Datos sobre el alineamiento:
    alphabet_length = len(alphabet) # Número de símbolos del alfabeto (4 para ADN y 20 para proteínas).
    sequence_number = len(motif.alignment) # Número de secuencias alineadas o filas de la matriz de posiciones específicas.
    alignment_length = motif.alignment.length # Longitud del alineamiento o columnas de la matriz de posiciones específicas.
    insertion_index = [column for column, gap in enumerate(consensus) if gap == '-'] # Columnas de inserción.
    match_columns_number = alignment_length-len(insertion_index) # Número total de columnas de coincidencia o delección.
    

    # Cálculo de la matriz de emisión:
    emissions_matrix = np.ones((2*match_columns_number+1, alphabet_length)) # {I0 M1 I1 ...} x {A C G T}
    inserts = np.zeros(alphabet_length)
    t = 0
    for column in range(alignment_length):
        if column in insertion_index: # Conteo de símbolos en columnas de inserción.
            for symbol in range(alphabet_length):
                inserts[symbol] += motif.counts[symbol, column]
        
        else: # Conteo de símbolos en columnas de coincidencia.
            emissions_matrix[t] += inserts
            inserts = np.zeros(alphabet_length)
            t += 2
            for symbol in range(alphabet_length):
                emissions_matrix[t-1][symbol] += motif.counts[symbol, column]
    
    emissions_matrix[t] += inserts
    # Normalización de cada fila:
    for row in range(2*match_columns_number+1):
        emissions_matrix[row, :] = emissions_matrix[row, :]/sum(emissions_matrix[row, :])
    
    
    # Matriz auxiliar de observaciones:
    # Donde los unos representan si hay un símbolo en cierta posición del alineamiento y los ceros representan si hay un
    # hueco. Para facilitar los cálculos, se incluye una última columna de unos representando el estado final.
    aux_matrix = np.ones((sequence_number, alignment_length+1), dtype=int)
    for sequence in range(sequence_number):
        for symbol in range(alignment_length):
            if alignment[sequence][symbol] == '-':
                aux_matrix[sequence][symbol] = 0
    
    
    # Cálculo de la matriz de transición reducida:
    # Contiene la misma información que la matriz de transición ocupando menos espacio.
    transitions_matrix = np.ones((3*(match_columns_number+1), 3))  # {M0 I0 D1 M1 I1 ...} x {I0 D1 M1 I1 ...}
    t = 0
    insertion_vector = np.zeros(sequence_number) # Cuenta las inserciones acumuladas entre dos columnas de coincidencia (o delección).
    previous_match_column = alignment_length
    for column in range(alignment_length+1):
        if column in insertion_index: # I -> I
            insertion_vector += aux_matrix[:, column]
        
        else: # ? -> ?
            insertion_row = (insertion_vector == False)
            transitions_matrix[t+1, 2] += sum(aux_matrix[:, previous_match_column] * insertion_row * aux_matrix[:, column]) # M -> ... -> M
            transitions_matrix[t+1, 1] += sum((aux_matrix[:, previous_match_column] * insertion_row) > aux_matrix[:, column]) # M -> ... -> D
            transitions_matrix[t, 2] += sum(aux_matrix[:, previous_match_column] < (insertion_row * aux_matrix[:, column])) # D -> ... -> M
            transitions_matrix[t, 1] += sum(np.logical_not(aux_matrix[:, previous_match_column]) * insertion_row * np.logical_not(aux_matrix[:, column])) # D -> ... -> D
            if column-1 in insertion_index: # I -> ?
                transitions_matrix[t+1, 0] += sum(np.logical_not(insertion_row) * aux_matrix[:, previous_match_column]) # I -> M
                transitions_matrix[t, 0] += sum(np.logical_not(insertion_row) > aux_matrix[:, previous_match_column]) # I -> D
                transitions_matrix[t+2, :] += [sum(insertion_vector - np.logical_not(insertion_row)),
                                               sum(np.logical_not(insertion_row) > aux_matrix[:, column]),
                                               sum(np.logical_not(insertion_row) * aux_matrix[:, column])] # I -> [I, D, M]
            
            t += 3
            insertion_vector = np.zeros(sequence_number) # Reinicio del conteo de inserciones.
            previous_match_column = column # Ubicación de la última columna de coincidencia (o delección).
    
    transitions_matrix = np.delete(transitions_matrix, 0, axis=0)
    transitions_matrix[-3:, 1] = np.zeros(3)
    # Normalización de cada fila:
    for row in range(3*match_columns_number+2):
        transitions_matrix[row, :] = transitions_matrix[row, :]/sum(transitions_matrix[row, :])
    
    
    return (consensus, transitions_matrix, emissions_matrix)


# Algoritmo forward para la evaluación de secuencias.
def forward(sequence, model, alphabet):
    # Selección de alfabeto:
    if alphabet == 'DNA' or 'ADN':
        alphabet = 'ACGT'
    
    else:
        alphabet = 'ACDEFGHIKLMNPQRSTVWY'
    
    
    # Conversión de la secuencia de entrada en una lista de números. Indicará la columna de la matriz de emisiones.
    num_seq = []
    for symbol in sequence:
        num_seq.append(alphabet.find(symbol))
    
    
    # Iniciación:
    # Primero se calcula el valor de alpha_0 como la probabilidad de que el sistema transite al primer estado de
    # delección nada más iniciarse, es decir, antes de emitir símbolo alguno.
    alpha_0 = model[1][0, 1] # alpha_0(D1)
    
    # A continuación se calcula la probabilidad de que el sistema se encuentre en los estados I0, D1, M1 y I1
    # respectivamente, habiendo emitido únicamente el primer símbolo de la secuencia de entrada en todos los casos.
    alpha_I0 = [model[1][0, 0] * model[2][0, num_seq[0]]]               # alpha_1(I0)
    alpha_D1 = [alpha_I0[0] * model[1][1, 1]]                           # alpha_1(D1)
    alpha_M1 = [model[1][0, 2] * model[2][1, num_seq[0]]]               # alpha_1(M1)
    alpha_I1 = [alpha_0 * model[1][2, 0] * model[2][2, num_seq[0]]]     # alpha_1(I1)
    
    # Ahora se calculan los componentes de una matriz (old_alpha), que se irá actualizando a lo largo de las iteraciones
    # del paso de recursión.
    for symbol in range(1, len(sequence)):
        alpha_I0.append(alpha_I0[symbol-1] * model[1][1, 0] * model[2][0, num_seq[symbol]]) # alpha(I0)
        alpha_D1.append(alpha_I0[symbol] * model[1][1, 1])                                      # alpha(D1)
        alpha_M1.append(alpha_I0[symbol-1] * model[1][1, 2] * model[2][1, num_seq[symbol]])     # alpha(M1)
        alpha_I1.append(sum([alpha_D1[symbol-1], alpha_M1[symbol-1], alpha_I1[symbol-1]] *
                             model[1][2:5, 0]) * model[2][2, num_seq[symbol]])                  # alpha(I1)
    
    old_alpha = np.array([alpha_D1, alpha_M1, alpha_I1])
    # Matriz con las probabilidades alpha_1 del sistema. Las filas representan la probabilidad de que el sistema se
    # encuentre en un estado D1, M1 o I1, y la columna t indica la probabilidad de que el sistema haya emitido t
    # símbolos de la secuencia de entrada.
    
    
    # Recursión:
    # Primero se genera una nueva matriz alpha, que será equivalente a la matriz old_alpha para el caso tau+1.
    alpha = np.ones([3, len(sequence)])
    
    # Bucle externo. Donde tau indica que el sistema ha llegado a tau estados de coincidencia (o delección).
    for tau in range(1, (len(model[2])-1)//2):
        alpha[0, 0] = sum(old_alpha[:, 0] * model[1][3*tau-1:3*tau+2, 1])               # alpha_1(D)
        alpha[1, 0] = alpha_0 * model[1][3*tau-1, 2] * model[2][2*tau+1, num_seq[0]]    # alpha_1(M)
        alpha_0 = alpha_0 * model[1][3*tau-1, 1]                                        # alpha_0(D)
        alpha[2, 0] = alpha_0 * model[1][3*tau+2, 0] * model[2][2*tau+2, num_seq[0]]    # alpha_1(I)
        
        # Bucle interno. El cual se recorre para cada símbolo de la secuencia de entrada, resultando en el cálculo de la
        # nueva matriz alpha, que contendrá la información actualizada de las probabilidades del sistema.
        for symbol in range(1, len(sequence)):
            alpha[:, symbol] = [sum(old_alpha[:, symbol] * model[1][3*tau-1:3*tau+2, 1]),
                                sum(old_alpha[:, symbol-1] * model[1][3*tau-1:3*tau+2, 2]) * model[2][2*tau+1, num_seq[symbol]],
                                sum(alpha[:, symbol-1] * model[1][3*(tau+1)-1:3*(tau+1)+2, 0]) * model[2][2*tau+2, num_seq[symbol]]]
        
        # Se reinicia alpha, para poder calcular las probabilidades para el siguiente valor de tau.
        old_alpha = alpha
        alpha = np.ones([3, len(sequence)])
    

    # Finalización:
    # Se suman las probabilidades de que el sistema acabe en cada uno de los 3 estados posibles (alpha_T(D_tau),
    # alpha_T(M_tau) y alpha_T(I_tau)), habiendo emitido la secuencia completa O_T, multiplicando cada termino por la
    # probabilidad correspondiente a transicionat de ese estado al estado terminal E.
    return sum(old_alpha[:, -1] * model[1][-3:, 2])


# Algoritmo de Viterbi para la estimación de la secuencias de estados.
def Viterbi_log(sequence, model, alphabet):
    # Selección de alfabeto:
    if alphabet == 'DNA' or 'ADN':
        alphabet = 'ACGT'

    else:
        alphabet = 'ACDEFGHIKLMNPQRSTVWY'
    
    
    # Se expresan como logaritmos las probabilidades contenidas en las matrices de transición y emisión.
    model[1][-3:, 1] = [1, 1, 1]
    A = np.log(model[1])    # Matriz de transición.
    B = np.log(model[2])    # Matriz de emisión.


    # Conversión de la secuencia de entrada en una lista de números.
    num_seq = []
    for symbol in sequence:
        num_seq.append(alphabet.find(symbol))


    # Matriz de rutas de estados. Donde cada celda, asociada a un estado del sistema y un elemento de la secuencia de
    # entrada, contiene una tupla con el nombre del estado del que es más probable que preceda y los índices de dicho
    # estado en esta misma matriz. Inicia vacía y se irá completando de acuerdo a las probabilidades ('delta') de Viterbi.
    path = np.empty((len(A)-1, len(num_seq)+1), dtype='object')


    # Iniciación:
    # El estado D1 para t=0 y los estados I0, D1, M1 y I1 para t=1 son sólo accesibles desde un único estado, por tanto,
    # en estos casos, el valor de la probabilidad de Viterbi ('delta') será idéntico a los valores 'alpha'
    # correspondientes a cada caso, calculados en el algoritmo forward.
    delta_0_D = [A[0, 1]]                         # alpha_0(D1) = delta_0(D1)
    path[1, 0] = 'B'                              # Sólo accesible desde el inicio.
    delta_1_I0 = A[0, 0] + B[0, num_seq[0]]                     # delta_1(I0)
    path[0, 1] = 'B'
    delta_1_D1 = delta_1_I0 + A[1, 1]                           # delta_1(D1)
    path[1, 1] = ('0I', 0, 1)
    delta_1_M1 = model[1][0, 2] + B[1, num_seq[0]]              # delta_1(M1)
    path[2, 1] = 'B'
    delta_1_I1 = delta_0_D[0] + A[2, 0] + B[2, num_seq[0]]      # delta_1(I1)
    path[3, 1] = ('1D', 1, 0)

    # A continuación se crean dos vectores útiles para la iteración:
    # El vector delta_0_D contiene las probabilidades de acceder a todos los estados de delección que componen el modelo
    # sin emitir ningún símbolo de la secuencia. Este vector será útil para calcular delta.
    for column in range(1, (len(B)-1)//2):
        delta_0_D.append(delta_0_D[column-1] + A[3*column-1, 1])    # delta_0(D)
        path[3*column+1, 0] = (str(column)+'D', 3*column-2, 0)

    # El vector delta contiene las probabilidades de acceder a un determinado estado del sistema habiendo emitido el
    # primer símbolo de la secuencia de entrada.
    delta = [delta_1_I0, delta_1_D1, delta_1_M1, delta_1_I1]        # delta_1(X)

    # El diccionario de estados 'state' servirá para interpretar de qué tipo de estado procede el sistema de acuerdo a
    # cuál sea el mayor de los argumentos en el cálculo de la probabilidad de Viterbi.
    state = {0:'D', 1:'M', 2:'I'}

    # Bucle para calcular delta completo, es decir, el vector con todas las probabilidades de cada estado en t=1.
    for column in range(1, (len(B)-1)//2):
        delta.append(max(delta[3*column-2:3*column+1] + A[3*column-1:3*column+2, 1]))       # delta_1(D)
        # Para evitar calculos, se intruduce el escalar 'previous', que indica cuál de los 3 estados previos a aquel en
        # el que se encuentra el sistema es más probable que le preceda.
        previous = np.argmax(delta[3*column-2:3*column+1] + A[3*column-1:3*column+2, 1])
        path[3*column+1, 1] = (str(column)+state[previous], 3*column+previous-2, 1)

        delta.append(delta_0_D[column-1] + A[3*column-1, 2] + B[2*column+1, num_seq[0]])    # delta_1(M)
        path[3*column+2, 1] = (str(column)+'D', 3*column-2, 0)
        delta.append(delta_0_D[column] + A[3*column+2, 0] + B[2*column+2, num_seq[0]])      # delta_1(I)
        path[3*column+3, 1] = (str(column+1)+'D', 3*column+1, 0)


    # Recursión:
    # Bucle externo. Donde 'symbol' indica el número de símbolos emitidos por el sistema.
    for symbol in range(1, len(num_seq)):
        
        # Bucle interno.
        # Primero, se guarda el vector delta actual en old_delta y se calcula el siguiente vector alpha para t+1.
        old_delta = delta

        # Luego se calculan los 4 primeros miembros del vector delta, que contendrá las probabilidades de Viterbi para los
        # estados I0, D1, M1 y I1 al tiempo de la emisión del t-ésimo observable.
        delta = [old_delta[0] + A[1, 0] + B[0, num_seq[symbol]]]                # delta_t(I0)
        path[0, symbol+1] = ('0I', 0, symbol)
        delta.append(delta[0] + A[1, 1])                                        # delta_t(D1)
        path[1, symbol+1] = ('0I', 0, symbol+1)
        delta.append(old_delta[0] + A[1, 2] + B[1, num_seq[symbol]])            # delta_t(M1)
        path[2, symbol+1] = ('0I', 0, symbol)
        delta.append(max(old_delta[1:4] + A[2:5, 0]) + B[2, num_seq[symbol]])   # delta_t(I1)
        previous = np.argmax(old_delta[1:4] + A[2:5, 0])
        path[3, symbol+1] = ('1'+state[previous], 1, symbol)

        # Cada iteración del bucle calcula los 3 siguientes elentos del vector delta, es decir, las probabilidades para los
        # estados D, M e I posteriores, dado el mismo valor de t (symbol).
        for column in range(1, (len(B)-1)//2):
            delta.append(max(delta[3*column-2:3*column+1] + A[3*column-1:3*column+2, 1]))               # delta_t(D)
            previous = np.argmax(delta[3*column-2:3*column+1] + A[3*column-1:3*column+2, 1])
            path[3*column+1, symbol+1] = (str(column)+state[previous], 3*column+previous-2, symbol+1)

            delta.append(max(old_delta[3*column-2:3*column+1] + A[3*column-1:3*column+2, 2]) +
                        B[2*column+1, num_seq[symbol]])                                                 # delta_t(M)
            previous = np.argmax(old_delta[3*column-2:3*column+1] + A[3*column-1:3*column+2, 2])
            path[3*column+2, symbol+1] = (str(column)+state[previous], 3*column+previous-2, symbol)

            delta.append(max(old_delta[3*column+1:3*column+4] + A[3*column+2:3*column+5, 0]) +
                        B[2*column+2, num_seq[symbol]])                                                 # delta_t(I)
            previous = np.argmax(old_delta[3*column+1:3*column+4] + A[3*column+2:3*column+5, 0])
            path[3*column+3, symbol+1] = (str(column+1)+state[previous], 3*column+previous+1, symbol)


    # Finalización:
    # Búsqueda de la ruta más probable:
    previous = np.argmax(delta[-3:] + A[-3:, 2])
    last_path = (str((len(B)-1)//2)+state[previous], 3*(len(B)-1)//2+previous-2, len(num_seq))
    best_path = ''  # Cadena que registra los estados transitados.
    while last_path != 'B': # Para cuando un estado remita al inicio 'B'. (Se recorre en sentido inverso.)
        best_path += last_path[0]+' '
        last_path = path[last_path[1],last_path[2]] # Estado anterior al actual.

    # Se invierte la cadena para obtener una ruta legible.
    best_path = best_path[::-1]
    best_path = best_path[1:]


    # Por último, se calcula la probabilidad de Viterbi para la secuencia de estados más probable, es decir, la
    # probabilidad de que el sistema transite por una serie de estados determinados, de forma tal que emita la secuencia
    # de entrada.
    return (math.exp(max(delta[-3:] + A[-3:, 2])), best_path)


# Función para la alineación de una secuencia con un alineamiento de múltiples secuencias.
def realign(sequence, alignment, model, estimation, name):
    # El primer bucle alinea la secuencia de entrada con la secuencia de consenso del alineamiento múltiple, de acuerdo
    # a la secuencia de estados predicha por el algoritmo de Viterbi.
    realignment = [[], []]  # Secuencia de consenso y nueva secuencia.
    consensus = list(model[0].replace('-', '')) # Se parte de la secuencia de consenso sin huecos.
    sequence = list(sequence)
    for state in estimation[1].split():
        # En un estado de delección hay un símbolo en el consenso que no aparece en la secuencia de entrada.
        if state[0] == 'D':
            realignment[0].append(consensus.pop(0))
            realignment[1].append('-')
        
        # En un estado de coincidencia aparecce un símbolo en el consenso y otro en la secuencia de entrada.
        if state[0] == 'M':
            realignment[0].append(consensus.pop(0))
            realignment[1].append(sequence.pop(0))
        
        # En un estado de inserción aparece un símbolo en la secuencia de entrada que no se incluye en el consenso.
        if state[0] == 'I':
            realignment[0].append('-')
            realignment[1].append(sequence.pop(0))
    
    
    # El segundo bucle genera espacios adicionales en la secuencia de entrada o columnas de huecos en el alineamiento
    # múltiple, para que el nuevo alineamiento encaje bien.
    consensus = list(model[0])  # Se parte de la secuencia de consenso con huecos.
    new_alignment = np.empty((len(alignment), 0)).tolist()  # Lista que contendrá el nuevo alineamiento.
    for seq in range(len(alignment)):
        new_alignment[seq] = list(alignment.alignment[seq])
    

    for symbol in range(len(model[0])):
        # Si no aparece un hueco en la columna de inserción, se incluye a la secuencia.
        if consensus[symbol] == '-' and realignment[0][symbol] != '-':
            realignment[0].insert(symbol, '-')
            realignment[1].insert(symbol, '-')
        
        # Si falta una columna de inserción pero hay un hueco en la secuencia, se inserta una columna vacía.
        if consensus[symbol] != '-' and realignment[0][symbol] == '-':
            consensus.insert(symbol, '-')
            for seq in range(len(alignment)):
                new_alignment[seq].insert(symbol, '-')
    
    # Como se ha recorrido la secuencia de símbolos contenidos en el consenso, puede ocurrir que la nueva secuencia sea
    # mayor. En este caso, los símbolos de la secuencia que quedan se interpretan como inserciones, y se generan nuevas
    # columnas de huecos para el nuevo alineamiento.
    
    while len(realignment[1]) < len(new_alignment[0]):
        realignment[1].append('-')
    
    # Por último se genera un archivo .fas que contendrá el nuevo alineamiento de secuencias.
    MSA = open(name + '.fas', 'w')
    for seq in range(len(alignment)):
        print('>'+alignment[seq].description, file = MSA)
        print(''.join(new_alignment[seq]), file = MSA)
    
    print('>'+name, file=MSA)
    print(''.join(realignment[1]), file = MSA)
    MSA.close()
