from Bio import AlignIO # Lectura de alineamientos.
from Bio import SeqIO # Lectura de secuencias.
# Documentación de biopython: https://biopython.org/wiki/Documentation
from pymsaviz import MsaViz # Visualización del alineamiento.
# Documentación de pymsaviz: https://github.com/moshi4/pyMSAviz
from os import walk
from profileHMM import HMM
from profileHMM import forward
from profileHMM import Viterbi_log
from profileHMM import realign
# Donde 'profileHMM.py' es un script de python que contiene funciones para estudiar secuencias biológicas con HMMs.


########################################################################################################################
### Evaluación de la similitud biológica entre una secuencia y un modelo.
########################################################################################################################

# Lectura del alineamiento múltiple de secuencias a partir del cual se generará el modelo:
MSA = AlignIO.read("Florida dentist/alineaciones/Alineamiento_W.fas", 'fasta')

# Representación del alineamiento: (No recomendado para secuencias largas.)
# for record in MSA:
#     print("%s\t%s" % (record.description, record.seq))

# Visualización del alineamiento: (Se genera un archivo .png.)
mv = MsaViz(MSA, wrap_length=60, show_count=True, show_consensus=True)
mv.savefig("contagio.png")


# Obtención del profile-HMM:
# La función HMM toma como entrada un alineamiento múltiple de secuencias y el alfabeto ('DNA' o 'protein') en el que se
# expresan las secuencias, y devuelve una tupla formada por la secuencia de consenso, la matriz transición de estados y
# la matriz de emisiones de símbolos.
print('Generating profile...')
profile = HMM(MSA.alignment, 'DNA')


# Lista de secuencias estrechamente relacionadas (en base a las cuales se supone contagio).
related_seq = ["Dentista", "Paciente A", "Paciente B", "Paciente C", "Paciente E", "Paciente G"]

# Se calcula el rango de probabilidades para secuencias estrechamente relacionadas.
print('Obtaining range of probabilities...')
prob_range = [1, 0]
for file in related_seq:
    filenames = next(walk("Florida dentist/genbank/" + file), (None, None, []))[2]
    for seq in filenames:
        new = SeqIO.read("Florida dentist/genbank/" + file + '/' + seq, 'fasta')
        # Cálculo de la probabilidad de emisión de una secuencia de acuerdo a un modelo:
        # La función forward toma como entrada una secuencia, un modelo y el alfabeto ('DNA' o 'protein') en el que se
        # expresa la secuencia, y devuelve la probabilidad de que el modelo emita la secuencia de observables introducida.
        evaluation = forward(new.seq, profile, 'DNA')
        if evaluation < prob_range[0]:
            prob_range[0] = evaluation
        
        if evaluation > prob_range[1]:
            prob_range[1] = evaluation

print(prob_range)   # [5.6246801163918036e-42, 4.0763771219361905e-28]


# Lista de secuencias no tan relacionadas. (Casos para los que se descarta el contagio).
unrelated_seq = ["Paciente D", "Paciente F", "Control 1", "Control 2", "Control 3", "Control 4", "Control 5", "Control 6", "Control 7"]

# Se calcula el rango de probabilidades para secuencias poco relacionadas.
print('Obtaining range of probabilities...')
prob_range = [1, 0]
for file in unrelated_seq:
    filenames = next(walk("Florida dentist/genbank/" + file), (None, None, []))[2]
    for seq in filenames:
        new = SeqIO.read("Florida dentist/genbank/" + file + '/' + seq, 'fasta')
        # Cálculo de la probabilidad de emisión de una secuencia de acuerdo a un modelo:
        # La función forward toma como entrada una secuencia, un modelo y el alfabeto ('DNA' o 'protein') en el que se
        # expresa la secuencia, y devuelve la probabilidad de que el modelo emita la secuencia de observables introducida.
        evaluation = forward(new.seq, profile, 'DNA')
        if evaluation < prob_range[0]:
            prob_range[0] = evaluation
        
        if evaluation > prob_range[1]:
            prob_range[1] = evaluation
 
print(prob_range)   # [4.854626613337559e-90, 2.0538712506943715e-53]


# Lectura de la secuencia de entrada:
new = SeqIO.read("Florida dentist/genbank/Paciente H/M90907.fasta", 'fasta')
# Representación de la secuencia de entrada:
#print("%s\t%s" % (new.description, new.seq))

# Evaluación de la secuencia de ADN viral del Paciente H para el modelo.
evaluation = forward(new.seq, profile, 'DNA')
print(evaluation)   # 9.435369511472031e-85


########################################################################################################################
### Agregar una secuencia a un alineamiento múltiple.
########################################################################################################################

# Estimación de la secuencia de estados ocutos que, de acuerdo al modelo, mejor explica la secuencia observable.
estimation = Viterbi_log(new.seq, profile, 'DNA')

# Alineamiento de la secuencia.
realign(new.seq, MSA, profile, estimation, 'Paciente_H')

# Visualización del alineamiento: (Se genera otro archivo .png.)
MSA = AlignIO.read("Paciente_H.fas", 'fasta')
mv = MsaViz(MSA, wrap_length=60, show_count=True, show_consensus=False)
mv.savefig("Paciente_H.png")
