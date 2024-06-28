# profile-HMM
HMM-based tools for the analysis of biological sequences.

Herramientas basadas en modelos ocultos de Markov para el análisis de secuencias biológicas.

## El módulo "profileHMM.py" contiene cuatro funciones:
  - "HMM": que genera un modelo probabilístico tomando un alineamiento de múltiples secuencias y declarando el alfabeto (ADN o proteína).
  - "forward": una implementación de alforitmo forward para perfiles que evalúa una secuencia respecto a un HMM.
  - "Viterbi_log": una implementación del algoritmo de Viterbi (en su versión logarítmica) que estima los estados subyacentes más probables para una secuencia biológica, dado un HMM.
  - "realign": que utiliza una secuencia biológica y sus estados ocultos asociados para integrarla en un alineamiento de múltiples secuencias.

## El programa "dentista.py" consiste en una práctica que utiliza las herramientas de "profileHMM.py" para realizar un estudio sobre la transmisión del VIH en una colsulta clínica.

Las secuencias genéticas utilizadas (y que se encuentran en la carpeta "Florida dentist") fueron extraídas de la base de datos GenBank del instituto nacional de salud de los Estados Unidos de América, y pueden encontrarse usando los códigos de acceso que van desde ![M90907](https://moshi4.github.io/pyMSAviz/](https://www.ncbi.nlm.nih.gov/nuccore/M90847)) hasta ![M90966](https://moshi4.github.io/pyMSAviz/](https://www.ncbi.nlm.nih.gov/nuccore/M90966).
