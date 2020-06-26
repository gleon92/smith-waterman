import pandas as pd
from Bio import SeqIO
from Bio.SubsMat import MatrixInfo


def maximo(a, b, c):
    Max = a
    if a <= 0 and b < 0 and c < 0:
        return 0
    elif b > Max:
        Max = b
    elif c > Max:
        Max = c
        if b > c:
            Max = b
    return Max


def convert_csv(matrix_file):
    """
    Función que lee un archivo csv y retorna un diccionario conteniendo el valor para cada par de proteínas
    por ejemplo: matriz['A', 'A'] = 1
    :param matrix_file:
    :return: matriz de sustitución
    """
    df = pd.read_csv(matrix_file, sep=',', header=None)
    data = df.values
    bases = ['A', 'C', 'G', 'T']
    matriz = {}
    for i, base_i in enumerate(bases):
        for j, base_j in enumerate(bases):
            matriz[base_i, base_j] = data[i][j]
    return matriz


def smith_waterman(seq_1, seq_2, gap_penalty, matriz_sustitucion=None):

    def backtrack(t, r, str1, str2, x, y, s1='', s2=''):
        """
        Función recursiva para generar las cadenas alineadas
        :param t: matriz de alineamiento (score matrix)
        :param r: salida de resultados
        :param str1: primera secuencia
        :param str2: segunda secuencia
        :param x: última posición horizontal de la matriz
        :param y: última posición vertical en la matriz
        :param s1: alineamiento de la secuencia 1
        :param s2: alineamiento de la secuencia 2
        """

        if x > 0 or y > 0:  # condición base
            c = t[x][y]
            if c == 0 and x > 0 and y > 0:
                arriba = (t[x - 1][y])
                izquierda = (t[x][y - 1])
                diagonal = (t[x - 1][y - 1])
                if diagonal == maximo(arriba, izquierda, diagonal):  # si el valor actual vino de la diagonal
                    backtrack(t, r, str1, str2,
                              x - 1, y - 1, str1[y - 1] + s1, str2[x - 1] + s2)
                if izquierda == maximo(arriba, izquierda, diagonal):  # si el valor actual vino de la izquierda
                    backtrack(t, r, str1, str2,
                              x, y - 1, str1[y - 1] + s1, '-' + s2)
                if arriba == maximo(arriba, izquierda, diagonal):  # si el valor actual vino de arriba
                    backtrack(t, r, str1, str2,
                              x - 1, y, '-' + s1, str2[x - 1] + s2)
            else:
                arriba = c == (t[x - 1][y] + gap_penalty)
                izquierda = c == (t[x][y - 1] + gap_penalty)
                diagonal = c == (t[x - 1][y - 1] + (
                    matrix[str1[y - 1], str2[x - 1]] if (str1[y - 1], str2[x - 1]) in matrix else matrix[
                        str2[x - 1], str1[y - 1]]))
                if diagonal:  # si el valor actual vino de la diagonal
                    backtrack(t, r, str1, str2,
                              x - 1, y - 1, str1[y - 1] + s1, str2[x - 1] + s2)
                if izquierda:  # si el valor actual vino de la izquierda
                    backtrack(t, r, str1, str2,
                              x - 1, y, str1[y - 1] + s1, '-' + s2)
                if arriba:  # si el valor actual vino de arriba
                    backtrack(t, r, str1, str2,
                              x, y - 1, '-' + s1, str2[x - 1] + s2)
                if x == 0 or y == 0:
                    if x == 0:
                        backtrack(t,r,str1,str2,0,0, str1[0:y] + s1, '-'*y + s2)
                    else:
                        backtrack(t, r, str1, str2, 0, 0, '-' * x + s1, str2[0:x] + s2)

        else:
            r.append((s1, s2))

    if matriz_sustitucion:
        matrix = convert_csv(matriz_sustitucion)
    else:
        matrix = MatrixInfo.blosum62

    # se llena el score matrix inicialmente con 0
    score_matrix = [[0 for j in range(len(seq_2) + 1)] for i in range(len(seq_1) + 1)]

    # se calcula los valores del resto de la matriz usando el algoritmo
    for i, val_i in enumerate(seq_2):
        for j, val_j in enumerate(seq_1):
            score_matrix[i + 1][j + 1] = maximo(
                score_matrix[i][j] + (matrix[val_i, val_j] if ((val_i, val_j) in matrix) else matrix[val_j, val_i]),
                score_matrix[i][j + 1] + gap_penalty,
                score_matrix[i + 1][j] + gap_penalty)

    results = []

    # arriba = score_matrix[i+1][j]
    # izquierda = score_matrix[i][j + 1]
    # diagonal = score_matrix[i+1][j+1]
    # if diagonal > arriba and diagonal > izquierda:
    #     pass
    # elif arriba > izquierda and arriba > diagonal:
    #     j = j - 1
    # else:
    #     i = i -1
    # se llama a la función backtrack para obtener las cadenas alineadas
    score = score_matrix[i + 1][j + 1]
    backtrack(score_matrix, results, seq_1, seq_2, i + 1, j + 1)

    for i in score_matrix:
        print(i)

    for n, r in enumerate(results):
        print('Resultado N°: ' + str(n + 1))
        print(r[1])
        print(r[0])
        print('Score: ' + str(score))


in_handle = open('flna.fasta')
record_iterator = SeqIO.parse(in_handle, "fasta")
rec_1 = next(record_iterator).upper()
rec_2 = next(record_iterator).upper()
smith_waterman('AAG', 'AGC', -5, 'matriz_sustitucion.csv')
# needleman_wunsch(rec_1, rec_2, -5)
