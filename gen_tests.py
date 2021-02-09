import numpy as np
import math
import random

MIN_VAL = -20
MAX_VAL = 20

d = {
    'int' : np.intc,
    'float' : np.single,
    'double' : np.double  
}

def write(mat, rows, cols, dtype, fname):
    with open(fname, 'a+') as f:
        f.write(str(rows) + ',' + str(cols) + ',' + dtype + '\n')
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                f.write(str(mat[i][j]) + ',')
            f.write('\n')
        f.write('\n')

def get_matrix(row, col, mtype, dtype):
    if mtype == 'zeros':
        mat =  np.zeros(row, col)
    elif mtype == 'ones':
        mat = np.ones(row, col)
    elif mtype == 'iden':
        mat = np.identity(row)
    elif mtype == 'tril':
        mat = np.random.uniform(low=MIN_VAL, high=MAX_VAL, size=(row, col))
        mat.tril()
    elif mtype == 'sym':
        mat = np.random.uniform(low=MIN_VAL, high=MAX_VAL, size=(row, col))
        mat = (mat + mat.T)/2
    else:
        mat = np.random.uniform(low=MIN_VAL, high=MAX_VAL, size=(row, col))
    return np.round(mat, 4).astype(dtype)

def tp_test(row, col, dtype='int', mtype='', fname=''):
    a = get_matrix(row, col, mtype, dtype)
    write(a, row, col, dtype, fname)
    write(a.T, col, row, dtype, fname)

def mm_test(a_row, a_col, b_row, b_col, a_type='', a_dtype='int', b_type='', b_dtype='int', fname=''):
    a = get_matrix(a_row, a_col, a_type, d[a_dtype])
    b = get_matrix(b_row, b_col, b_type, d[b_dtype])
    write(a, a_row, a_col, a_dtype, fname)
    write(b, b_row, b_col, b_dtype, fname)
    write(a@b, a_row, b_col, a_dtype, fname)

#mm tests
#Tests 0-5
for i in range(6):
    val = int(math.pow(2, i))
    f_num = f"{i:03}"
    fname = 'mm_' + f_num + '.txt'
    print(fname)
    mm_test(val, val, val, val, fname=fname)

#Tests 6-8
for i, key in enumerate(d.keys()):
    f_num = f"{i+6:03}"
    fname = 'mm_' + f_num + '.txt'
    print(fname, key)
    mm_test(32, 32, 32, 32, a_dtype=key, b_dtype=key, fname=fname)

#Test 9-15 rectangular matrix
for i in range(9, 15):
    a_row = random.randint(1, 30)
    k = random.randint(1, 10)
    b_row = random.randint(1, 30)
    f_num = f"{i:03}"
    fname = 'mm_' + f_num + '.txt'
    mm_test(a_row, k, k, b_row, a_dtype='float', b_dtype='float', fname=fname)

#tp tests
for i in range(6):
    val = int(math.pow(2, i))
    f_num = f"{i:03}"
    fname = 'tp_' + f_num + '.txt'
    print(fname)
    tp_test(val, val, fname=fname)

for i in range(6, 11):
    f_num = f"{i:03}"
    fname = 'tp_' + f_num + '.txt'
    print(fname)
    a = random.randint(0, 100)
    b = random.randint(0, 100)
    tp_test(a, b, fname=fname)

