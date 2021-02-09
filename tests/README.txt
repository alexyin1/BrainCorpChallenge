Naming Convention:
	mm ~ matrix multiplication
	tp ~ transpose
	'_XXX' ~ test number

MM test file will contain matrices A, B, and C, each seperated by a newline
TP test file will contain matrices A and At, each seperated by a newline

matrix will have following format:
~example.txt~
,datatype,n_rows,n_cols,
,matrix[0][0],,matrix[0][1],...,,matrix[0][n_cols-1],
...
,matrix[n_rows-1][0],...,,matrix[n_rows-1][n_cols-1],

,matrixtype,n_rows,n_cols,
,matrix[0][0],,matrix[0][1],...,,matrix[0][n_cols-1],
...
,matrix[n_rows-1][0],...,,matrix[n_rows-1][n_cols-1],