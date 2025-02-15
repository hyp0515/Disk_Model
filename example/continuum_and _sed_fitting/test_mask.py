import numpy as np

a = np.array([[1, 2, 3],
             [4, 5, 6],
             [7, 8, 9]])

mask = a < 5
a[mask] = 0    
print(mask)
print(a)