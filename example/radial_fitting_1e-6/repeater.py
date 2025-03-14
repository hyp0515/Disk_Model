import time
import os

for i in range(100000):
    os.system('python plot.py')
    time.sleep(100)
    print(f"{i}")