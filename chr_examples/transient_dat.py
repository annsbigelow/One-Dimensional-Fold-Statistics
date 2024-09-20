import numpy as np

x=np.linspace(0,10,101)

for xx in x:
    print(xx,2+0.1/(xx+1)+np.random.random()*0.01)
