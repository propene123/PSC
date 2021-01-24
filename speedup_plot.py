import numpy as np
from matplotlib import pyplot
from scipy.optimize import curve_fit


x = np.array([[1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24]])
y = np.array([[1868.15, 1039.689, 679.609, 548.631, 467.033, 442.838,
               418.656, 414.006, 398.188, 381.936, 370.069, 374.413, 375.387]])
y = 1868.15/y

def objective(s, p):
    return 1/(p+((1-p)/s))


args, _ = curve_fit(objective, list(x.flatten()), list(y.flatten()))
p = args

pyplot.scatter(x,y, label="Speedup")
new_x = np.arange(1,28)
new_y = objective(new_x, p)
pyplot.xticks([1]+list(range(2,30,2)))
pyplot.plot(new_x, new_y, '--', color='red', label=f"Fitted Amdahl's curve\nf={p}")
pyplot.grid(b=True, which='major', color='#666666')
pyplot.minorticks_on()
pyplot.grid(b=True, which='minor', color='#999999', alpha=0.2)
pyplot.legend(loc="upper left")
pyplot.xlabel('Number of Cores')
pyplot.ylabel('Speedup')
pyplot.show()


