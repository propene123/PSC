import numpy as np
from pandas import read_csv
from matplotlib import pyplot
from scipy.optimize import curve_fit


def objective(x, m, c):
    return m*x+c;



dataframe = read_csv('./step1conv', header=None)
data = dataframe.values
x, y = data[:,0], data[:,-1]
x_ar = list(x)
y_ar = list(y)
x_ar = np.array([x_ar])
y_ar = np.array([y_ar])
y_ar = 0.01-y_ar
x_ar = 10.0/x_ar
y_ar = np.log(y_ar)
x_ar = np.log(x_ar)
foo, _ = curve_fit(objective, list(x_ar.flatten()), list(y_ar.flatten()))
m,c = foo
print(m)
print(c)

pyplot.plot(x_ar.flatten(), y_ar.flatten())
pyplot.show()



