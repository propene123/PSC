import numpy as np
from pandas import read_csv
from matplotlib import pyplot
from scipy.optimize import curve_fit


def objective(x, m, c):
    return m*x+c;



dataframe = read_csv('./step2conv', header=None)
data = dataframe.values
x, y = data[:,0], data[:,-1]
x_ar = list(x)
y_ar = list(y)
for i in range(0, len(y_ar)-1):
    y_ar[i] = abs(y_ar[i] - y_ar[i+1])

y_ar = y_ar[:len(y_ar)-1]
x_ar = x_ar[:len(x_ar)-1]

x_ar = np.array([x_ar])
y_ar = np.array([y_ar])
x_ar = 10.0/x_ar
y_ar = np.log(y_ar)
x_ar = np.log(x_ar)
foo, _ = curve_fit(objective, list(x_ar.flatten()), list(y_ar.flatten()))
m,c = foo
print(m)
print(c)
trend = np.arange(9,16)
trend_y = objective(trend, m, c)
pyplot.scatter(x_ar.flatten(), y_ar.flatten())
pyplot.plot(trend, trend_y, '--', color='red')
pyplot.show()



