import numpy as np
from pandas import read_csv
from matplotlib import pyplot
from scipy.optimize import curve_fit


def objective(h,p,c):
    return c*(h**p)



dataframe = read_csv('./step1conv', header=None)
data = dataframe.values
x, y = data[:,0], data[:,-1]
x_ar = np.array(list([x]))
y_ar = np.array(list([y]))

y_ar = 0.01 - y_ar
foo, _ = curve_fit(objective, list(x_ar.flatten()), list(y_ar.flatten()))
m,c = foo
print(m)
x_ar = 10.0/x_ar
y_ar = np.log(y_ar)
x_ar = np.log(x_ar)




pyplot.scatter(x_ar, y_ar)
pyplot.show()



