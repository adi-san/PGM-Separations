import numpy as np

def sigmoid(x):
    return 1 / (1 + np.exp(-x))
def sigmoid_derivative(x):
    return sigmoid(x) * (1 - sigmoid(x))
def relu(x):
    return np.maximum(0, x)
def relu_derivative(x):
    return np.where(x > 0, 1, 0)
def mse(y_true, y_pred):
    return np.mean((y_true - y_pred) ** 2)
def mse_derivative(y_true, y_pred):
    print("y_true:")

print(1)

