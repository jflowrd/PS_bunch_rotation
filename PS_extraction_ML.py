import matplotlib.pyplot as plt
import numpy as np
import tensorflow as tf

from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, mean_squared_error
from tensorflow.keras import layers, models
from scipy.optimize import curve_fit


# Load your data
data0 = np.load('results/HTcondor_test/bunch_lengths_0.npz')
data1 = np.load('results/HTcondor_test/bunch_lengths_1.npz')
data2 = np.load('results/HTcondor_test/bunch_lengths_2.npz')
data3 = np.load('results/HTcondor_test/bunch_lengths_3.npz')
data4 = np.load('results/HTcondor_test/bunch_lengths_4.npz')

bunch_len_string = 'fwhm_save' # 'gaussian_4sig_save', 'fwhm_save'

time = (data0['cycle_time'] *1e6, data1['cycle_time'] *1e6, data2['cycle_time'] *1e6, data3['cycle_time'] *1e6, data4['cycle_time'] *1e6)
length = (data0[bunch_len_string].reshape(-1) * 1e9, data1[bunch_len_string].reshape(-1) * 1e9, data2[bunch_len_string].reshape(-1) * 1e9, data3[bunch_len_string].reshape(-1) * 1e9, data4[bunch_len_string].reshape(-1) * 1e9)

extraction_time = np.concatenate(time) 
bunch_length = np.concatenate(length)  

x_data = (extraction_time) / extraction_time.max()
y_data = (bunch_length) - bunch_length.mean()


plt.plot(extraction_time, bunch_length, '.')
plt.xlabel('Extraction time [us]')
plt.ylabel('FWHM [ns]')
plt.show()

# Build a simple MLP model with nonlinear activation functions for nonlinear regression
model = tf.keras.Sequential([
    layers.Dense(128, input_shape=[1], activation='relu'),  # First hidden layer with ReLU activation
    layers.Dense(128, activation='relu'),  # Second hidden layer with ReLU activation
    layers.Dense(128, activation='relu'),  # Second hidden layer with ReLU activation
    layers.Dense(1)  # Output layer
])

# Compile the model specifying the optimizer and loss function
model.compile(optimizer='adam', loss='mean_squared_error')

# Train the model with the data
model.fit(x_data, y_data, epochs=300, verbose=0)

min_extraction_time = np.min(x_data) #np.min(extraction_time)  # Or set your minimum
print("Minimum Extraction Time:", min_extraction_time)
max_extraction_time = np.max(x_data) #np.max(extraction_time)  # Or set your maximum
print("Maximum Extraction Time:", max_extraction_time)
extraction_time_grid = np.linspace(min_extraction_time, max_extraction_time, 10000).reshape(-1, 1)

x_test = extraction_time_grid
y_pred = model.predict(x_test)

optimal_index = np.argmin(y_pred)
optimal_extraction_time = x_test[optimal_index][0] * extraction_time.max()
optimal_bunch_length = y_pred[optimal_index][0] + bunch_length.mean()
opt_bunch_length_str = "Bunch Length: " + str(round(optimal_bunch_length,2)) + " ns"
opt_extraction_time_str = "Ext. time: " + str(round(optimal_extraction_time,2)) + " us"
title = opt_extraction_time_str + " " + opt_bunch_length_str

print("Optimal Extraction Time:", optimal_extraction_time)
print("Optimal Bunch Length:", optimal_bunch_length)

x_data *= extraction_time.max()
x_test *= extraction_time.max()
y_data += bunch_length.mean()
y_pred += bunch_length.mean()

plt.scatter(x_data, y_data, color='blue', label='Data')
plt.plot(x_test, y_pred, color='red', label='Model')
plt.xlabel('Extraction time [us]')
plt.ylabel('Bunch length (FWHM) [ns]')
plt.title(title)
plt.legend()
plt.show()

'''
def cos_func(x, A, B, C, K):
    y = (A*np.cos(K*x) + C)*np.exp(-B * x)
    return y

parameters, covariance = curve_fit(cos_func, x_data, y_data, p0=(8, 0.001, 6, 2e-3),  maxfev = 5000)


fit_cosine = cos_func(x_data, *parameters)

plt.plot(x_data, y_data, 'o', label='data')
plt.plot(x_data, fit_cosine, '-', label='fit')
plt.show()
'''
'''
# Split the dataset
X_train, X_test, y_train, y_test = train_test_split(extraction_time_norm.reshape(-1, 1), bunch_length_norm, test_size=0.5, random_state=42)

plt.plot(X_train, y_train, '.')
plt.plot(X_test, y_test, '.')
plt.xlabel('Extraction time [us]')
plt.ylabel('FWHM [ns]')
plt.show()

# Build the model
model = models.Sequential([
    layers.Dense(128, input_shape=(1,), activation='relu'),
    layers.Dense(128, activation='relu'),
    layers.Dense(1)
])

opt = tf.keras.optimizers.Adam(0.01)
model.compile(optimizer=opt, loss='mse') #, metrics=['mae']
# Train the model
history = model.fit(X_train, y_train, epochs=100) #validation_split=0.2, verbose=0, callbacks=[tf.keras.callbacks.EarlyStopping(patience=10)]

y_pred = model.predict(X_test)

plt.plot(X_test, y_pred, '.', label='Predicted')
plt.show()

mae = mean_absolute_error(y_test, y_pred) # average absolute difference between predicted and actual values
mse = mean_squared_error(y_test, y_pred) # penalizes larger errors more than MAE does
rmse = np.sqrt(mse)

print(f"MAE: {mae}")
print(f"MSE: {mse}")
print(f"RMSE: {rmse}")


min_extraction_time = np.min(extraction_time)  # Or set your minimum
print("Minimum Extraction Time:", min_extraction_time)
max_extraction_time = np.max(extraction_time)  # Or set your maximum
print("Maximum Extraction Time:", max_extraction_time)
extraction_time_grid = np.linspace(min_extraction_time, max_extraction_time, 10000).reshape(-1, 1)

predicted_bunch_lengths = model.predict(extraction_time_grid).flatten()

plt.plot(extraction_time_grid, predicted_bunch_lengths, '.', label='Predicted')
plt.show()

optimal_index = np.argmin(predicted_bunch_lengths)
optimal_extraction_time = extraction_time_grid[optimal_index] #* extraction_time.std() + extraction_time.mean()
optimal_bunch_length = predicted_bunch_lengths[optimal_index] #* bunch_length.std() + bunch_length.mean()

print("Optimal Extraction Time:", optimal_extraction_time)
print("Optimal Bunch Length:", optimal_bunch_length)

# Plot the loss and MAE
plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
plt.plot(history.history['loss'], label='Loss')
#plt.plot(history.history['val_loss'], label='Validation Loss')
plt.title('Loss Over Epochs')
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.legend()

# plt.subplot(1, 2, 2)
# plt.plot(history.history['mae'], label='MAE')
# plt.plot(history.history['val_mae'], label='Validation MAE')
# plt.title('MAE Over Epochs')
# plt.xlabel('Epoch')
# plt.ylabel('MAE')
# plt.legend()

plt.tight_layout()
plt.show()


opt2 = tf.keras.optimizers.Adam(0.01)
model.compile(optimizer=opt, loss='mse')
r = model.fit(extraction_time, bunch_length, epochs=100)

# Plot the loss
plt.plot(r.history['loss'], label='loss')
plt.show()

predicted_bunch_lengths = model.predict(extraction_time_grid).flatten()
plt.plot(extraction_time_grid, predicted_bunch_lengths, '.', label='Predicted')
plt.plot(extraction_time, bunch_length, '.')
plt.xlabel('Extraction time [us]')
plt.ylabel('FWHM [ns]')
plt.show()

'''