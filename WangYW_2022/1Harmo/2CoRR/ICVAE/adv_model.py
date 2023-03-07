from tensorflow.keras.layers import Input, Dense
from tensorflow.keras.models import Model
import tensorflow.keras.callbacks
import tensorflow as tf
import os
from matplotlib import pyplot as plt
from IPython.display import clear_output
import tensorflow.keras.backend as K
from utils import onehot_test_label
import warnings
warnings.filterwarnings('ignore')
class LossHistory(tensorflow.keras.callbacks.Callback):
    def on_train_begin(self, logs={}):
        self.losses = []
 
    def on_batch_end(self, batch, logs={}):
        self.losses.append(logs.get('loss'))

class PlotLosses(tensorflow.keras.callbacks.Callback):
    def on_train_begin(self, logs={}):
        self.i = 0
        self.x = []
        self.losses = []
        self.val_losses = []
        
        self.fig = plt.figure()
        
        self.logs = []

    def on_epoch_end(self, epoch, logs={}):
        
        self.logs.append(logs)
        self.x.append(self.i)
        self.losses.append(logs.get('loss'))
        self.val_losses.append(logs.get('val_loss'))
        self.i += 1
        
        clear_output(wait=True)
        plt.plot(self.x, self.losses, label="loss")
        plt.plot(self.x, self.val_losses, label="val_loss")
        plt.legend()
        if epoch%100==0:
            plt.show();
        
plot_losses = PlotLosses()
#(VAEoutput,labels,epoch)
def adv_model():
   # This returns a tensor
    inputs = Input(shape=(512,))
    # a layer instance is callable on a tensor, and returns a tensor
    x = Dense(32, activation='tanh')(inputs)
    x = Dense(32, activation='tanh')(x)
    predictions = Dense(6, activation='softmax')(x)
    # This creates a model that includes
    # the Input layer and three Dense layers
    # adv = Model(inputs=inputs, outputs=predictions)
    # adv.compile(optimizer='Adam',
    #             loss='binary_crossentropy',
    #             metrics=['accuracy'])
    return Model(inputs,predictions)

def train_adv(samples,labels,epoch,advh5_name,dis_trainable):     
#history = LossHistory()
    adv = adv_model()
    adv.compile(optimizer='Adam',
                loss='categorical_crossentropy',
                metrics=['accuracy'])

    if dis_trainable == 'True':
        if not os.path.exists(advh5_name):  
            adv.fit(x = samples, y = labels, validation_split = 0.3,epochs=epoch,verbose=0)
            adv.save_weights(advh5_name)     
        else:
            adv.load_weights(advh5_name)
            adv.fit(x = samples, y = labels, validation_split = 0.3,epochs=epoch,verbose=0)
            adv.save_weights(advh5_name)
    else:
        if not os.path.exists(advh5_name): 
            print("currenctly there is no pre-trained discriminator, you may train one first or load from others")
            dis_trainable = 'True'
            print('begin to train adv')
            adv.fit(x = samples, y = labels, validation_split = 0.3,epochs=100,verbose=0)
            adv.save_weights(advh5_name) 
            
        
        adv.load_weights(advh5_name)
        
        # predict 
        loss,accuracy = adv.evaluate(samples,labels,verbose=0)
        print('loss:',loss,'accuracy:',accuracy)
        return loss





