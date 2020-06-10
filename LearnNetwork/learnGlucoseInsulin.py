import tensorflow as tf
import numpy as np
import sys
import random
import scipy.io

# Taisa Kushner     taisa.kushner@colorado.edu
# learn NN code
# example command line arg to run:
# python3 learnGlucoseInsulin.py TensorFlowInput_DallaMan_ptid1_first250.csv 37 37 DallaMan_pt1

# LearningData is a class that maintains the data from the CSV file.
# The format of the CSV file should be as follows
# No header, use commas to separate values
# Insulin(t-180),  Insulin(t-175), ... ,Insulin(t), Glucose(t-180), .. Glucose(t), Glucose(t+Tg)
# Specify about 20% of the file  size as the test_set. This will be used to evaluate accuracy at the very end.

# Please reference:
# S. Dutta*, T. Kushner*, S. Sankaranarayanan, Robust Data-Driven Control of Artificial
# Pancreas Systems using Neural Networks. Computational Methods in Systems Biology.
# (CMSB). Lecture Notes in Computer Science, vol 11095. (2018) Springer

class LearningData:
    # Class Constructor
    def __init__(self, csv_filename, num_insulin_inputs, num_gluc_inputs, test_set_size):
        # Load the data from the csv_filename
        self.__all_rows = [] # This is a list of all the contents of the CSV file
        self.__test_size = test_set_size
        self.__num_inputs = num_insulin_inputs + num_gluc_inputs
        self.__num_insulin_inputs = num_insulin_inputs
        self.__num_gluc_inputs = num_gluc_inputs
        # Read the file and store the data into self.__all_rows list
        with open(csv_filename, 'r') as fhandle:
            for line in fhandle:
                lst = line.split(',')
                row = []
                for x in lst:
                    row.append(float(x))
                self.__all_rows.append(row)
                # print(row)
            fhandle.close()

    def get_num_inputs(self):
        return self.__num_inputs

    # Randomly select "batch_size" number of rows from the CSV and return them in two
    # matrices: insulin_x_values, gluc_x_values and outputs stored in y_values
    def get_random_inputs_and_outputs(self, batch_size):
        num_inputs = self.__num_inputs
        num_insulin_inputs = self.__num_insulin_inputs
        num_gluc_inputs = self.__num_gluc_inputs
        n = len(self.__all_rows) - self.__test_size
        # choose batch_size random indices from 0 to n-1
        lst = random.sample(range(n), batch_size)
        insulin_x_values = np.zeros(shape=(num_insulin_inputs, batch_size))
        gluc_x_values = np.zeros(shape=(num_gluc_inputs, batch_size))
        y_values = np.zeros(shape=(1, batch_size))
        count = 0
        for i in lst:
            row_lst = self.__all_rows[i]
            insulin_x_values[:, count] = np.matrix(row_lst[0:num_insulin_inputs]).reshape(1, num_insulin_inputs)
            gluc_x_values[:, count] = np.matrix(row_lst[num_insulin_inputs:num_inputs]).reshape(1, num_gluc_inputs)
            y_values[0, count] = np.matrix([row_lst[num_inputs]])
            count = count + 1
        return insulin_x_values, gluc_x_values, y_values

    # Get the test set for evalating accuracy of the final model
    def get_test_set(self):
        num_inputs = self.__num_inputs
        n = len(self.__all_rows)
        count = 0
        batch_size = self.__test_size
        insulin_x_values = np.zeros(shape=(num_insulin_inputs, batch_size))
        gluc_x_values = np.zeros(shape=(num_gluc_inputs, batch_size))
        y_values = np.zeros(shape=(1, batch_size))
        count = 0
        for i in range(n - self.__test_size, n):
            row_lst = self.__all_rows[i]
            insulin_x_values[:, count] = np.matrix(row_lst[0:num_insulin_inputs]).reshape(1, num_insulin_inputs)
            gluc_x_values[:, count] = np.matrix(row_lst[num_insulin_inputs:num_inputs]).reshape(1, num_gluc_inputs)
            y_values[0, count] = np.matrix([row_lst[num_inputs]])
            count = count + 1
        return insulin_x_values, gluc_x_values, y_values


# This is a worker function that sets up a single hidden layer of the network.
# DO not call directly
def setup_multilayer_network(num_inputs, inp_vars, num_outputs, layer_sizes, make_output_layer):
    n0 = num_inputs
    prev_layer_in = inp_vars
    weights_list = []
    for n in layer_sizes:
        W1 = tf.Variable(tf.random_normal(shape=(n, n0), mean=0.0, stddev=1.0, dtype=tf.float32))
        b1 = tf.Variable(tf.random_normal(shape=(n, 1), mean=0.1, stddev=0.5, dtype=tf.float32))
        layer1_in = tf.matmul(W1, prev_layer_in) + b1
        layer1_out = tf.nn.relu(layer1_in)
        weights_list.append((W1, b1))
        prev_layer_in = layer1_out
        n0 = n
    if make_output_layer:
        # Output Layer
        W3 = tf.Variable(tf.random_normal(shape=(num_outputs, n0), mean=0.0, stddev=1.0, dtype=tf.float32))
        b3 = tf.Variable(tf.random_normal(shape=(num_outputs, 1), mean=0.1, stddev=0.5, dtype=tf.float32))
        weights_list.append((W3, b3))
        y_mdl = tf.matmul(W3, prev_layer_in) + b3
    else:
        y_mdl = prev_layer_in
    return y_mdl, weights_list



def setup_insulin_glucose_nn(fn, insulin_layer_sizes, gluc_layer_sizes, combined_layer_sizes,
                            num_insulin_inputs, num_gluc_inputs, output_fstem):


    # First setup the input place holders
    # Placeholder for insulin input
    insulin_x = tf.placeholder(shape=(num_insulin_inputs, None), dtype=tf.float32)
    # Placeholder for glucose input
    gluc_x = tf.placeholder(shape=(num_gluc_inputs, None), dtype=tf.float32)
    # Place holder for the overall output prediction
    y = tf.placeholder(shape=(1, None), dtype=tf.float32)
    # Set up the network for insulin
    insulin_out, weights_list_insulin = setup_multilayer_network(num_insulin_inputs, insulin_x, 0, insulin_layer_sizes, False)
    # Set up the glucose network
    gluc_out, weights_list_gluc = setup_multilayer_network(num_gluc_inputs, gluc_x, 0, gluc_layer_sizes, False)
    # Now make a combined layer with insulin and gluc outputs
    joint_layer_in = tf.concat([gluc_out, insulin_out], 0)
    joint_layer_size = joint_layer_in.get_shape().as_list()
    #print('joint_layer_size', joint_layer_size)
    # Setup the overall network
    y_mdl, weights_list_joint = setup_multilayer_network(joint_layer_size[0], joint_layer_in,
                                                         1, combined_layer_sizes, True)
    # This is the loss function -- mean
    loss_fn = tf.reduce_mean(tf.square(y - y_mdl))
    # example loss fn for lower quantile, flip y & y_mdl for upper 
    #loss_fn = tf.reduce_max(tf.maximum(tf.scalar_mul(0.7,(y_mdl-y)), tf.scalar_mul(0.3,(y-y_mdl))))


    # Setting up the training
    opt = tf.train.AdamOptimizer(1e-4)
    train = opt.minimize(loss_fn)

    # DO the batch training: Stochastic gradient descent or backpropagation
    # Parameters for training
    num_training_steps = 100000 #This is the total number of training steps taken
    batch_size = 200  #This is the size of each batch in the training
    with tf.Session() as sess:
        sess.run(tf.global_variables_initializer())
        for i in range(num_training_steps): # train for 100000 steps -- hopefully it is sufficient
            # Get a random batch of batch_size inputs
            batch_X_insulin, batch_X_gluc, batch_Y = fn.get_random_inputs_and_outputs(batch_size)
            # Run training -- appropriate data added to placeholders
            train.run(feed_dict={insulin_x: batch_X_insulin, gluc_x: batch_X_gluc, y: batch_Y})

            if i % 1000 == 0: # every 1000 steps, let us evaluate accuracy on the training data itself.
                test_X_insulin, test_X_gluc, test_Y = fn.get_random_inputs_and_outputs(10*batch_size)
                test_accuracy = loss_fn.eval(feed_dict={insulin_x: test_X_insulin, gluc_x: test_X_gluc, y: test_Y})
                print('Step: ', i, 'Accuracy on training data samples:', test_accuracy)
        # DONE: Let s do the test
        final_test_X_insulin, final_test_X_gluc , final_test_Y = fn.get_test_set()
        test_accuracy = loss_fn.eval(feed_dict={insulin_x: final_test_X_insulin, gluc_x: final_test_X_gluc, y: final_test_Y})
        print('Final Accuracy on TEST data:', test_accuracy)
        md = {}
        # Storing it in a .mat file for loading from matlab for further simulations
        i = 1
        for (W,b) in weights_list_insulin:
            md['W_insulin_%d'%(i)] = W.eval(sess)
            md['b_insulin_%d'%(i)] = b.eval(sess)
            i = i + 1
        i = 1
        for (W, b) in weights_list_gluc:
            md['W_gluc_%d'%(i)] = W.eval(sess)
            md['b_gluc_%d'%(i)] = b.eval(sess)
            i = i + 1
        i = 1
        for (W, b) in weights_list_joint:
            md['W_joint_%d'%(i)] = W.eval(sess)
            md['b_joint_%d'%(i)] = b.eval(sess)
            i = i + 1
        scipy.io.savemat('%s.mat' % output_fstem, md)
        print('Results saved to file: %s' % output_fstem)


if __name__ == '__main__':
    if len(sys.argv) < 5:
        print('Usage: ', sys.argv[0], '[name of csv] [num_insulin_inputs] [num_gluc_inputs] [outputfilestem]')
        sys.exit(2)

    csv_filename = sys.argv[1]
    num_insulin_inputs = int(sys.argv[2])
    num_gluc_inputs = int(sys.argv[3])
    output_fstem = sys.argv[4]
    # SETUP the number of layers and sizes for each network
    insulin_layer_sizes = [4] # Insulin network is a single layer with 4 neurons
    gluc_layer_sizes = [4] # Glucose network is a single layer with 4 neurons
    combined_layer_sizes = [8] # The combined network is a single layer with 8 neurons.

    fn = LearningData(csv_filename, num_insulin_inputs, num_gluc_inputs, 6000) # This will read the CSV file and have a test set of size 6000 -- use 20% of data
    setup_insulin_glucose_nn(fn, insulin_layer_sizes, gluc_layer_sizes,
                            combined_layer_sizes, num_insulin_inputs,
                            num_gluc_inputs, output_fstem)
