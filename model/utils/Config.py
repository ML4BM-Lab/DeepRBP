class Config(): 
    """Class RBP settings.
    """
    def __init__(self,
                 max_node = 1024,
                 num_hidden_layers = 2,
                 learning_rate = 0.001,
                 batch_size = 128,
                 epochs = 1000,
                 same_num_nodes = False, 
                 node_divition_factor = 4,
                 activation_layer = 'relu',
                 optimizer = 'adamW',
                 source_train = 'TCGA',
                 train_tumor_types = 'all',
                 test_tumor_types = 'all'):
                
        """Setting deepsf NN characteristics.

        Args:
            max_node(int): Number of nodes in the first hidden layer of the neural network.
            num_hidden_layers(int): Number of hidden layers in the neural network.
            learning_rate (float): learning rate. Defaults to 1e-4.
            batch_size (int): The number of samples used in each iteration of training. Defaults to 128.
            epochs(int): The number of times the entire dataset is passed forward and backward through the neural network during training.
            same_num_nodes(bool): A boolean indicating whether all hidden layers except the last one should have the same number of nodes.
            node_divition_factor(int): Factor by which the number of nodes decreases from one hidden layer to the next..
            activation_layer(string): The activation function used in the hidden layers of the neural network, such as ReLU (Rectified Linear Unit).
            optimizer(string): The optimization algorithm used to update the weights of the neural network to minimize the loss function, such as AdamW.
            source_train (string): string with the name of the source we are using to the train the model. If 'all' both TCGA and GTEX samples
            are used to train the model.
            train_tumor_types(string): string with the cancer types used for training. If more than one de names must be delimited by , without spaces    
            test_tumor_types(string): string with the cancer types used for testing. If more than one de names must be delimited by , without spaces       
        """
        self.max_node = max_node
        self.num_hidden_layers = num_hidden_layers
        self.learning_rate = learning_rate
        self.batch_size = batch_size
        self.epochs = epochs
        self.same_num_nodes = same_num_nodes
        self.node_divition_factor = node_divition_factor
        self.activation_layer = activation_layer
        self.optimizer = optimizer
        self.source_train = source_train
        self.train_tumor_types = train_tumor_types
        self.test_tumor_types = test_tumor_types

    def get_config(self):
        config = dict(
            max_node = self.max_node,
            num_hidden_layers = self.num_hidden_layers,
            learning_rate = self.learning_rate,
            batch_size = self.batch_size,
            epochs = self.epochs,
            same_num_nodes = self.same_num_nodes, 
            node_divition_factor = self.node_divition_factor,
            activation_layer = self.activation_layer,
            optimizer = self.optimizer,
            source_train = self.source_train,
            train_tumor_types = self.train_tumor_types,
            test_tumor_types = self.test_tumor_types,
        )
        return config