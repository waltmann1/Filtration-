import numpy as np
import copy
class KMeans():
    def __init__(self, n_clusters, distance_measure="euclidean"):
        """
        This class implements the traditional KMeans algorithm with hard assignments:

        https://en.wikipedia.org/wiki/K-means_clustering

        The KMeans algorithm has two steps:

        1. Update assignments
        2. Update the means

        While you only have to implement the fit and predict functions to pass the
        test cases, we recommend that you use an update_assignments function and an
        update_means function internally for the class.

        Use only numpy to implement this algorithm.

        Args:
            n_clusters (int): Number of clusters to cluster the given data into.

        """
        self.distance_measure = self.euclidean_distance_measure

        if distance_measure == "sqrt":
            self.distance_measure = self.sqrt_distance_measure

        if distance_measure == "convolution":
            self.distance_measure = self.convolutional_distance_measure

        if distance_measure == "cover":
            self.distance_measure = self.covered_unconvered_distance_measure

        self.n_clusters = n_clusters
        self.means = None
        self.features = None
        self.labels = None

    def fit(self, features):
        """
        Fit KMeans to the given data using `self.n_clusters` number of clusters.
        Features can have greater than 2 dimensions.

        Args:
            features (np.ndarray): array containing inputs of size
                (n_samples, n_features).
        Returns:
            None (saves model - means - internally)
        """
        self.labels = np.zeros(features.shape[0])
        store_labels = np.ones(features.shape[0])

        distances = np.zeros((features.shape[0], self.n_clusters))

        average = np.average(np.array(features), axis=0)
        self.means = np.array([average for _ in range(self.n_clusters)])

        #self.means = np.array([features[i] for i in range(self.n_clusters)])
        self.means = np.add(self.means, np.multiply(1, np.subtract(np.random.rand(self.n_clusters, features.shape[1]), .5)))
        self.means = np.where(self.means < 0, 0, self.means)
        #print(self.means.shape)
        #print(self.means[0][:10])
        #quit()

        self.features = features
        iter_count = 0
        while not np.all(store_labels == self.labels) and iter_count < 10:
            #print("loop")
            #print(self.labels)
            store_labels = copy.deepcopy(self.labels)
            for i in range(features.shape[0]):
                for j in range(self.n_clusters):
                    #istances[i][j] = np.linalg.norm(np.subtract(features[i], self.means[j]))
                    distances[i][j] = self.distance_measure(features[i], self.means[j])
                self.labels[i] = np.argmin(distances[i])

            #store_means = self.means
            #print(np.all(store_labels == self.labels))
            self.means = [np.average([features[i] for i in range(features.shape[0]) if self.labels[i] == k], axis=0) for k in range(len(self.means))]
            #print(np.all(store_means[7] == self.means[7]))
            iter_count = iter_count +1


    def predict(self, features):
        """
        Given features, an np.ndarray of size (n_samples, n_features), predict cluster
        membership labels.

        Args:
            features (np.ndarray): array containing inputs of size
                (n_samples, n_features).
        Returns:
            predictions (np.ndarray): predicted cluster membership for each features,
                of size (n_samples,). Each element of the array is the index of the
                cluster the sample belongs to.
        """

        distances = np.zeros((features.shape[0], self.n_clusters))
        labels = np.zeros(features.shape[0])

        for i in range(features.shape[0]):
            for j in range(self.n_clusters):
                distances[i][j] = np.linalg.norm(np.subtract(features[i], self.means[j]))
            labels[i] = np.argmin(distances[i])

        return labels


    def euclidean_distance_measure(self, array_1, array_2):

        return np.array(np.linalg.norm(np.subtract(array_1, array_2)))

    def sqrt_distance_measure(self, array_1, array_2):

        temp_1 = np.sqrt(array_1)

        temp_2 = np.sqrt(array_2)

        return self.euclidean_distance_measure(temp_1, temp_2)

    def convolutional_distance_measure(self, array_1, array_2, window_size=5):

        temp_1 = np.zeros(len(array_1) - window_size)
        temp_2 = np.zeros(len(array_2) - window_size)
        for i in range(len(array_1) - window_size):
            temp_1[i] = np.sum(array_1[i:i+5])
            temp_2[i] = np.sum(array_2[i:i+5])

        return self.euclidean_distance_measure(temp_1, temp_2)

    def covered_unconvered_distance_measure(self, array_1, array_2):

        return self.euclidean_distance_measure(
            self.zero_one_filter(array_1), self.zero_one_filter(array_2))

    def zero_one_filter(self, array, cut=1):

        return np.where(np.array((array)) >= cut, 1, 0)

    def rms_error(self):

        #print(self.features.shape)
        #print(self.labels)
        total = 0
        for ind, feature in enumerate(self.features):
            error = self.distance_measure(feature, self.means[int(self.labels[ind])])
            #print(error)
            total += error * error
        return np.sqrt(total)


    #def write_means(self):

     #   for ind, mean in enumerate(self.means):



