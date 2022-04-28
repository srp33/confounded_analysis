import numpy as np
import random
from sklearn import preprocessing
from sklearn.datasets import make_classification
from sklearn.model_selection import cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
import sys

out_file_path = sys.argv[1]

def relU(x) :
  return np.maximum(0.3*x, x)

def leakyrelU(x) :
  return np.maximum(0.3*x, x)

class Noise(object):
    def __init__(
            self, shape, order=2, discount_factor=0.01, activation=relU, min_= -10.0, max_=10.0
        ):
        self._layers = [
            self._Layer(shape, discount_factor, activation, min_, max_)
            for _ in range(order)
        ]

    def adjust(self, x):
        y = x
        for layer in self._layers:
            y = layer.adjust(y)
        return y
    class _Layer(object):
        def __init__(self, shape, discount_factor, activation=relU, min_=-10.0, max_=10.0):
            self.weights = np.random.normal(size=shape)
            self.bias = np.random.normal()
            self.discount_factor = discount_factor
            self.activation = activation
            self.max = max_
            self.min = min_

        def adjust(self, array):
            adjusted = array + (
                self.activation(
                    array * self.weights + self.bias
                ) * self.discount_factor
            )
            return self.threshold(adjusted)

        def threshold(self, array):
            too_low = array < self.min
            too_high = array > self.max
            array[too_low] = self.min
            array[too_high] = self.max
            return array

n_samples = 200 # Must be an even number
n_random_features = 800
#n_random_features = 19800
n_informative_features = 25
n_redundant_features = 175
n_features = n_random_features + n_informative_features + n_redundant_features
# This controls the ability to differentiate between the batches and class labels.
class_sep = 0.0

random_state = 0
random.seed(random_state)

n1 = Noise((n_features,), order=15, discount_factor=0.04, activation=relU, min_= -8.0, max_=8.0) #batch
n2 = Noise((n_features,), order=15, discount_factor=0.04, activation=relU, min_= -8.0, max_=8.0) #batch
n3 = Noise((n_features,), order=15, discount_factor=0.07, activation=leakyrelU) #labels
n4 = Noise((n_features,), order=15, discount_factor=0.07, activation=leakyrelU) #labels

X, labels = make_classification(n_samples, n_features, n_informative_features, n_redundant_features, n_clusters_per_class=1, random_state=random_state, shuffle=False, flip_y=0.0, class_sep = class_sep)
batches = np.array([0 for i in range(int(n_samples / 4))] + [1 for i in range(int(n_samples / 4))] + [0 for i in range(int(n_samples / 4))] + [1 for i in range(int(n_samples / 4))])

random_labels = labels.copy()
random_batches = batches.copy()
np.random.shuffle(random_labels)
np.random.shuffle(random_batches)

batch = []
images = []
imagesforlabel = []

#Make Noise Matrix
i = 0
for row in X:
    add_noise = np.random.uniform(-0.15, 0.15, size=(n_features,)) #batch
    add_noise1 = np.random.uniform(-0.1, 0.1, size=(n_features,))
    mult_noise = np.random.uniform(1.85, 1.95, size=(n_features,)) #batch
    mult_noise1 = np.random.uniform(0.85, 0.95, size=(n_features,)) #labels

    image = row.tolist()
    image2 = row.tolist()
    image = (image + add_noise) * mult_noise
    image2 = (image2 + add_noise) * mult_noise1

    if batches[i] == 0:
        images.append(n1.adjust(image))
    else:
        images.append(n2.adjust(image))
    if labels[i] == 0:
        imagesforlabel.append(n3.adjust(image2))
    else:
        imagesforlabel.append(n4.adjust(image2))

    i+= 1

images = np.array(images)
imagesforlabel = np.array(imagesforlabel)

X = images + imagesforlabel

X = preprocessing.scale(X)

# Verify that we can differentiate the classes and batches
clf = RandomForestClassifier(max_depth=5, n_estimators=100, random_state=random_state)

print("Predicting true labels:")
print(round(np.mean(cross_val_score(clf, X, labels, cv=5)), 3))
print("Predicting randomized labels:")
print(round(np.mean(cross_val_score(clf, X, random_labels, cv=5)), 3))
print("Predicting batches:")
print(round(np.mean(cross_val_score(clf, X, batches, cv=5)), 3))
print("Predicting randomized batches:")
print(round(np.mean(cross_val_score(clf, X, random_batches, cv=5)), 3))

with open(out_file_path, 'w') as out_file:
    out_file.write(",".join(["ID", "Class", "Batch"] + ["Gene%i" % i for i in range(1, n_features + 1)]) + "\n")

    i = 0
    for row in X:
        batch = str(batches[i])
        label = str(labels[i])
        out_file.write(",".join(["Sample{}".format(i + 1), label, batch] + [str(x) for x in row]) + "\n")
        i += 1
