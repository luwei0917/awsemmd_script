#  Copyright 2016 The TensorFlow Authors. All Rights Reserved.
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
"""DNNRegressor with custom input_fn for Housing dataset."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import itertools
import argparse
import pandas as pd
import tensorflow as tf

# parser = argparse.ArgumentParser(description="This is my TensorFlow program")
# parser.add_argument("protein", help="The name of the protein")
# args = parser.parse_args()
tf.logging.set_verbosity(tf.logging.INFO)

# COLUMNS = ["crim", "zn", "indus", "nox", "rm", "age",
        #    "dis", "tax", "ptratio", "medv"]
FEATURES = ["Rama", "DSSP", "P_AP", "Water", "Burial",
            "Helix", "Frag_Mem", "QGO", "VTotal", "Rw"]
LABEL = "Good"


def input_fn(data_set):
    continuous_cols = {k: tf.constant(data_set[k].values, shape=[data_set[k].size,1]) for k in FEATURES}
    feature_cols = dict(continuous_cols)
    label_keys = tf.constant(data_set[LABEL].values, shape=[data_set[LABEL].size,1])
    return feature_cols, label_keys


def main(unused_argv):
    # Load datasets
    training_set = pd.read_csv("~/Research/data/tensorFlow/train.csv", skipinitialspace=True)
    test_set = pd.read_csv("~/Research/data/tensorFlow/test.csv", skipinitialspace=True)

    # Set of 6 examples for which to predict median house values
    # name = args.protein
    # prediction_set = pd.read_csv("~/Research/data/tensorFlow/{}.csv".format(name), skipinitialspace=True)
    # prediction_set = pd.read_csv("~/Research/data/tensorFlow/test.csv", skipinitialspace=True)
    # Feature cols
    feature_cols = [tf.contrib.layers.real_valued_column(k)
                  for k in FEATURES]

    # Build 2 layer fully connected DNN with 10, 10 units respectively.
    classifier = tf.contrib.learn.DNNClassifier(feature_columns=feature_cols,
                                            hidden_units=[10, 20, 10],
                                            n_classes=2,
                                            model_dir="/tmp/awsem_classifier")

    # Fit
    classifier.fit(input_fn=lambda: input_fn(training_set), steps=5000)

    # Score accuracy
    # ev = classifier.evaluate(input_fn=lambda: input_fn(test_set), steps=1)
    # loss_score = ev["loss"]
    # print("Loss: {0:f}".format(loss_score))

    # accuracy_score = classifier.evaluate(input_fn=lambda: input_fn(test_set))["accuracy"]
    # print('Accuracy: {0:f}'.format(accuracy_score))
    # Print out predictions
    # y = regressor.predict(input_fn=lambda: input_fn(prediction_set))
    # y = classifier.predict(input_fn=lambda: input_fn(prediction_set), as_iterable=True)
    #y = classifier.predict(input_fn=lambda: input_fn(prediction_set))
    # .predict() returns an iterator; convert to a list and print predictions
    # predictions = list(itertools.islice(y, 6))
    name_list = ["1MBA", "T251", "T0766", "T0784", "T0792", "T0803", "T0815", "T0833"]
    for name in name_list:
        prediction_set = pd.read_csv("~/Research/data/tensorFlow/{}.csv".format(name), skipinitialspace=True)
        y = classifier.predict_classes(input_fn=lambda: input_fn(prediction_set))
        predictions = list(y)
        with open("/Users/weilu/Research/data/tensorFlow/{}_results.csv".format(name), "w") as f:
            f.write("Result\n")
            for i in predictions:
                f.write(str(i) + "\n")

    # with open("/Users/weilu/Research/data/tensorFlow/{}_results.csv".format(name), "w") as f:
    #     f.write("Result\n")
    #     for i in predictions:
    #         f.write(str(i) + "\n")
if __name__ == "__main__":
    tf.app.run()
