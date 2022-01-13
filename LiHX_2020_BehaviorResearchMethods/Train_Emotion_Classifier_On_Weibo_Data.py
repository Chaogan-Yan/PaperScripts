#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 08:54:59 2019

@author: -Bin Lu lub@psych.ac.cn
Adopted from https://github.com/strongio/keras-bert/blob/master/keras-bert.py
See https://towardsdatascience.com/bert-in-keras-with-tensorflow-hub-76bcbc9417b for more information
"""


#!pip install pandas
#!pip install tensorflow_hub
#!pip install bert-tensorflow
#!pip install tqdm

import tensorflow as tf
import pandas as pd
import tensorflow_hub as hub
import os
import re
import numpy as np
from bert.tokenization import FullTokenizer
from tqdm import tqdm
from tensorflow.keras import backend as K
import scipy.io as sio
from keras.utils import multi_gpu_model
from keras.utils.np_utils import to_categorical


# Initialize session
sess = tf.Session()


class PaddingInputExample(object):
    """Fake example so the num input examples is a multiple of the batch size.
  When running eval/predict on the TPU, we need to pad the number of examples
  to be a multiple of the batch size, because the TPU requires a fixed batch
  size. The alternative is to drop the last batch, which is bad because it means
  the entire output data won't be generated.
  We use this class instead of `None` because treating `None` as padding
  battches could cause silent errors.
  """


class InputExample(object):
    """A single training/test example for simple sequence classification."""

    def __init__(self, guid, text_a, text_b=None, label=None):
        """Constructs a InputExample.
    Args:
      guid: Unique id for the example.
      text_a: string. The untokenized text of the first sequence. For single
        sequence tasks, only this sequence must be specified.
      text_b: (Optional) string. The untokenized text of the second sequence.
        Only must be specified for sequence pair tasks.
      label: (Optional) string. The label of the example. This should be
        specified for train and dev examples, but not for test examples.
    """
        self.guid = guid
        self.text_a = text_a
        self.text_b = text_b
        self.label = label


def create_tokenizer_from_hub_module(bert_path):
    """Get the vocab file and casing info from the Hub module."""
    bert_module = hub.Module(bert_path)
    tokenization_info = bert_module(signature="tokenization_info", as_dict=True)
    vocab_file, do_lower_case = sess.run(
        [tokenization_info["vocab_file"], tokenization_info["do_lower_case"]]
    )

    return FullTokenizer(vocab_file=vocab_file, do_lower_case=do_lower_case)


def convert_single_example(tokenizer, example, max_seq_length=256):
    """Converts a single `InputExample` into a single `InputFeatures`."""

    if isinstance(example, PaddingInputExample):
        input_ids = [0] * max_seq_length
        input_mask = [0] * max_seq_length
        segment_ids = [0] * max_seq_length
        label = 0
        return input_ids, input_mask, segment_ids, label

    tokens_a = tokenizer.tokenize(example.text_a)
    if len(tokens_a) > max_seq_length - 2:
        tokens_a = tokens_a[0 : (max_seq_length - 2)]

    tokens = []
    segment_ids = []
    tokens.append("[CLS]")
    segment_ids.append(0)
    for token in tokens_a:
        tokens.append(token)
        segment_ids.append(0)
    tokens.append("[SEP]")
    segment_ids.append(0)

    input_ids = tokenizer.convert_tokens_to_ids(tokens)

    # The mask has 1 for real tokens and 0 for padding tokens. Only real
    # tokens are attended to.
    input_mask = [1] * len(input_ids)

    # Zero-pad up to the sequence length.
    while len(input_ids) < max_seq_length:
        input_ids.append(0)
        input_mask.append(0)
        segment_ids.append(0)

    assert len(input_ids) == max_seq_length
    assert len(input_mask) == max_seq_length
    assert len(segment_ids) == max_seq_length

    return input_ids, input_mask, segment_ids, example.label


def convert_examples_to_features(tokenizer, examples, max_seq_length=256):
    """Convert a set of `InputExample`s to a list of `InputFeatures`."""

    input_ids, input_masks, segment_ids, labels = [], [], [], []
    for example in tqdm(examples, desc="Converting examples to features"):
        input_id, input_mask, segment_id, label = convert_single_example(
            tokenizer, example, max_seq_length
        )
        input_ids.append(input_id)
        input_masks.append(input_mask)
        segment_ids.append(segment_id)
        labels.append(label)
    return (
        np.array(input_ids),
        np.array(input_masks),
        np.array(segment_ids),
        np.array(labels).reshape(-1, 1),
    )


def convert_text_to_examples(texts, labels):
    """Create InputExamples"""
    InputExamples = []
    for text, label in zip(texts, labels):
        InputExamples.append(
            InputExample(guid=None, text_a=" ".join(text), text_b=None, label=label)
        )
    return InputExamples


class BertLayer(tf.keras.layers.Layer):
    def __init__(
        self,
        n_fine_tune_layers=10,
        pooling="mean",
        # The default bert_path used to be a online package, changed to a local package
        bert_path="/mnt/Data3/RfMRILab/Lubin/ContainerVolumes/Bert/Model/bert_chinese_L-12_H-768_A-12",
        **kwargs
    ):
        self.n_fine_tune_layers = n_fine_tune_layers
        self.trainable = True
        self.output_size = 768
        self.pooling = pooling
        self.bert_path = bert_path
        if self.pooling not in ["first", "mean"]:
            raise NameError(
                "Undefined pooling type (must be either first or mean, but is {self.pooling}"
            )

        super(BertLayer, self).__init__(**kwargs)

    def build(self, input_shape):
        self.bert = hub.Module(
            self.bert_path, trainable=self.trainable, name="{}_module".format(self.name)
        )

        # Remove unused layers
        trainable_vars = self.bert.variables
        if self.pooling == "first":
            trainable_vars = [var for var in trainable_vars if not "/cls/" in var.name]
            trainable_layers = ["pooler/dense"]

        elif self.pooling == "mean":
            trainable_vars = [
                var
                for var in trainable_vars
                if not "/cls/" in var.name and not "/pooler/" in var.name
            ]
            trainable_layers = []
        else:
            raise NameError(
                "Undefined pooling type (must be either first or mean, but is {self.pooling}"
            )

        # Select how many layers to fine tune
        for i in range(self.n_fine_tune_layers):
            trainable_layers.append("encoder/layer_{str(11 - i)}")

        # Update trainable vars to contain only the specified layers
        trainable_vars = [
            var
            for var in trainable_vars
            if any([l in var.name for l in trainable_layers])
        ]

        # Add to trainable weights
        for var in trainable_vars:
            self._trainable_weights.append(var)

        for var in self.bert.variables:
            if var not in self._trainable_weights:
                self._non_trainable_weights.append(var)

        super(BertLayer, self).build(input_shape)

    def call(self, inputs):
        inputs = [K.cast(x, dtype="int32") for x in inputs]
        input_ids, input_mask, segment_ids = inputs
        bert_inputs = dict(
            input_ids=input_ids, input_mask=input_mask, segment_ids=segment_ids
        )
        if self.pooling == "first":
            pooled = self.bert(inputs=bert_inputs, signature="tokens", as_dict=True)[
                "pooled_output"
            ]
        elif self.pooling == "mean":
            result = self.bert(inputs=bert_inputs, signature="tokens", as_dict=True)[
                "sequence_output"
            ]

            mul_mask = lambda x, m: x * tf.expand_dims(m, axis=-1)
            masked_reduce_mean = lambda x, m: tf.reduce_sum(mul_mask(x, m), axis=1) / (
                    tf.reduce_sum(m, axis=1, keepdims=True) + 1e-10)
            input_mask = tf.cast(input_mask, tf.float32)
            pooled = masked_reduce_mean(result, input_mask)
        else:
            raise NameError("Undefined pooling type (must be either first or mean, but is {self.pooling}")

        return pooled

    def compute_output_shape(self, input_shape):
        return (input_shape[0], self.output_size)


# Build model
def build_model(max_seq_length):
    in_id = tf.keras.layers.Input(shape=(max_seq_length,), name="input_ids")
    in_mask = tf.keras.layers.Input(shape=(max_seq_length,), name="input_masks")
    in_segment = tf.keras.layers.Input(shape=(max_seq_length,), name="segment_ids")
    bert_inputs = [in_id, in_mask, in_segment]

    bert_output = BertLayer(n_fine_tune_layers=3)(bert_inputs) #n_fine_tune_layers=3
    dense = tf.keras.layers.Dense(256, activation="relu")(bert_output)
    dense = tf.keras.layers.Dense(128, activation="relu")(dense)
    dense = tf.keras.layers.Dense(64, activation="relu")(dense)
    dropout = tf.keras.layers.Dropout(0.5)(dense) ##Bin
    pred = tf.keras.layers.Dense(9, activation="softmax")(dropout)

    model = tf.keras.models.Model(inputs=bert_inputs, outputs=pred)
    model.compile(loss='sparse_categorical_crossentropy', optimizer='adam', metrics=['accuracy']) #'sparse_categorical_accuracy'
    model.summary()

    return model


def initialize_vars(sess):
    sess.run(tf.local_variables_initializer())
    sess.run(tf.global_variables_initializer())
    sess.run(tf.tables_initializer())
    K.set_session(sess)


def main():
    # Params for bert model and tokenization
    bert_path = "/mnt/Data3/RfMRILab/Lubin/ContainerVolumes/Bert/Model/bert_chinese_L-12_H-768_A-12"
    max_seq_length = 256


    # Load data
    Content = '/mnt/Data3/RfMRILab/Lubin/ContainerVolumes/Bert/Data/Content_Banl.mat'
    Content = sio.loadmat(Content)
    Content = Content['ContentBanl'] 
    Label = '/mnt/Data3/RfMRILab/Lubin/ContainerVolumes/Bert/Data/Emotion_1_Num_Banl.mat'
    Label = sio.loadmat(Label)
    Label = Label['EmoBanl'] 
    #Label = to_categorical(Label)
    
    train_text = Content[0:20000,0]
    train_label = Label[0:20000,0]
    test_text = Content[20000:24000,0]
    test_label = Label[20000:24000,0]

#    UniList = np.unique(train_label)
    UniList = (1,2,3,4,5,6,7,8)  # manually assign for later class_weight assign
    class_weight = {1:0,2:0,3:0,4:0,5:0,6:0,7:0,8:0}
    nEmo = np.zeros(shape=(len(UniList),1)) 
    for i in range(len(UniList)):
        index = np.where(train_label==UniList[i])
        nEmo[i,0] = len(index[0])
    for i in range(len(UniList)):
        class_weight[UniList[i]] = nEmo[5,0]/nEmo[i,0]
            
    # Instantiate tokenizer
    tokenizer = create_tokenizer_from_hub_module(bert_path)

    # Convert data to InputExample format
    train_examples = convert_text_to_examples(train_text, train_label)
    test_examples = convert_text_to_examples(test_text, test_label)

    # Convert to features
    (
        train_input_ids,
        train_input_masks,
        train_segment_ids,
        train_labels,
    ) = convert_examples_to_features(
        tokenizer, train_examples, max_seq_length=max_seq_length
    )
    (
        test_input_ids,
        test_input_masks,
        test_segment_ids,
        test_labels,
    ) = convert_examples_to_features(
        tokenizer, test_examples, max_seq_length=max_seq_length
    )
    
#    tf.keras.backend.clear_session()
    model = build_model(max_seq_length)

    # Instantiate variables
    initialize_vars(sess)
    
    model.fit(
        [train_input_ids, train_input_masks, train_segment_ids],
        train_labels,
        validation_data=(
            [test_input_ids, test_input_masks, test_segment_ids],
            test_labels,
        ),
        epochs=50,
        batch_size=32,
        class_weight = class_weight,
        shuffle=True
    )
        
#    model.save('/mnt/Data3/RfMRILab/Lubin/ContainerVolumes/Model/EmotionPrediction.h5')

if __name__ == "__main__":
    main()
