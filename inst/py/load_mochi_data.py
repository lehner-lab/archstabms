#!/usr/bin/env python

import pandas as pd

import bz2
import pickle
import _pickle as cPickle

# import torch
# from torch import Tensor
# from torch.nn import functional as F

import pymochi
from pymochi.data import *


def load_mochi_data(file):
    with bz2.BZ2File(file, 'rb') as f:
        load_dict = cPickle.load(f)
    return {
        'model_design' : load_dict['data'].model_design,
        'max_interaction_order' : load_dict['data'].max_interaction_order,
        'vtable' : load_dict['data'].fdata.vtable,
        'sequenceType' : load_dict['data'].fdata.sequenceType,
        'variantCol' : load_dict['data'].fdata.variantCol,
        'mutationOrderCol' : load_dict['data'].fdata.mutationOrderCol,
        'wildtype' : load_dict['data'].fdata.wildtype,
        'wildtype_split' : load_dict['data'].fdata.wildtype_split,
        'additive_trait_names' : load_dict['data'].additive_trait_names,
        'phenotype_names' : load_dict['data'].phenotype_names,
        'fitness' : load_dict['data'].fitness,
        'phenotypes' : load_dict['data'].phenotypes,
        'Xohi' : load_dict['data'].Xohi,
        'cvgroups' : load_dict['data'].cvgroups,
        'coefficients' : load_dict['data'].coefficients,
        'coefficients_userspec' : load_dict['data'].coefficients_userspec,
        'feature_names' : load_dict['data'].feature_names}