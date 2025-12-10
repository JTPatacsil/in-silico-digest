# -*- coding: utf-8 -*-
def init_model(engine):
    """Call me before using any of the tables or classes in the model."""

# Import your model modules here.
from .Seq import Seq
from .Enzyme import Enzyme
from .Fragment import Fragment
from .custom_io import access_uniprot
from .digest import get_missed_cleavages, enzyme_digest


__all__ = ["Seq","digest","Enzyme","Fragment","custom_io"]


