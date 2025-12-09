# -*- coding: utf-8 -*-
"""Main Controller"""

from tg import expose, flash, require, url, lurl
from tg import request, redirect, tmpl_context
from tg.i18n import ugettext as _, lazy_ugettext as l_
from tg.exceptions import HTTPFound

from digest.lib.base import BaseController
from digest.controllers.error import ErrorController

__all__ = ['RootController']

from digest.model.digest import *

import tw2.forms as twf
from tg import validate
from formencode import validators, compound, schema, Invalid

# used to validate the search form. Make sure arguments are valid
class SearchFormValidator(schema.Schema):
    query = compound.All(validators.PlainText(),
                         validators.Regex("^[AaCcDdEeFfGgHhIiKkLlMmNnPpQqRrSsTtVvWwXxYy0-9 ]+$"),
                         strip = True)

    enzyme = validators.OneOf(["Tryspin", "Arg C", "Asp N",
                               "Lys N", "Lys C",  "CNBr",
                               "Protein Kinase K", "Pepsin (pH 1.3)", "Pepsin (pH > 2)"])

    min_l = validators.Int(min = 0)
    max_l = validators.Int(min = 0)
    min_w = validators.Int(min = 0)
    max_w = validators.Int(min = 0)

    misses = validators.Int(min = 0)

    # This function was given by Claude
    # used to make sure that the minimums are less than the maximums
    def validate_python(self, value_dict, state):
        if value_dict.get('min_l') and value_dict.get('max_l') and int(value_dict['min_l']) >= int(value_dict['max_l']):
            raise Invalid('Minimum length must be less than maximum length', value_dict, state)
        if value_dict.get('min_w') and value_dict.get('max_w') and int(value_dict['min_w']) >= int(value_dict['max_w']):
            raise Invalid('Minimum weight must be less than maximum weight', value_dict, state)


class SearchForm(twf.Form):
    class child(twf.TableLayout):
        query  = twf.TextArea(label = "Search Term", rows = 3, cols = 50)
        enzyme = twf.SingleSelectField(lable = "Enzyme",
                                       options = ["Arg C", "Asp N",
                                                  "Lys N", "Lys C",
                                                  "CNBr", "Protein Kinase K",
                                                  "Pepsin (pH 1.3)", "Pepsin (pH >2)"],
                                       prompt_text = "Trypsin")
        
        min_l = twf.TextField(label = "Min Fragment length")
        max_l = twf.TextField(label = "Max Fragment length")
        min_w = twf.TextField(label = "Min Fragment Weight (Da)")
        max_w = twf.TextField(label = "Max Fragment Weight (Da)")

        misses = twf.SingleSelectField(label = "Number of missed cleavages allowed",
                                       options = ["1","2","3","4","5"],
                                       prompt_text = "0")


        # for presentation only...
        css_class = 'table'
        attrs = {'style': 'width: 600px;'}
    
    action = '/digest'
    submit = twf.SubmitButton(value="Search")
    validator = SearchFormValidator

class RootController(BaseController):
    ### Navigating the websites
    @expose('digest.templates.index')
    def index(self, **kwargs):
        """ Handle the front page """
        return dict(form = SearchForm # form for data
                    )

    @expose('digest.templates.enzymes')
    def enzymes(self):
        """Enzyme information page"""
        return dict(page = 'enzymes')

    @expose('digest.templates.about')
    def about(self):
        """Project Information Page"""
        return dict(page = 'about')
    
    ### Adding digest functionality
    @expose("digest.templates.digest")
    @validate(SearchForm, error_handler = index)
    def digest(self, query, enzyme = None, min_l = None, max_l = None,
               min_w = None, max_w = None, misses = None):

        ### FORMATTING PARAMATERS
        # sequence
        if any(c in "0123456789" for c in query): #Checking to see if it has a UniProt ID
            seq = Seq(UniProt_acc = query)
        else:
            seq = Seq(seq = query, name = "User Provided Protein")
        
        # setting enzyme to enzyme class
        if enzyme == None or enzyme == "Trypsin":
            enzyme = Trypsin()
        elif enzyme == "Arg C":
            enzyme = ArgC()
        elif enzyme == "Asp N":
            enzyme = AspN()
        elif enzyme == "Lys N":
            enzyme = Lys_n()
        elif enzyme == "Lys C":
            enzyme = Lys_c()
        elif enzyme == "CNBr":
            enzyme = CNBr()
        elif enzyme == "Protein Kinase K":
            enzyme = PtKinase_K()
        elif enzyme == "Pepsin (pH 1.3)":
            enzyme = Pepsin_1_3()
        elif enzyme == "Pepsin (pH > 2)":
            enzyme = Pepsin_gt2()


        # Other paramaters
        if min_l == None:
            min_l = 0
        if max_l == None:
            max_l = 10000
        if min_w == None:
            min_w = 0
        if max_w == None:
            max_w = 1000000
        if misses == None:
            misses = 0

        #### running the program
        valid_fragments = enzyme_digest(seq, enzyme, min_l, max_l, min_w, max_w, misses)

        # returning results to the webpage
        return dict(seq=seq,
                    enzyme = enzyme,
                    v_frags = valid_fragments)
    
    


