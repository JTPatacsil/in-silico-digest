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
from formencode import validators, compound, schema

# used to validate the search form. Make sure arguments are valid
class SearchFormValidator(schema.Schema):
    query = validators.PlainText()

    enzyme = validators.OneOf(["Tryspin", "Arg","Lys N"])

    min_l = validators.Int(min = 0)
    max_l = validators.Int(min = 0)
    min_w = validators.Int(min = 0)
    max_w = validators.Int(min = 0)

    misses = validators.Int(min = 0)


class SearchForm(twf.Form):
    class child(twf.TableLayout):
        # args
        query  = twf.TextArea(label = "Search Term",
                              rows = 3, cols = 50)
        # file   = twf.FileField(label = "Example File")
        enzyme = twf.SingleSelectField(lable = "Enzyme",
                                       options = ["Arg C",
                                                  "Lys N"],
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
    def index(self):
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
    
    ### Adding search functionality
    @expose("digest.templates.digest")
    @validate(SearchForm, error_handler = index)
    def digest(self, query, enzyme = None, min_l = None, max_l = None,
               min_w = None, max_w = None, misses = None):

        ### FORMATTING PARAMATERS
        # sequence
        seq = Seq(UniProt_acc = "P12345")
        
        # setting enzyme to enzyme class
        if enzyme == None or enzyme == "Trypsin":
            enzyme = Trypsin()
        ### add other implemented enzymes
        elif enzyme == "Arg C":
            pass

        
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
        valid_fragments = digest(seq, enzyme, min_l, max_l, min_w, max_w, misses)

        # returning results to the webpage
        return dict(seq=seq,
                    enzyme = enzyme,
                    v_frags = valid_fragments)
    
    
