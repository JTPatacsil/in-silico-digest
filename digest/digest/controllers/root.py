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
all_enz = set()
for e in all_enzymes:
    n = e.name
    all_enz.add(e.name)
    
    if n.find("(") >=0:
        n = n.replace("(","")
    if n.find(")") >=0:
        n = n.replace(")", "")
    if n.find(" ") >=0:
        n = n.replace(" ", "_")

    all_enz.add(n.lower())


class SearchFormValidator(schema.Schema):
    query = compound.All(validators.PlainText(),
                         validators.Regex(r'^[ \n]*[A-Za-z0-9\n]*[ \n]*$'),
                         validators.String(min=6, strip = True))

    enzyme = validators.OneOf(all_enz)

    min_l = validators.Int(min = 0)
    max_l = validators.Int(min = 0)
    min_w = validators.Int(min = 0)
    max_w = validators.Int(min = 0)

    misses = validators.Int(min = 0)

    sort_method   = validators.OneOf(["Position in sequence","Weight"])

    # This function was given by Claude
    # used to make sure that the minimums are less than the maximums
    def validate_python(self, value_dict, state):
        if value_dict.get('min_l') and value_dict.get('max_l') and int(value_dict['min_l']) >= int(value_dict['max_l']):
            raise Invalid('Minimum length must be less than maximum length', value_dict, state)
        if value_dict.get('min_w') and value_dict.get('max_w') and int(value_dict['min_w']) >= int(value_dict['max_w']):
            raise Invalid('Minimum weight must be less than maximum weight', value_dict, state)


class SearchForm(twf.Form):
    class child(twf.TableLayout):
        query  = twf.TextArea(label = "Sequence", rows = 3, cols = 50)

        enzyme = twf.SingleSelectField(lable = "Enzyme",
                                       options = [x.name for x in all_enzymes[1:]],
                                       prompt_text = "Trypsin")
        
        min_l = twf.TextField(label = "Min Fragment length")
        max_l = twf.TextField(label = "Max Fragment length")
        min_w = twf.TextField(label = "Min Fragment Weight (Da)")
        max_w = twf.TextField(label = "Max Fragment Weight (Da)")

        misses = twf.SingleSelectField(label = "Number of missed cleavages allowed",
                                       options = ["1","2","3","4","5"],
                                       prompt_text = "0")

        sort_method = twf.SingleSelectField(label = "Display results by",
                                    options = ["Weight"],
                                    prompt_text = "Position in sequence")


        # for presentation only...
        css_class = 'table'
        attrs = {'style': 'width: 600px;'}
    
    action = '/digest'
    submit = twf.SubmitButton(value="Digest")
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
        enzyme_list = all_enzymes # list in the Enzymes module
        return dict(page = 'enzymes',
                    enzyme_list = enzyme_list)

    @expose('digest.templates.about')
    def about(self):
        """Project Information Page"""
        return dict(page = 'about')
    
    ### Adding digest functionality
    @expose("digest.templates.digest")
    @expose("digest.templates.digestxml",
            content_type = "text/xml")
    
    @validate(SearchForm, error_handler = index)
    def digest(self, query, enzyme = None, min_l = None, max_l = None,
               min_w = None, max_w = None, misses = None, sort_method = None):

        ### FORMATTING PARAMATERS
        # sequence
        if any(c in "0123456789" for c in query): #Checking to see if it has a UniProt ID
            seq = Seq(UniProt_acc = query)
        else:
            print("Manual setting sequence")
            seq = Seq(seq = query, name = "User Provided Protein")
        
        # setting enzyme to enzyme class
        if enzyme!= None:
            enzyme = enzyme.lower()
        
        if enzyme == None or enzyme == "trypsin":
            enzyme = Trypsin()
        elif enzyme in ["argc","arg_c","arg c"]:
            enzyme = ArgC()
        elif enzyme in ["aspn", "asp_n", "asp n"]:
            enzyme = AspN()
        elif enzyme in ["lysn", "lys_n", "lys n"]:
            enzyme = Lys_n()
        elif enzyme in ["lysc", "lys_c", "lys c"]:
            enzyme = Lys_c()
        elif enzyme == "cnbr":
            enzyme = CNBr()
        elif enzyme == "Proteinase K".lower() or enzyme == "proteinase_k":
            enzyme = Proteinase_K()
        elif enzyme == "Pepsin (pH 1.3)".lower() or enzyme == "pepsin_ph_1_3":
            enzyme = Pepsin_1_3()
        elif enzyme == "Pepsin (pH > 2)".lower() or enzyme == "pepsin_ph_gt2":
            enzyme = Pepsin_gt2()
        elif enzyme == "Thermolysin".lower():
            enzyme = Thermolysin()
        else:
            raise Exception("Invalid enzyme given, Please input a correct enyzme")

        # Other paramaters. Setting defaults
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

        if sort_method == None:
            # print("storing default value")
            sort_method = "Position"

        #### running the program
        enzyme_digest(seq, enzyme, min_l, max_l, min_w, max_w, misses)
        seq.sort_fragments(sort_method)
        print("sort method:",sort_method)

        # returning results to the webpage
        return dict(seq=seq,
                    enzyme = enzyme,
                    v_frags = seq.valid_fragments,
                    params = [min_l,max_l,min_w,max_w,misses])
    
