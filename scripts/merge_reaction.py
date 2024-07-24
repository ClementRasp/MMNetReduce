#!/usr/bin/env python3
# coding: utf-8 
"""
This file contains functions for merging reactions that are in the same direction or in reverse directions.
The purpose is to create reversible reactions to facilitate michaelis-menten analysis.
"""
import libsbml, numpy, re
#------------------------------------------------------------------------------------------------------------------------
numpy.set_printoptions(threshold=numpy.inf) # Display all the elements of a numpy array
#------------------------------------------------------------------------------------------------------------------------
#-----                                                  Constants                                                   -----
#------------------------------------------------------------------------------------------------------------------------
BOLD = "\033[1m"
RESET = "\033[0m"
MAGENTA = "\033[35m"
#------------------------------------------------------------------------------------------------------------------------
def get_stoichiometric_matrix(model: libsbml.Model, verbose=0) -> numpy.ndarray:
    """Get the stoichiometric matrix of the model.
    
    Parameter
    ---------
    model : libsbml.Model
        The model from which the stoichiometric matrix is extracted.
        
    Return
    ------
    numpy.ndarray :
        The stoichiometric matrix of the model.
    """
    reactions = model.getListOfReactions()
    metabolites = model.getListOfSpecies()
    #---
    metabolite_indices = {met.getId(): i for i, met in enumerate(metabolites)}
    num_metabolites = len(metabolites)
    num_reactions = len(reactions)
    #---
    stoichiometric_matrix = numpy.zeros((num_metabolites, num_reactions))
    for i, reaction in enumerate(reactions):
        for reactant in reaction.getListOfReactants():
            metabolite_id = reactant.getSpecies()
            coeff = reactant.getStoichiometry()
            stoichiometric_matrix[metabolite_indices[metabolite_id], i] -= coeff
        #---
        for product in reaction.getListOfProducts():
            metabolite_id = product.getSpecies()
            coeff = product.getStoichiometry()
            stoichiometric_matrix[metabolite_indices[metabolite_id], i] += coeff
    #---        
    if verbose:
        print(f"{BOLD}Stoichiometric matrix:{RESET}\n{stoichiometric_matrix}\n")
    #---
    return stoichiometric_matrix
#------------------------------------------------------------------------------------------------------------------------
def new_formula(model: libsbml.Model, klaw_1: str, klaw_2: str, op: str) -> str:
    """Create a new formula by merging the formulas of two kinetic laws.
    This function rename parameters that are common to both formulas to avoid conflicts.
    
    Parameters
    ----------
    model : libsbml.Model
        The SBML model object.
        
    klaw_1 : str
        The formula of the kinteic law of the first reaction.
        
    klaw_2 : str
        The formula of the kinteic law of the second reaction.
        
    op : str
        The operator to be used for merging the formulas ('+' or '-').
        
    Return
    ------
    str :
        The new formula.
    """
    e1 = klaw_1.getFormula()
    e2 = klaw_2.getFormula()
    #---
    # Get the list of parameters in the model and reactions
    lst_of_global_param = [p.getId() for p in model.getListOfParameters()] 
    lst_of_param_r1 = [p.getId() for p in klaw_1.getListOfParameters()]
    lst_of_param_r2 = [p.getId() for p in klaw_2.getListOfParameters()]
    #---
    # List the parameters that are common to both reactions to be renamed
    common_param = list(set(lst_of_param_r1) & set(lst_of_param_r2))
    all_param = list(set(lst_of_param_r1) | set(lst_of_param_r2)) + lst_of_global_param
    #---
    for cp in common_param:
        new_name = f"{cp}_"
        while new_name in all_param: # Add an underscore to the name if it already exists
            new_name += "_"
        #---
        # Change the name of the parameter in the kinetic law of the second reaction
        p = klaw_2.getParameter(cp)
        p.setId(new_name)
        #---
        # Add the new name to the list of all parameters
        all_param.append(new_name)
        #---
        # Replace the old name by the new name in the formula of the second reaction
        e2 = re.sub(rf'\b{cp}\b', new_name, e2)
    #---
    if op == '+': # Merge the formulas
        return e1 + " + " + e2
    else:
        return e1 + " - (" + e2 + ")"
#------------------------------------------------------------------------------------------------------------------------
def is_reaction_in_a_rule(model: libsbml.Model, reaction: libsbml.Reaction) -> bool:
    """Check if the reaction is in a rule.
    
    Parameters
    ----------
    model : libsbml.Model
        The SBML model object.
        
    reaction : libsbml.Reaction
        The reaction object.
        
    Return
    ------
    bool :
        True if the reaction is in a rule, False otherwise.
    """
    for rule in model.getListOfRules():
        if re.search(r'\b' + re.escape(reaction.getId()) + r'\b', rule.getFormula()):
            return True
    return False
#
#------------------------------------------------------------------------------------------------------------------------
#-----                                       Functions for merging reactions                                        -----
#------------------------------------------------------------------------------------------------------------------------
# 
def merge_reaction(model: libsbml.Model, verbose = 0) -> None:
    """Modify the model by merging reactions that are in the same direction or in reverse directions if they have the same
    reactants and products.
    
    Parameter
    ---------
    model : libsbml.Model
        The model to be modified.
    """
    # Get all reactions and the stoichiometric matrix
    all_reactions = model.getListOfReactions() 
    matrix = get_stoichiometric_matrix(model)
    #---
    # Transpose the matrix to iterate over reactions instead of metabolites
    matrix = matrix.T 
    #---
    reaction_merged = []
    suppressed = set()
    #---
    for i in range(len(matrix)): # Iterate over reactions
        if i in suppressed:
            # Skip if the reaction has already been merged
            continue
        #---
        if all_reactions[i].getKineticLaw() == None or is_reaction_in_a_rule(model, all_reactions[i]):
            # If the reaction has no kinetic law or is in a rule, skip it.
            continue
        #---
        for j in range(i+1, len(matrix)): # Iterate over other reactions
            if j in suppressed:
                # Skip if the reaction has already been merged
                continue
            if all_reactions[j].getKineticLaw() == None or is_reaction_in_a_rule(model, all_reactions[j]):
                # If the reaction has no kinetic law or is in a rule, skip it.
                continue
            #---
            species_i = [reactant.getSpecies() for reactant in all_reactions[i].getListOfReactants()] + [product.getSpecies() for product in all_reactions[i].getListOfProducts()]
            species_j = [reactant.getSpecies() for reactant in all_reactions[j].getListOfReactants()] + [product.getSpecies() for product in all_reactions[j].getListOfProducts()]
            #---
            if numpy.array_equal(matrix[i], matrix[j]) and sorted(species_i) == sorted(species_j): 
                # If reactions are the same, merge them
                op = '+'
            elif numpy.array_equal(matrix[i], -matrix[j]) and sorted(species_i) == sorted(species_j): 
                # If reactions are the reverse of each other, merge them
                op = '-'
                all_reactions[i].setReversible(True)
            else:
                continue
            #---
            # Get the kinetic law of the reaction and create a new formula from the formulas of the two reactions
            k_law = all_reactions[i].getKineticLaw()
            expression = new_formula(model, k_law, all_reactions[j].getKineticLaw(), op)
            #---
            t = k_law.setFormula(f"{expression}")
            if t != libsbml.LIBSBML_OPERATION_SUCCESS:
                raise ValueError(f"Error while setting the formula for the reaction {all_reactions[i].getId()}.")
            #---
            # Add the parameters of the second reaction to the first reaction
            for param in all_reactions[j].getKineticLaw().getListOfParameters():
                p = k_law.createParameter()
                p.setId(param.getId())
                p.setValue(param.getValue())
                p.setUnits(param.getUnits())
                p.setMetaId(param.getMetaId())
            #---
            # Add the modifiers of the second reaction to the first reaction
            for modifier in all_reactions[j].getListOfModifiers():
                m = all_reactions[i].createModifier()
                m.setSpecies(modifier.getSpecies())
                m.setId(modifier.getId())
                m.setName(modifier.getName())
            #---
            # Name the new reaction by merging the names of the two reactions
            all_reactions[i].setName(all_reactions[i].getName() + '//' + all_reactions[j].getName())
            all_reactions[i].setId(all_reactions[i].getId() + '__' + all_reactions[j].getId())
            #---
            # Add the index of the second reaction to the list of suppressed reactions
            suppressed.add(j) 
            #---
            if verbose:
                reaction_merged.append((op, all_reactions[i].getId(), all_reactions[j].getId()))    
    #---
    # Remove the suppressed reactions
    for i in sorted(suppressed, reverse=True):
        model.removeReaction(i)
    #---
    if verbose:
        if reaction_merged:
            print(f"{BOLD}Reaction modification:{RESET}\n")
            for op, r1, r2 in reaction_merged:
                if op == '+':
                    print(f"Reactions {MAGENTA}{r1}{RESET} and {MAGENTA}{r2}{RESET} "
                    f"are in the same direction, thus they have been merged into a single reaction {MAGENTA}{r1 + '__' + r2}{RESET}, while maintaining their original directions.\n")
                else:
                    print(f"Reactions {MAGENTA}{r1}{RESET} and {MAGENTA}{r2}{RESET} "
                    f"are are in reverse directions, so they have been amalgamated into a single reversible reaction {MAGENTA}{r1 + '__' + r2}{RESET}.\n")
            print()