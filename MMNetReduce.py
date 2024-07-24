#!/usr/bin/env python3
# coding: utf-8 
"""
This Python program can be used to  ..... 
Upon execution, it  initiates a dialog prompting the user to select the desired SBML file for reading. 
The output comprises the SBML file containg the reduced Michaelis-Menten model.
"""
#
#------------------------------------------------------------------------------------------------------------------------
#-----                                                   Imports                                                    -----
#------------------------------------------------------------------------------------------------------------------------
# 
import time, warnings, re, glob
from useful import *
from merge_reaction import merge_reaction
#
#------------------------------------------------------------------------------------------------------------------------
#-----                                           Michaelis-Menten Schemes                                           -----
#------------------------------------------------------------------------------------------------------------------------
#
def get_model_functions (model: libsbml.Model, verbose = 0) -> tuple[list[str], list[list[str]], list[str]]:
    """ Get the list of functions, their arguments and expressions in the model.
    
    Parameter
    ---------
    model : libsbml.Model
        The model object.
        
    Returns
    -------
    list[str]
        The list of function id.
        
    list[list[str]]
        The list of function arguments.
        
    list[str]
        The list of function expressions.
    """
    #---
    def extract_arguments(Inside):
        Arguments = []
        #---
        #Function to check if a character is a mathematical symbol
        def is_math_symbol(char):
            return char in ('+', '-', '*', '/', '^', '(', ')')
        #---
        #Loop until a mathematical symbol is encountered
        while Inside:
            #Find the index of the first comma
            comma_index = Inside.find(',')
            #---
            #If there is no comma or the comma is at the start, break the loop
            if comma_index == -1 or comma_index == 0:
                break
            #---
            #Extract the part before the first comma
            part = Inside[:comma_index].strip()
            #If the part contains a mathematical symbol, break the loop
            if any(is_math_symbol(c) for c in part):
                break
            #---
            #Add the part to Arguments and remove it from Inside along with the comma
            Arguments.append(part)
            Inside = Inside[comma_index+1:].strip()
        return Arguments, Inside
    #---
    #Get the list of functions in the model
    function_id_list = []
    function_list = []
    #---
    for i in range(model.getNumFunctionDefinitions()):        
        function = model.getFunctionDefinition(i)
        function_id = function.getId()
        math_ast = function.getMath()
        math_string = libsbml.formulaToString(math_ast)
        function_id_list.append(function_id)
        function_list.append(math_string)
        if verbose > 1:
            print(f"{BOLD}Function{RESET} {i+1}:")
            print("ID:", function_id)
            print("Expression:", math_string)
    #---
    arguments = []
    expressions = []
    for i in range(len(function_id_list)):
        #---
        inside_parentheses_i = re.search(r'lambda\((.*)\)', function_list[i]).group(1)
        #print("Inside parantheses:", inside_parentheses_i)
        #---
        arguments_i, expression_i = extract_arguments(inside_parentheses_i)
        #---
        arguments.append(arguments_i)
        expressions.append(expression_i)
        #---
    return function_id_list, arguments, expressions
#------------------------------------------------------------------------------------------------------------------------
def get_rate(model: libsbml.Model, reaction: libsbml.Reaction) -> str:
    """ Get the rate of a reaction.
    Change functions by their expressions in the rate.
    
    Parameters
    ----------
    model : libsbml.Model
        The model object.
        
    reaction : libsbml.Reaction
        The reaction object.
        
    Return
    ------
    str
        The rate of the reaction.
    """
    function_id, arguments, expressions = get_model_functions(model,0)
    #---
    rate = reaction.getKineticLaw().getFormula()
    f_j =  [[func, re.search(r'\b{}\b'.format(re.escape(func)) + r'\((.*)\)', rate).group(1)] for func in function_id if re.search(r'\b{}\b'.format(re.escape(func)) + r'\((.*)\)', rate)]
    #---
    for i in range(len(f_j)):
        elem_ji = f_j[i]
        #---
        fun_id = elem_ji[0]
        new_args = elem_ji[1].replace(" ","").split(',') # AJOUT : remove spaces
        new_name = fun_id+ '(' + ', '.join(new_args) + ')'
        #---
        old_arg = arguments[function_id.index(fun_id)]
        old_expres = expressions[function_id.index(fun_id)]
        #---
        new_expess = old_expres
        for i in range(len(old_arg)):
            new_expess = re.sub(r'\b{}\b'.format(re.escape(old_arg[i])), new_args[i], new_expess) # Changement de la fonction replace
        #---
        rate = rate.replace(new_name, new_expess)
    #---
    return rate
#------------------------------------------------------------------------------------------------------------------------
def get_parameters(model: libsbml.Model, reaction_1: libsbml.Reaction, reaction_2: libsbml.Reaction) -> tuple[libsbml.Parameter, libsbml.Parameter, libsbml.Parameter]:
    """ Get the parameters k+, k- and kcat of the Michaelis-Menten scheme.
    
    Parameters
    ----------
    model : libsbml.Model
        The model object.
        
    reaction_1 : libsbml.Reaction
        The reversible reaction object.
        
    reaction_2 : libsbml.Reaction
        The irreversible reaction object.
        
    Return
    ------
    tuple[libsbml.Parameter, libsbml.Parameter, libsbml.Parameter]
        The parameters k+, k- and kcat of the Michaelis-Menten scheme.
    """
    # List all parameters that can be used in the rate of the reaction
    list_of_parameters = model.getListOfParameters()
    list_of_reaction_1_parameters = [p for p in reaction_1.getKineticLaw().getListOfParameters()]
    list_of_reaction_2_parameters = [p for p in reaction_2.getKineticLaw().getListOfParameters()]
    #---
    if len(list_of_reaction_2_parameters) == 1:
        # If there's only one local parameter in the reaction, it's the kcat parameter
        k_cat = list_of_reaction_2_parameters[0]
    elif len(list_of_reaction_2_parameters) > 1 :
        # If there's more than one local parameter in the reaction, it's not possible to determine the kcat parameter
        return None, None, None
    else:
        # If there's no local parameter in the reaction, the kcat parameter is global, we need to find it in the
        # rate of the reaction.
        k_law_2 = get_rate(model, reaction_2)
        parameters_founded = []
        for p in list_of_parameters:
            if re.search(r'\b' + re.escape(p.getId()) + r'\b', k_law_2):
                parameters_founded.append(p)
        #---
        if len(parameters_founded) == 1:
            k_cat = parameters_founded[0]
        else:
            return None, None, None
    #---
    # For the first reaction, we need to determine the k+ and k- parameters,
    # we can find them by spliting the rate of the reaction.
    k_law_1 = get_rate(model, reaction_1)
    split = k_law_1.split('-')
    if len(split) != 2:
        # If there's more or less than two parts after the split, then it's not possible to determine
        # the k+ and k- parameters. 
        return None, None, None
    k_law_1_positive, k_law_1_negative = split[0], split[1]
    #---
    if len(list_of_reaction_1_parameters) == 2:
        # If there's only two parameters in the reaction, we can determine the k+ and k- parameters.
        # If the first parameter is in the positive part of the rate, then it's the k+ parameter,
        # otherwise it's the k- parameter.
        if re.search(r'\b' + re.escape(list_of_reaction_1_parameters[0].getId()) + r'\b', k_law_1_positive):
            k_plus = list_of_reaction_1_parameters[0]
            k_minus = list_of_reaction_1_parameters[1]
        else:
            k_plus = list_of_reaction_1_parameters[1]
            k_minus = list_of_reaction_1_parameters[0]
    elif len(list_of_reaction_1_parameters) > 2:
        # If there's more than two parameters in the reaction, it's not possible to determine the k+ and k- parameters.
        return None, None, None
    else:
        # Else, we need to find the k+ and k- parameters in the rate of the reaction.
        parameters_founded = []
        for p in list_of_parameters:
            if re.search(r'\b' + re.escape(p.getId()) + r'\b', k_law_1_positive):
                parameters_founded.append(p)
        #---
        if len(parameters_founded) == 1:
            k_plus = parameters_founded[0]
        else:
            return None, None, None
        #---
        parameters_founded = []
        for p in list_of_parameters:
            if re.search(r'\b' + re.escape(p.getId()) + r'\b', k_law_1_negative):
                parameters_founded.append(p)
        #---
        if len(parameters_founded) == 1:
            k_minus = parameters_founded[0]
        else:
            return None, None, None
    #---
    return k_plus, k_minus, k_cat
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
def get_Michaelis_Menten_scheme(model: libsbml.Model, verbose = 0) -> list[dict]:
    """ Get the Michaelis-Menten schemes from the model.
    
    Parameters
    ----------
    model : libsbml.Model
        The model object.
    
    Return
    ------
    list[dict]
        The list of Michaelis-Menten schemes.
    """
    michaelis_menten_scheme = []
    id = 0
    #---
    list_of_reaction = model.getListOfReactions()
    buffered_species_id = [species.getId() for species in model.getListOfSpecies() if species.getConstant()]
    #---
    for reaction_1 in list_of_reaction:
        if not reaction_1.getReversible() or reaction_1.getKineticLaw() == None:
            # Skip if the first reaction is not reversible or if the reaction doesn't have rate. 
            continue
        for reaction_2 in list_of_reaction:
            if reaction_2.getReversible() or reaction_2.getKineticLaw() == None:
                # Skip if the second reaction is reversible to keep only this type of reaction A <---> B, C ---> D
                # or if the reaction doesn't have rate. 
                continue
            #---
            c_1 = set([species for species in [reactant.getSpecies() for reactant in reaction_1.getListOfReactants()] if species not in buffered_species_id])
            c_2 = set([species for species in [product.getSpecies() for product in reaction_1.getListOfProducts()] if species not in buffered_species_id])
            d_2 = set([species for species in [reactant.getSpecies() for reactant in reaction_2.getListOfReactants()] if species not in buffered_species_id])
            #---
            if c_2 == d_2:
                d_1 = c_1
            elif c_1 == d_2:
                d_1 = c_2
            else:
                # Skip if the reaction 1 reactant or product is not the same as the reaction 2 reactant
                # to keep only this type of reaction B = C : A <---> C ---> D or A = C : B <---> C ---> D
                continue
            #---
            d_3 = set([species for species in [product.getSpecies() for product in reaction_2.getListOfProducts()] if species not in buffered_species_id])
            #---
            if len(d_1) != 2 or len(d_2) != 1 or len(d_3) < 1:
                # Skip if the number of reactants or products is not correct
                continue
            #---
            if len(d_1 & d_3) != 1:
                # Skip if there's more than one species in common between d_1 and d_3 (should have one enzyme)
                continue
            #---
            k_plus, k_minus, k_cat = get_parameters(model, reaction_1, reaction_2)
            if k_plus is None:
                # Skip if it's not possible to determine the parameters k+, k- and kcat.
                continue
            #---
            # Create a dictionary to store the Michaelis-Menten scheme
            scheme = dict() 
            enzyme = d_1 & d_3
            product = d_3 - enzyme
            substrate = d_1 - enzyme
            scheme['reaction_1'] = reaction_1
            scheme['reaction_2'] = reaction_2
            scheme['enzyme'] = enzyme.pop()
            scheme['substrate'] = list(substrate)[0]
            scheme['product'] = list(product)
            scheme['intermediate'] = d_2.pop()
            scheme['id'] = id = id + 1
            scheme['k_plus'] = k_plus
            scheme['k_minus'] = k_minus
            scheme['k_cat'] = k_cat
            michaelis_menten_scheme.append(scheme)
    #---
    if verbose:
        for scheme in michaelis_menten_scheme:
            print(f"{BOLD}Michaelis-Menten scheme {scheme['id']}:{RESET}")
            print(f"\t{BOLD}Reversible reaction ID:{RESET} {scheme['reaction_1'].getId()}")
            print(f"\t{BOLD}Irreversible reaction ID:{RESET} {scheme['reaction_2'].getId()}")
            full_scheme = f"{BLUE}" + scheme['substrate'] + f"{RESET} + {GREEN}{scheme['enzyme']}{RESET} <---> {ORANGE}{scheme['intermediate']}{RESET} ---> {GREEN}{scheme['enzyme']}{RESET} + {RED}" + " + ".join(scheme['product']) + f"{RESET}"
            print(f"\t{BOLD}Scheme:{RESET} {full_scheme}")
            print(f"\t{BOLD}Enzyme:{RESET} {GREEN}{scheme['enzyme']}{RESET}")
            print(f"\t{BOLD}Substrates:{RESET} {BLUE}{scheme['substrate']}{RESET}")
            print(f"\t{BOLD}Products:{RESET} {RED}{', '.join(scheme['product'])}{RESET}")
            print(f"\t{BOLD}Intermediate:{RESET} {ORANGE}{scheme['intermediate']}{RESET}")
            print()
        print()
    #---     
    return michaelis_menten_scheme
#
#------------------------------------------------------------------------------------------------------------------------
#-----                                            Michaelis-Menten Pools                                            -----
#------------------------------------------------------------------------------------------------------------------------
#
def get_units_for_parameter(model: libsbml.Model, param: libsbml.Parameter, reaction: libsbml.Reaction) -> list[libsbml.Unit]:
    """Return the list of units for a parameter in a reaction.
    
    Parameters
    ----------
    model : libsbml.Model
        The SBML model object.
        
    param : libsbml.Parameter
        The parameter object.
    
    reaction : libsbml.Reaction
        The reaction object.
    
    Return
    ------
    list
        The list of units for the parameter in the reaction.
    """
    unit_def = model.getUnitDefinition(param.getUnits())
    if not unit_def:
        warnings.warn(f"No unit definition was found for the parameter {param.getId()} in the reaction {reaction.getId()}.", UserWarning)
        return None
    #---
    l = []
    for unit in unit_def.getListOfUnits():
        l.append(libsbml.Unit(unit))
    #---
    return l
#------------------------------------------------------------------------------------------------------------------------
def find_unit_name(list_of_unit: list[libsbml.Unit], default: str, Km_value: float) -> tuple[str, float]:
    """
    Find the corresponding name of the unit from the list of units to have an explicit name.
    
    Parameters
    ----------
    list_of_unit : list
        The list of units.
        
    default : str
        The default name of the unit.
        
    Km_value : float
        The value of the Km parameter.
        
    Return
    ------
    tuple
        The name of the unit and the Km value.   
    """
    name = ""
    #---
    list_of_unit = sorted(list_of_unit, key=lambda unit: (unit.getExponent() < 0, abs(unit.getExponent())))
    #---
    for unit in list_of_unit:
        if unit.getKind() == libsbml.UNIT_KIND_DIMENSIONLESS:
            return default, Km_value
        #---
        if unit.getExponent() < 0:
            name += "per_"
        if abs(unit.getExponent()) == 2:
                name += "square_"
        #---    
        if unit.getScale() >= -3 and unit.getScale() < 0:
            name += "m"
        elif unit.getScale() >= -6 and unit.getScale() < -3:
            name += "u"
        elif unit.getScale() >= -9 and unit.getScale() < -6:
            name += "n"
        elif unit.getScale() >= -12 and unit.getScale() < -9:
            name += "p"
        elif unit.getScale() >= 3 and unit.getScale() < 6:
            name += "k"
        elif unit.getScale() >= 6 and unit.getScale() < 9:
            name += "M" 
        #---
        diff = unit.getScale() % 3
        if diff != 0:
            Km_value *= 10**diff
        #---
        if unit.getKind() == libsbml.UNIT_KIND_MOLE:
            name += "mol"
        elif unit.getKind() == libsbml.UNIT_KIND_SECOND:
            name += "second"        
        elif unit.getKind() == libsbml.UNIT_KIND_METRE:
            name += "metre"
        elif unit.getKind() == libsbml.UNIT_KIND_METER:
            name += "meter"
        elif unit.getKind() == libsbml.UNIT_KIND_GRAM:
            name += "gram"
        elif unit.getKind() == libsbml.UNIT_KIND_KILOGRAM:
            name += "kilogram"
        elif unit.getKind() == libsbml.UNIT_KIND_LITRE:
            name += "litre"
        elif unit.getKind() == libsbml.UNIT_KIND_LITER:
            name += "liter"
        else:
            return default
        #---
        name += "_"
    #---
    return name[:-1], Km_value
#------------------------------------------------------------------------------------------------------------------------
def unitDefinition_from_units(model: libsbml.Model, unit_id: str, units_list: list[libsbml.Unit]) -> libsbml.UnitDefinition:
    """
    Get the unit definition from the list of units.
    If the unit definition already exists, return it, otherwise create it.
    
    Parameters
    ----------
    model : libsbml.Model
        The SBML model object.
        
    unit_id : str
        The id of the unit definition.
        
    units_list : list
        The list of units.
        
    Return
    ------
    libsbml.UnitDefinition
        The unit definition.
    """
    def is_unit_definition(model: libsbml.Model, list_of_unit: list[libsbml.Unit]) -> libsbml.UnitDefinition:
        for unit in model.getListOfUnitDefinitions():
            if sorted([u.getKind() for u in unit.getListOfUnits()]) == sorted([u.getKind() for u in list_of_unit]):
                return unit
        return None 
    #---
    # Check if the unit definition already exists
    unit_def = is_unit_definition(model, units_list)
    # If the unit definition doesn't exist, create it
    if unit_def == None:
        unit_def = model.createUnitDefinition()
        unit_def.setId(unit_id)
        #---
        for unit in units_list:
            u = unit_def.createUnit()
            u.setKind(unit.getKind())
            u.setExponent(unit.getExponent())
            u.setScale(unit.getScale())
            u.setMultiplier(unit.getMultiplier())
    #---
    return unit_def
#------------------------------------------------------------------------------------------------------------------------
def get_Km_Kcat(model: libsbml.Model, reaction_reversible: libsbml.Reaction, reaction_irreversible: libsbml.Reaction,
                Kplus: libsbml.Parameter, Kminus: libsbml.Parameter, list_of_Kcat: list[libsbml.Parameter], i: int) -> tuple[float, libsbml.Unit, libsbml.Parameter]:
    """
    Get the Km and Kcat values and units from the reversible and irreversible reactions.
    
    Parameters
    ----------
    model : libsbml.Model
        The SBML model object.
        
    reaction_reversible : libsbml.Reaction
        The reversible reaction object.
    
    reaction_irreversible : libsbml.Reaction
        The irreversible reaction object.
    
    Return
    ------
    tuple
        The Km value, the Km unit, and the Kcat parameter.    
        
    Warnings
    --------
    If one of the parameters has no unit definition, a warning is issued.\\
    If the unit of Kcat and Kminus are not the same, a warning is issued.\\
    If the unit of Kcat parameter is not in the form '1/second', a warning is issued.
    """
    list_Kplus_unit = get_units_for_parameter(model, Kplus, reaction_reversible)
    list_Kminus_unit = get_units_for_parameter(model, Kminus, reaction_reversible)
    list_of_Kcat_unit = [get_units_for_parameter(model, Kcat, reaction_irreversible) for Kcat in list_of_Kcat]
    #---
    if list_Kplus_unit == None or list_Kminus_unit == None or None in list_of_Kcat_unit:
        # If one of the parameters has no unit definition, Km will be dimensionless.
        warnings.warn(f"One of the parameters has no unit definition.\nUnits will be ignore.", UserWarning)
        #---
        list_of_Kcat[i].setUnits("dimensionless")
        if Kplus.getValue() == 0:
            Km_value = 0
        else:
            Km_value = (sum([Kcat.getValue() for Kcat in list_of_Kcat]) + Kminus.getValue()) / Kplus.getValue()
        #---
        Km_unit = "dimensionless"
        #---
        return Km_value, Km_unit, list_of_Kcat[i]
    #---
    if all(element == list_of_Kcat[0] for element in list_of_Kcat):
        # If one Kcat doesn't have the same value as the others, a warning is issued.
        warnings.warn(f"Kcat parameters for the complex {[s.getSpecies() for s in reaction_irreversible.getListOfReactants()]} have different values.", UserWarning)
    #---
    if len(list_of_Kcat_unit[i]) != 1:
        # If there's more than one unit for the Kcat parameter, we'll be searching for the unit of time.
        warnings.warn(f"Kcat parameter {list_of_Kcat[i].getId()} in the reaction {reaction_irreversible.getId()} has more than one unit.", UserWarning)
        Kcat_unit = list_of_Kcat_unit[i][0]
        #---
        for u in list_of_Kcat_unit[i]:
            if u.getKind() == libsbml.UNIT_KIND_SECOND:
                # If the unit of time is found, we'll chose it as the unit of all Kcat parameters.
                Kcat_unit = u
    else:
        # Get the unit of Kcat, we'll chose the i unit as the unit of all Kcat parameters.
        Kcat_unit = list_of_Kcat_unit[i][0]
    #---
    if len(list_Kminus_unit) != 1:
        # Same things, Kminus need to be a unit of time. So we'll be searching for the unit of time.
        warnings.warn(f"Kminus parameter {Kminus.getId()} in the reaction {reaction_reversible.getId()} has more than one unit.", UserWarning)
        Kminus_unit = list_Kminus_unit[0]
        for u in list_Kminus_unit:
            if u.getKind() == libsbml.UNIT_KIND_SECOND:
                Kminus_unit = u
    else:
        Kminus_unit = list_Kminus_unit[0]
    #---
    if (not libsbml.Unit.areEquivalent(Kcat_unit, Kminus_unit)):
        # Kcat and Kminus should be the same unit, a unit of time.
        warnings.warn(f"The unit of Kcat and Kminus are not the same.", UserWarning)
    #---
    if Kcat_unit.getKind() != libsbml.UNIT_KIND_SECOND or Kcat_unit.getExponent() != -1:
        # Verify if the unit of Kcat is in the form '1/second'.
        warnings.warn(f"The unit of Kcat parameter {list_of_Kcat[i].getId()} in the reaction {reaction_irreversible.getId()} is not in the form '1/second'.", UserWarning) 
    #---
    # If there's a multiplier for the unit of Kcat, we'll multiply all Kcat parameters by this multiplier.
    # to get a metric unit.
    for j in range(len(list_of_Kcat)):
        if list_of_Kcat[j].getMultiplier() != 1:
            list_of_Kcat[j].setValue(list_of_Kcat[j].getValue() * list_of_Kcat[j].getMultiplier())
            list_of_Kcat[j].setMultiplier(1)
    #---
    # If there's a multiplier for the unit of Kminus, we 'll multiply Kminus by this multiplier.
    if Kminus_unit.getMultiplier() != 1:
        Kminus.setValue(Kminus.getValue() * Kminus_unit.getMultiplier())
        Kminus_unit.setMultiplier(1)
    #--- 
    # Changing the scale of Kcat and Kminus to the same scale.
    if Kcat_unit.getScale() > Kminus_unit.getScale():
        list_of_Kcat[i].setValue(list_of_Kcat[i].getValue() * 10**(Kcat_unit.getScale() - Kminus_unit.getScale()))
        Kcat_unit.setScale(Kminus_unit.getScale())
    elif Kcat_unit.getScale() < Kminus_unit.getScale():
        Kminus.setValue(Kminus.getValue() * 10**(Kminus_unit.getScale() - Kcat_unit.getScale()))
        Kminus_unit.setScale(Kcat_unit.getScale())
    #---
    # For the unit of Kplus, we searching for the unit of time. If we find it, we'll change the scale of Kplus 
    # to the same scale as Kcat.
    for u in list_Kplus_unit:
        if u.getKind() == libsbml.UNIT_KIND_SECOND:
            if u.getMultiplier() != 1:
                Kplus.setValue(Kplus.getValue() * u.getMultiplier())
                u.setMultiplier(1)
            if u.getScale() > Kcat_unit.getScale():
                Kplus.setValue(Kplus.getValue() * 10**(u.getScale() - Kcat_unit.getScale()))
                u.setScale(Kcat_unit.getScale())
            elif u.getScale() < Kcat_unit.getScale():
                Kplus.setValue(Kplus.getValue() * 10**(Kcat_unit.getScale() - u.getScale()))
                u.setScale(Kcat_unit.getScale())
            break
    else:
        warnings.warn(f"Kplus parameter {Kplus.getId()} in the reaction {reaction_reversible.getId()} has no unit of time.", UserWarning)
    #---
    # Math part to get the Km value.
    Kplus_value = Kplus.getValue()
    Kminus_value = Kminus.getValue()
    Kcat_value = sum([Kcat.getValue() for Kcat in list_of_Kcat])
    #---
    if Kplus_value == 0:
        Km_value = 0
    else:
        Km_value = (Kcat_value + Kminus_value) / Kplus_value
    #---
    if len(list_Kplus_unit) == 1:
        # If there's only one unit for Kplus, we'll check if it's a unit of time, if yes Km is dimensionless.
        if list_Kplus_unit[0].getKind() == libsbml.UNIT_KIND_SECOND:
            return Km_value, "dimensionless", list_of_Kcat[i]
    #---
    # Create Km unit
    list_Km_unit = []
    for u in list_Kplus_unit:
        if u.getKind() != libsbml.UNIT_KIND_SECOND:
            new_unit = libsbml.Unit(u)
            new_unit.setExponent(-new_unit.getExponent())
            list_Km_unit.append(new_unit)
        else:
            if u.getExponent() != -1:
                warnings.warn(f"Kplus parameter {Kplus.getId()} in the reaction {reaction_reversible.getId()} has not the right unit of time.", UserWarning)
    #---
    if len(list_Km_unit) == 0:
        warnings.warn(f"Kplus parameter {Kplus.getId()} in the reaction {reaction_reversible.getId()} has no unit.", UserWarning) 
        return Km_value, "dimensionless", list_of_Kcat[i]
    #---
    # Find the name of the unit for Km
    Km_name, Km_value =  find_unit_name(list_Km_unit, f"Km_unit_{reaction_reversible.getId()}", Km_value)
    # Search if the unit definition already exists, else create it.
    Km_unit = unitDefinition_from_units(model, Km_name, list_Km_unit)
    Km_name = Km_unit.getId()
    model.addUnitDefinition(Km_unit)
    #--- 
    return Km_value, Km_name, list_of_Kcat[i]
#------------------------------------------------------------------------------------------------------------------------
def schemes_modifications(model: libsbml.Model, list_of_schemes: list[dict]) -> dict:
    """ Remove schemes with complex used in reactions that are not in some schemes.
    Compute also the Km by pooling the Kcat values.
    
    Parameters
    ----------
    model : libsbml.Model
        The SBML model object.
        
    list_of_schemes : list
        The list of Michaelis-Menten schemes.
        
    Return
    ------
    dict
        The dictionary of Michaelis-Menten schemes stored by pool corresponding to the complex of the scheme.
    """
    pool = dict()
    l = []
    #---
    # Create the pool
    for scheme in list_of_schemes:
        if scheme['intermediate'] not in pool:
            pool[scheme['intermediate']] = [scheme]
        else:
            pool[scheme['intermediate']].append(scheme)
    #--- 
    for schemes in pool.values():
        " Listing all the Kcat values for the pool of schemes."
        list_of_k_cat = [kcat['k_cat'] for kcat in schemes]
        #---
        for i, scheme in enumerate(schemes):
            # For each scheme, we compute the Km value and unit.
            k_cat = scheme['k_cat']
            k_plus = scheme['k_plus']
            k_minus = scheme['k_minus']
            reaction_1 = scheme['reaction_1']
            reaction_2 = scheme['reaction_2']
            #---
            k_m_value, k_m_unit, k_cat = get_Km_Kcat(model, reaction_1, reaction_2, k_plus, k_minus, list_of_k_cat, i)
            #---
            scheme['k_m_value'] = k_m_value
            scheme['k_m_unit'] = k_m_unit
            #---
            l.append(scheme)
    #---
    # We remove the schemes with complex used in reactions that are not in some schemes.
    intermediate_to_remove = []
    #---
    for intermediate, schemes in pool.items():
        list_of_reactions = [scheme['reaction_1'] for scheme in schemes]
        list_of_reactions += [scheme['reaction_2'] for scheme in schemes]
        for reaction in model.getListOfReactions():
            # If the complex is in a reaction that isn't in the list of reactions of the pool
            # we remove the complex.
            if reaction not in list_of_reactions:
                species = [r.getSpecies() for r in reaction.getListOfReactants()]
                species += [p.getSpecies() for p in reaction.getListOfProducts()]
                #---
                if intermediate in species and intermediate not in intermediate_to_remove:
                    intermediate_to_remove.append(intermediate)
                    break
        #---
        for rule in model.getListOfRules():
            # If the complex is in a rule that isn't in the list of reactions of the pool
            if re.search(r'\b{}\b'.format(re.escape(intermediate)), rule.getFormula()):
                if intermediate not in intermediate_to_remove:
                    intermediate_to_remove.append(intermediate)
                    break
            #---
            if rule.getVariable() == intermediate:
                if intermediate not in intermediate_to_remove:
                    intermediate_to_remove.append(intermediate)
                    break
    #---
    # Remove schemes that have the same complex but with to different reversible reactions.
    for schemes in pool.values():
        all_reaction = []
        enzymes = []
        for scheme in schemes:
            if scheme['reaction_1'] not in all_reaction:
                all_reaction.append(scheme['reaction_1'])
            if len(all_reaction) > 1:
                if scheme['intermediate'] not in intermediate_to_remove:
                    intermediate_to_remove.append(scheme['intermediate'])
            if scheme['enzyme'] not in enzymes:
                enzymes.append(scheme['enzyme'])
            if len(enzymes) > 1:
                if scheme['intermediate'] not in intermediate_to_remove:
                    intermediate_to_remove.append(scheme['intermediate'])
    #---
    for intermediate in intermediate_to_remove:
        pool.pop(intermediate)
    #---
    return pool
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
def get_Michaelis_Menten_pools(model: libsbml.Model, verbose=0) -> dict:
    """Get the Michaelis-Menten schemes stored by pool corresponding to the enzyme of the scheme.
    
    Parameter
    ---------
    model : libsbml.Model
        The model from which the Michaelis-Menten schemes are extracted.
        
    Return
    ------
    dict
        The dictionary of Michaelis-Menten schemes stored by pool corresponding to the complex of the scheme
        stored by pool corresponding to the enzyme of the scheme.
    {enzyme1: [{complex1: [scheme1, scheme2]}, {complex2: [scheme3]}], enzyme2: [{complex3: [scheme4]}]}
    
    Display exemple
    ---------------
    Michaelis-Menten pool corresponding to the enzyme D.
        Michaelis-Menten scheme 1 for the complex A.
        Michaelis-Menten scheme 2 for the complex A.
        Michaelis-Menten scheme 3 for the complex B.
        
    Michaelis-Menten pool corresponding to the enzyme E.
        Michaelis-Menten scheme 4 for the complex C.
    """
    list_of_scheme = get_Michaelis_Menten_scheme(model, verbose)
    pool_of_intermediate = schemes_modifications(model, list_of_scheme)
    #---
    pools = dict()
    for key, values in pool_of_intermediate.items():
        if values[0]['enzyme'] not in pools:
            pools[values[0]['enzyme']] = [{key: values}]
        else:
            pools[values[0]['enzyme']].append({key: values})
    #---  
    if verbose:
        for enzyme, schemes in pools.items():
            print(f"{BOLD}Michaelis-Menten pool corresponding to the enzyme {GREEN}{enzyme}{RESET}.")
            for scheme in schemes:
                for intermediate, values in scheme.items():
                    for nb in values:
                        print(f"\tMichaelis-Menten scheme {BOLD}{nb['id']}{RESET} for the complex {intermediate}.")

            print()
        print()
    #---
    return pools
#
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
#
def add_intermediates_to_rule(model: libsbml.Model, intermediate: str, formula: str) -> None:
    """
    Add the intermediate species to the assignment rule.
    
    Parameters
    ----------
    model : libsbml.Model
        The SBML model object.
    
    intermediate : str
        The id of the intermediate species.
        
    formula : str
        The formula of the intermediate species.
    """
    rule = model.createAssignmentRule()
    rule.setVariable(intermediate)
    rule.setFormula(formula)
#------------------------------------------------------------------------------------------------------------------------
def delet_reactions(model: libsbml.Model, pools_of_schemes: dict):
    """
    Delete the original reactions from the model.
    
    Parameters
    ----------
    model : libsbml.Model
        The SBML model object.
        
    pools_of_schemes : dict
        The dictionary of Michaelis-Menten schemes stored by pool corresponding to the enzyme of the scheme.
    """
    for list_of_pool in pools_of_schemes.values():
        for pool in list_of_pool:
            for list_of_schemes in pool.values():
                for scheme in list_of_schemes:
                    model.removeReaction(scheme['reaction_2'].getId())
                model.removeReaction(list_of_schemes[0]['reaction_1'].getId())
#------------------------------------------------------------------------------------------------------------------------
def remove_not_used_parameters(model: libsbml.Model) -> None:
    """
    Remove the parameters that are not used in the model.
    
    Parameters
    ----------
    model : libsbml.Model
        The SBML model object.
    """
    not_used_parameters = [p for p in model.getListOfParameters()]
    #---
    # If a parameter is used in a reaction, a rule, an event, an initial assignment, or a constraint,
    # we remove it from the list of not used parameters. Parameters not used will be removed from the model.
    for rule in model.getListOfRules():
        i = 0
        while i < len(not_used_parameters):
            p = not_used_parameters[i]
            if re.search(r'\b{}\b'.format(re.escape(p.getId())), rule.getFormula()):
                not_used_parameters.remove(p)
            elif rule.getVariable() == p.getId():
                not_used_parameters.remove(p)
            else:
                i += 1
    #---
    for reaction in model.getListOfReactions():
        i = 0
        while i < len(not_used_parameters):
            p = not_used_parameters[i]
            if re.search(r'\b{}\b'.format(re.escape(p.getId())), reaction.getKineticLaw().getFormula()):
                not_used_parameters.remove(p)
            else:
                i += 1
    #---      
    for event in model.getListOfEvents():
        i = 0
        while i < len(not_used_parameters):
            p = not_used_parameters[i]
            if re.search(r'\b{}\b'.format(re.escape(p.getId())), libsbml.formulaToString(event.getTrigger().getMath())):
                not_used_parameters.remove(p)
            else:
                i += 1
        #---
        for action in event.getListOfEventAssignments():
            i = 0
            while i < len(not_used_parameters):
                p = not_used_parameters[i]
                if re.search(r'\b{}\b'.format(re.escape(p.getId())), libsbml.formulaToString(action.getMath())):
                    not_used_parameters.remove(p)
                else:
                    i += 1
    #---
    for assignment in model.getListOfInitialAssignments():
        i = 0
        while i < len(not_used_parameters):
            p = not_used_parameters[i]
            if re.search(r'\b{}\b'.format(re.escape(p.getId())), libsbml.formulaToString(assignment.getMath())):
                not_used_parameters.remove(p)
            else:
                i += 1
    #---  
    for constraint in model.getListOfConstraints():
        i = 0
        while i < len(not_used_parameters):
            p = not_used_parameters[i]
            if re.search(r'\b{}\b'.format(re.escape(p.getId())), libsbml.formulaToString(constraint.getMath())):
                not_used_parameters.remove(p)
            else:
                i += 1
    #---
    for p in not_used_parameters:
        model.removeParameter(p.getId())
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
def get_Michaelis_Menten_reducted_model(model: libsbml.Model, verbose=0) -> bool:
    """Create the new Michaelis-Menten reactions and remove the original reactions from the model.
    
    Parameter
    ---------
    model : libsbml.Model
        The model from which the Michaelis-Menten schemes are extracted.
        
    Return
    ------
    bool
        True if the model has been modified, False otherwise.
        
    Output display
    --------------
    Michaelis-Menten pool corresponding to the enzyme D.
        Original Michaelis-Menten scheme: E + S1 + S2 + S3  <---> ES ---> S1 + S2 + S3
        Original reaction rates: ...; ...
        Reduced Michaelis-Menten scheme: S1 + S2 + S3 ---> S1 + S2 + S3
        Reduced reaction rate: ...
    """
    pools_of_schemes = get_Michaelis_Menten_pools(model, verbose)
    #---
    nb = 1
    nb_kcat = 1
    #---
    for enzyme, list_of_pools in pools_of_schemes.items():
        # For each pool of enzyme
        if verbose:
            print(f"{BOLD}Michaelis-Menten pool corresponding to the enzyme {GREEN}{enzyme}{RESET}.")
        #---
        # Create the liste of substrate for each pool of enzyme
        list_of_substrate = [list_of_schemes[0]['substrate'] for pool in list_of_pools for list_of_schemes in pool.values()]
        #---
        idt = 1
        for iiiii, pool_of_intermediate in enumerate(list_of_pools):
            # For each pool of intermediates
            #---
            nb_pool = 0
            for list_of_schemes in pool_of_intermediate.values():
                # For each pool of schemes
                #---
                for i, scheme in enumerate(list_of_schemes):
                    # For each scheme
                    #---
                    # We create the new Michaelis-Menten reaction
                    reaction = model.createReaction()
                    reaction.setId(f"MM_{enzyme}_{idt}")
                    idt += 1
                    reaction.setName(f"Michaelis-Menten-reduced conversion of {scheme['substrate']} to {' and '.join(scheme['product'])}")
                    reaction.setReversible(False)
                    reaction.setFast(False)
                    #---
                    # Adding the species to the reaction
                    for species in model.getListOfSpecies():
                        if species.getId() in scheme['product']:
                            reaction.addProduct(species)
                        elif species.getId() == scheme['substrate']:
                            reaction.addReactant(species)
                        elif species.getId() == enzyme:
                            reaction.addModifier(species)
                        elif species.getId() in [list_of_substrate[j] for j in range(len(list_of_pools))] and species.getId() != scheme['substrate']:
                            reaction.addModifier(species)
                    #---
                    # Adding the kinetic law to the reaction
                    kinetic_law = reaction.createKineticLaw()
                    #---
                    nom = f"{enzyme} * ({scheme['substrate']} / Km_{nb + iiiii})"
                    denom = f"1 + {' + '.join([f'{list_of_substrate[j]} / Km_{nb + j}' for j in range(len(list_of_pools))])}"
                    #---
                    intermediate_formula = f"{nom} / ({denom})"
                    #---
                    compartment = model.getSpecies(scheme['reaction_1'].getReactant(0).getSpecies()).getCompartment()
                    formula = f"{compartment} * Kcat_{nb_kcat + iiiii + nb_pool} * ({intermediate_formula})"
                    kinetic_law.setFormula(formula)
                    #---
                    # Adding the parameters to the reaction
                    p = kinetic_law.createParameter()
                    p.setId(f"Kcat_{nb_kcat + iiiii + nb_pool}")
                    p.setValue(scheme['k_cat'].getValue())
                    p.setUnits(scheme['k_cat'].getUnits())
                    #---
                    nb_pool += 1
                    #---
                    if verbose:
                        full_scheme = f"{BLUE}{scheme['substrate']}{RESET} + {GREEN}{enzyme}{RESET} <---> {ORANGE}{scheme['intermediate']}{RESET} ---> {GREEN}{enzyme}{RESET} + {RED}{' + '.join(scheme['product'])}{RESET}"
                        print(f"\t{BOLD}Original Michaelis-Menten scheme:{RESET} {full_scheme}")
                        print(f"\t{BOLD}Original reaction rates:{RESET} {MAGENTA}{scheme['reaction_1'].getKineticLaw().getFormula()}{RESET};  {MAGENTA}{scheme['reaction_2'].getKineticLaw().getFormula()}{RESET}")
                        new_scheme = f"{BLUE}{scheme['substrate']}{RESET} ---> {RED}{' + '.join(scheme['product'])}{RESET}"
                        print(f"\t{BOLD}Reduced Michaelis-Menten scheme:{RESET} {new_scheme}")
                        print(f"\t{BOLD}Reduced reaction rate:{RESET} {MAGENTA}{formula}{RESET}")
                        print()
                #---
                # Add the intermediate species to the assignment rule
                add_intermediates_to_rule(model, list_of_schemes[0]['intermediate'], intermediate_formula)
            #---
            # Add the Km parameters to the model
            p = model.createParameter()
            p.setId(f"Km_{nb + iiiii}")
            p.setValue(scheme['k_m_value'])
            p.setUnits(scheme['k_m_unit'])
            p.setConstant(True)
            #--- 
            nb_kcat += nb_pool - 1
        #---
        nb += len(list_of_pools)
        nb_kcat += len(list_of_pools)
    #---
    if pools_of_schemes:
        # Remove the original reactions from the model and the parameters that are not used now.
        delet_reactions(model, pools_of_schemes)
        remove_not_used_parameters(model)
        return True
    #---
    return False
#
#------------------------------------------------------------------------------------------------------------------------
#-----                                                 Main script                                                  -----
#------------------------------------------------------------------------------------------------------------------------
#
if __name__ == '__main__':
    # Default output directory
    output_dir = "michaelis_menten_reduce_models"  
    #---
    if(len(sys.argv) == 1):
        # If the program is start with no argument, open a pop up window to select files to be reduce.
        files_path = get_file()
        # If you want a pop up window for the output directory, uncomment the line below
        # output_dir = get_output_dir()
    else:
        # Else, check if the argument -o is present. If yes, take the next argument for the output directory
        # and the other files for biomodel to reduce, else take all the argument for biomodel. The output directory
        # will be the one by default.
        if "-o" in sys.argv:
            output_dir = sys.argv[sys.argv.index("-o")+1]
            if not os.path.exists(output_dir):
                # Create the directory if it doesn't exist.
                os.makedirs(output_dir)
            files_path = sys.argv[1:sys.argv.index("-o")] + sys.argv[sys.argv.index("-o")+2:]
        else:
            files_path = sys.argv[1:]          
    #---
    # For each files given as input
    for file_path in files_path:
        #---
        # Check if one is a directory, to take all the xml files in this directory.
        glob_parm = (file_path if ".xml" in file_path else (file_path + os.sep + f'*.xml'))
        for fname in glob.glob(glob_parm):
            # For each files, get the model, merge reaction and try the Michaelis-Menten Reduction
            # if there's a reduction, save the new model in the output directory
            #---
            file_name = os.path.splitext(os.path.basename(fname))[0]
            print(f"{BOLD}Reading the SBML file {BLUE}{file_name}.xml{RESET}...")
            #---
            start_time = time.time()
            #---
            # Read the sbml file and extract the model
            sbml_model, sbml_document = get_model(fname)
            #---
            # Modify reactions based on their directions
            merge_reaction(sbml_model, verbose = 0)
            #---
            # Get the new Michaelis-Menten reduced model
            mm_flag = get_Michaelis_Menten_reducted_model(sbml_model, verbose = 0)
            #---
            if mm_flag:
                # Save the new model in a new file of the same name
                save_model(sbml_document, file_name, output_dir)
                #---
                execution_time = time.time() - start_time
                #---
                print(f"{BOLD}Execution time is{RESET} {BLUE} {execution_time}{RESET} {BOLD}seconds.{RESET}\n")