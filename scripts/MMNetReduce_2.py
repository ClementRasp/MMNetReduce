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
import libsbml, time, sympy as sp, glob
from useful import *
from merge_reaction import merge_reaction
#
#------------------------------------------------------------------------------------------------------------------------
#-----                                           Michaelis-Menten Schemes                                           -----
#------------------------------------------------------------------------------------------------------------------------
#
def get_Michaelis_Menten_scheme (model: libsbml.Model, verbose=0) -> list[dict]:
    """Get the Michaelis-Menten schemes from the model given as input.
    All the schemes are stored in a list of dictionaries. Each dictionary contains the reaction IDs forming the Michaelis-Menten scheme,
    the enzyme species, the intermediate enzyme species, and the lists of products and substrates involved in the Michaelis-Menten scheme.
    
    Parameter
    ---------
    model : libsbml.Model
        The model from which the Michaelis-Menten schemes are extracted.
        
    Return
    ------
    list[dict]
        The list of dictionaries representing the Michaelis-Menten schemes.
    
    Return exemple
    --------------
    [{reaction_1: R1, reaction_2: R2, enzyme: D, substrate: [A, B, C], product: [E, F], intermediate: D, id: 1},
     {reaction_1: R3, reaction_2: R4, enzyme: J, substrate: [G, H, I], product: [K, L], intermediate: J, id: 2}]
     
    Display exemple
    ---------------
    Michaelis-Menten scheme 1:
        Reversible reaction ID: R1
        Irreversible reaction ID: R2
        Scheme: A + B + C <---> D ---> E + F
        Enzyme: D
        Substrates: A, B, C
        Products: E, F
        Intermediate: D
        
    Michaelis-Menten scheme 2:
        Reversible reaction ID: R3
        Irreversible reaction ID: R4
        Scheme: G + H + I <---> J ---> K + L
        Enzyme: J
        Substrates: G, H, I
        Products: K, L
        Intermediate: J
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
            if len(d_2) != 1 or len(d_1 & d_3) != 1:
                # Skip if the reactions do not respect the Michaelis-Menten scheme
                # to keep only this type of reaction S + E <---> ES ---> E + P (for A or B = S + E, C = ES, D = E + P)
                continue
            #---
            if (len(d_1) < 2 and len(d_3) < 1) or (len(d_1) < 1 and len(d_3) < 2):
                # Skip if the reactions do not respect the Michaelis-Menten scheme
                # to keep only this type of reaction S + E <---> ES ---> E + P (for A or B = S + E, C = ES, D = E + P)
                continue
            #---
            scheme = dict() # Create a dictionary to store the Michaelis-Menten scheme
            enzyme = d_1 & d_3
            product = d_3 - enzyme
            substrate = d_1 - enzyme
            scheme['reaction_1'] = reaction_1
            scheme['reaction_2'] = reaction_2
            scheme['enzyme'] = enzyme.pop()
            scheme['substrate'] = list(substrate)
            scheme['product'] = list(product)
            scheme['intermediate'] = d_2.pop()
            scheme['id'] = id = id + 1
            michaelis_menten_scheme.append(scheme)
    #---
    if verbose:
        for scheme in michaelis_menten_scheme:
            print(f"{BOLD}Michaelis-Menten scheme {scheme['id']}:{RESET}")
            print(f"\t{BOLD}Reversible reaction ID:{RESET} {scheme['reaction_1'].getId()}")
            print(f"\t{BOLD}Irreversible reaction ID:{RESET} {scheme['reaction_2'].getId()}")
            full_scheme = f"{BLUE}" + " + ".join(scheme['substrate']) + f"{RESET} + {GREEN}{scheme['enzyme']}{RESET} <---> {ORANGE}{scheme['intermediate']}{RESET} ---> {GREEN}{scheme['enzyme']}{RESET} + {RED}" + " + ".join(scheme['product']) + f"{RESET}"
            print(f"\t{BOLD}Scheme:{RESET} {full_scheme}")
            print(f"\t{BOLD}Enzyme:{RESET} {GREEN}{scheme['enzyme']}{RESET}")
            print(f"\t{BOLD}Substrates:{RESET} {BLUE}{', '.join(scheme['substrate'])}{RESET}")
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
def get_Michaelis_Menten_pools(model: libsbml.Model, verbose=0) -> dict:
    """Get the Michaelis-Menten schemes stored by pool corresponding to the enzyme of the scheme.
    
    Parameter
    ---------
    model : libsbml.Model
        The model from which the Michaelis-Menten schemes are extracted.
        
    Return
    ------
    dict
        The dictionary of Michaelis-Menten schemes stored by pool corresponding to the enzyme of the scheme.
        
    Return exemple
    --------------
    {'D': [{reaction_1: R1, reaction_2: R2, enzyme: D, substrate: [A, B, C], product: [E, F], intermediate: DX, id: 1}],
     'J': [{reaction_1: R3, reaction_2: R4, enzyme: J, substrate: [G, H, I], product: [K, L], intermediate: JX_1, id: 2},
              {reaction_1: R5, reaction_2: R6, enzyme: J, substrate: [M, N, O], product: [P, Q], intermediate: JX_2, id: 3}]
     
    Display exemple
    ---------------
    Michaelis-Menten pool corresponding to the enzyme D.
        Michaelis-Menten scheme 1
        
    Michaelis-Menten pool corresponding to the enzyme J.
        Michaelis-Menten scheme 2\n
        Michaelis-Menten scheme 3
    """
    list_of_scheme = get_Michaelis_Menten_scheme(model, verbose)
    #---
    pools = dict()
    for scheme in list_of_scheme:
        if scheme['enzyme'] not in pools:
            pools[scheme['enzyme']] = [scheme]
        else:
            pools[scheme['enzyme']].append(scheme)
    #---  
    if verbose:
        for enzyme, schemes in pools.items():
            print(f"{BOLD}Michaelis-Menten pool corresponding to the enzyme {GREEN}{enzyme}{RESET}.")
            for scheme in schemes:
                print(f"\tMichaelis-Menten scheme {BOLD}{scheme['id']}{RESET}")
            print()
        print()
    #---
    return pools
#
#------------------------------------------------------------------------------------------------------------------------
#-----                                          Michaelis-Menten Reduction                                          -----
#------------------------------------------------------------------------------------------------------------------------
#
def string_to_expression(expression: str) -> sp.Expr:
    """ Convert a string expression to a sympy expression.
    
    Parameter
    ---------
    expression : str
        The string expression to convert.
        
    Return
    ------
    sp.Expr
        The sympy expression.
    """
    d = dict()
    for e in expression.replace('*', ' ').replace('/', ' ').replace('+', ' ').replace('-', ' ').replace('(',' ').replace(')',' ').split():
        if e != 'pow':
            d[e] = sp.Symbol(e)
    expression = sp.sympify(expression, d)
    return expression
#------------------------------------------------------------------------------------------------------------------------
def rename_parameter(reaction: str, dict_of_param: dict, suffix: str) -> tuple[sp.Expr, dict]:
    """ Rename the parameters of the reaction and return the new expression.
    
    Parameters
    ----------
    reaction : str
        The reaction to rename the parameters.
    
    dict_of_param : dict
        The dictionary of parameters to rename.
        
    suffix : str
        The suffix to add to the parameter name.
        
    Returns
    -------
    sp.Expr
        The new expression of the reaction.
        
    dict
        The dictionary of parameters with the new name.
    """
    expression = string_to_expression(reaction.getFormula())
    #---
    lst_of_param = [p.getId() for p in reaction.getListOfParameters()]
    all_param = set()
    #---
    for symbol in expression.free_symbols:
        if str(symbol) in lst_of_param:
            if symbol not in all_param:
                all_param.add(symbol)
                new_symbol = f"{symbol}_{suffix}"
                expression = expression.subs(symbol, sp.Symbol(new_symbol))
                p = reaction.getParameter(str(symbol))
                p.setId(new_symbol)
                dict_of_param[new_symbol] = (p.getValue(), p.getUnits())
    #---
    return expression, dict_of_param
#------------------------------------------------------------------------------------------------------------------------
def remove_Species_Reaction_Rules(model: libsbml.Model, pools_of_schemes: dict) -> None:
    """ Remove the species, reactions and rules that are not used in the model.
    
    Parameters
    ----------
    model : libsbml.Model
        The model from which the species, reactions and rules are removed.
        
    pools_of_schemes : dict
        The dictionary of Michaelis-Menten schemes stored by pool corresponding to the enzyme of the scheme.
    """  
    def is_species_used(model: libsbml.Model, species_id: str) -> bool:
        """ Check if the species is used in the model. """
        for reaction in model.getListOfReactions():
            for reactant in reaction.getListOfReactants():
                if reactant.getSpecies() == species_id:
                    return True
            for product in reaction.getListOfProducts():
                if product.getSpecies() == species_id:
                    return True
            for modifier in reaction.getListOfModifiers():
                if modifier.getSpecies() == species_id:
                    return True
                
        for rule in model.getListOfRules():
            if rule.getFormula().find(species_id) != -1:
                return True
            
        return False
    removed_species = set()
    #---
    for schemes in pools_of_schemes.values():
        for scheme in schemes:
            model.removeReaction(scheme['reaction_1'].getId())
            model.removeReaction(scheme['reaction_2'].getId())
            if not is_species_used(model, scheme['intermediate']):
                model.removeSpecies(scheme['intermediate'])  
                removed_species.add(scheme['intermediate'])
    #---
    plugin = model.getPlugin("groups")
    if plugin is not None:
        for group in plugin.getListOfGroups():
            member_list = [member for member in group.getListOfMembers()]
            for member in member_list:
                if member.getIdRef() in removed_species:
                    member.removeFromParentAndDelete()
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
def get_Michaelis_Menten_reducted_model(model: libsbml.Model, verbose=0) -> None:
    """Create the new Michaelis-Menten reactions and remove the original reactions from the model.
    
    Parameter
    ---------
    model : libsbml.Model
        The model from which the Michaelis-Menten schemes are extracted.
        
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
    for enzyme, list_of_schemes in pools_of_schemes.items():
        if verbose:
            print(f"{BOLD}Michaelis-Menten pool corresponding to the enzyme {GREEN}{enzyme}{RESET}.")
        #---
        list_of_expressions = []
        list_of_reaction_2_rate = []
        dict_of_param = dict()
        intermediates = [scheme['intermediate'] for scheme in list_of_schemes]
        #---
        for i, scheme in enumerate(list_of_schemes):
            reaction_reversible = scheme['reaction_1'].getKineticLaw()
            reaction_irreversible = scheme['reaction_2'].getKineticLaw()
            #---
            e1, dict_of_param = rename_parameter(reaction_reversible, dict_of_param, f"1_{i+1}")
            e2, dict_of_param = rename_parameter(reaction_irreversible, dict_of_param, f"2_{i+1}")
            #---            
            expression = e1 - e2
            expression = sp.simplify(expression)
            #---
            e_expr = f"{enzyme} - ({'+'.join(intermediates)})"
            e_expr = string_to_expression(e_expr)
            #---
            expression = expression.subs(enzyme, e_expr)
            list_of_expressions.append(expression)
            list_of_reaction_2_rate.append(e2)
        #---
        new_intermediates = sp.solve(list_of_expressions, [sp.Symbol(inter) for inter in intermediates])
        #---
        list_of_param = []
        for scheme in list_of_schemes:
            list_of_param += [p.getId() for p in scheme['reaction_1'].getKineticLaw().getListOfParameters()]
            list_of_param += [p.getId() for p in scheme['reaction_2'].getKineticLaw().getListOfParameters()]
        #---
        for i, scheme in enumerate(list_of_schemes):
            Michalelis_Menten_reaction = model.createReaction()
            Michalelis_Menten_reaction.setId(f"MM_{enzyme}_{i + 1}")
            Michalelis_Menten_reaction.setName(f"Michaelis-Menten reaction {i + 1} of {scheme['enzyme']} enzyme")
            Michalelis_Menten_reaction.setReversible(False)
            Michalelis_Menten_reaction.setFast(False)
            #---
            rate = list_of_reaction_2_rate[i]
            rate = rate.subs(new_intermediates)
            rate = sp.simplify(rate)
            #---
            k_law = Michalelis_Menten_reaction.createKineticLaw()
            k_law.setMath(libsbml.parseL3Formula(str(rate)))
            #---
            symbols = [str(symbol) for symbol in rate.free_symbols]
            #---
            for species in scheme['reaction_1'].getListOfReactants():
                if species.getSpecies() != enzyme:
                    Michalelis_Menten_reaction.addReactant(species)
                    if species.getSpecies() in symbols:
                        symbols.remove(species.getSpecies())
            for species in scheme['reaction_2'].getListOfProducts():
                if species.getSpecies() != enzyme:
                    Michalelis_Menten_reaction.addProduct(species)
            #---
            for p in list_of_param:
                if p in symbols:
                    param = k_law.createParameter()
                    param.setId(p)
                    param.setValue(dict_of_param[p][0])
                    for unit in model.getListOfUnitDefinitions():
                        if unit.getId() == dict_of_param[p][1]:
                            param.setUnits(unit.getId())
                    # param.setUnits(dict_of_param[p][1])
                    symbols.remove(p)
            #---
            for species in model.getListOfSpecies():
                if species.getId() in symbols:
                    Michalelis_Menten_reaction.addModifier(species)
            #---
            if verbose:
                full_scheme = f"{BLUE}{' + '.join(scheme['substrate'])}{RESET} + {GREEN}{enzyme}{RESET} <---> {ORANGE}{scheme['intermediate']}{RESET} ---> {GREEN}{enzyme}{RESET} + {RED}{' + '.join(scheme['product'])}{RESET}"
                print(f"\t{BOLD}Original Michaelis-Menten scheme:{RESET} {full_scheme}")
                print(f"\t{BOLD}Original reaction rates:{RESET} {MAGENTA}{scheme['reaction_1'].getKineticLaw().getFormula()}{RESET};  {MAGENTA}{scheme['reaction_2'].getKineticLaw().getFormula()}{RESET}")
                new_scheme = f"{BLUE}{' + '.join(scheme['substrate'])}{RESET} ---> {RED}{' + '.join(scheme['product'])}{RESET}"
                print(f"\t{BOLD}Reduced Michaelis-Menten scheme:{RESET} {new_scheme}")
                print(f"\t{BOLD}Reduced reaction rate:{RESET} {MAGENTA}{rate}{RESET}")
                print()
    #---
    remove_Species_Reaction_Rules(model, pools_of_schemes)
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