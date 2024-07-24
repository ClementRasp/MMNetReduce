#!/usr/bin/env python3
# coding: utf-8 
"""
This file contains all the funtions used by MMNetReduce.py and MMNetReduce_2.py.
"""
#------------------------------------------------------------------------------------------------------------------------
#-----                                                  Includes                                                    -----
#------------------------------------------------------------------------------------------------------------------------
import libsbml, os, sys
from tkinter import Tk, filedialog
#------------------------------------------------------------------------------------------------------------------------
#-----                                                  Constants                                                   -----
#------------------------------------------------------------------------------------------------------------------------
BOLD = "\033[1m"     #ANSI escape code for bold text
RESET = "\033[0m"    #Reset ANSI escape code
GREEN = "\033[32m"              # Enzymes for MM
BLUE = "\033[34m"               # Substrates for MM  &&  Highlighting
RED = "\033[31m"                # Products for MM
ORANGE = "\033[38;5;214m"       # Intermediate for MM
MAGENTA = "\033[35m"
#------------------------------------------------------------------------------------------------------------------------
#-----                                         File manipulation functions                                          -----
#------------------------------------------------------------------------------------------------------------------------
def get_file() -> tuple[str]:
    """Open a dialog to select the SBML file and return the file path and name.
    
    Returns
    -------
    str :
        The file path and the file name.
    """
    root = Tk()
    root.withdraw()  # Open file dialog to select the SBML file
    files_path = filedialog.askopenfilenames( 
                parent=root,
                title="Select SBML file",
                filetypes=[("xml", "*.xml"),("All files", "*.*")]
    )
    #---
    if not files_path: # If no file is selected
        print(f"\n{BOLD}No file selected.{RESET}\n")
        sys.exit(1)
    #---
    return files_path
#------------------------------------------------------------------------------------------------------------------------
def get_output_dir() -> str:
    """Open a dialog to select the output directory and return the directory path.
    
    Returns
    -------
    str :
        The directory path.
    """
    root = Tk() 
    root.withdraw()  
    # Open file dialog to select the output directory
    output_dir = filedialog.askdirectory( 
                parent=root,
                title="Select output directory"
    ) 
    #---
    if not output_dir: 
        # If no directory is selected
        output_dir = "MMNetReduce_output"  # Default output directory
    #---
    return output_dir

#------------------------------------------------------------------------------------------------------------------------
def get_model(file_path: str) -> tuple[libsbml.Model, libsbml.SBMLDocument]:
    """Read the SBML file and extract the model and document containing the model for the saving process.
    
    Parameter
    ---------
    file_path : str
        The path to the SBML file.
    
    Returns
    -------
    libsbml.Model
        The model extracted from the SBML file.
        
    libsbml.SBMLDocument   
        The document containing the model.
    """
    reader = libsbml.SBMLReader()
    document = reader.readSBMLFromFile(file_path)
    #---
    errors = document.getNumErrors()
    if errors > 0:
        print(f"{BOLD}Error loading the SBML file:{RESET}")
        for e in range(errors):
            print(f"\t{document.getError(e).getMessage()}")
        # sys.exit(1)
    #---
    document.enablePackage(libsbml.GroupsExtension.getXmlnsL3V1V1(), 'groups', False)
    model = document.getModel()
    #---
    if model is None:
        print(f"{BOLD}No model found in the SBML file.\n{RESET}")
        sys.exit(1)
    #---
    return model, document
#------------------------------------------------------------------------------------------------------------------------
def save_model(sbml_document: libsbml.SBMLDocument, file_name: str, output_dir) -> None:
    """Save the modified SBML file in the output directory.
    
    Parameters
    ----------
    sbml_document : libsbml.SBMLDocument
        The SBML document containing the modified model.
    
    file_name : str
        The name of the SBML file to be saved.
    """
    output_file = output_dir + os.sep + file_name + ".xml" # Creation of the output file
    #---
    try: # Create a directory to store the output file
        os.mkdir(output_dir)
    except OSError:
        pass
    #---
    status = libsbml.writeSBML(sbml_document, output_file)
    #---
    if status:
        print(f"The modified SBML file {BLUE}{file_name}.xml{RESET} is successfully saved in {BLUE}{output_dir}{RESET}.\n") 
    else:
        print(f"{BOLD}Error encountered while saving the modified SBML file.\n{RESET}")
#------------------------------------------------------------------------------------------------------------------------
