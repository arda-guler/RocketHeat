import analysis
import autoanalysis
import time
import sys
import os

from ui import clear_cmd_terminal

help_text = """
Thermal analysis tool for liquid propellant rocket engines.
Written by H. Arda Güler for METU Rocket Society in 2022.

Single Analysis: Performs analysis on an engine with predetermined properties.
Auto Analysis: Loops over some cooling system properties provided at initialization.

All rights reserved. Use of this software is only permitted to
members of METU Rocket Society and people given permission
explicitly.
"""

def is_run_from_idle():
    return bool("idlelib" in sys.modules)

def init():
    if not os.path.exists("inputs"):
        os.makedirs("inputs")

    if is_run_from_idle():
        print("Please do not run the program on IDLE. Use the terminal/cmd instead.")
        return
        
    clear_cmd_terminal()
    print(" = = = LPRE THERMAL ANALYSIS TOOL = = = ")
    print("")
    print("Please make a choice:")
    print(" 1) (S)ingle Analysis")
    print(" 2) (A)uto Analysis")
    print(" 3) (H)elp & About")
    print("")
    selection = input(" > ")
    
    if selection == "1" or selection.lower() == "s":
        analysis.perform()
    elif selection == "2" or selection.lower() == "a":
        autoanalysis.auto_analyze()
    elif selection == "3" or selection.lower() == "h":
        print(help_text)
        input("Press Enter to return.")
        init()
    else:
        print("Invalid selection!")
        time.sleep(2)
        init()

init()
