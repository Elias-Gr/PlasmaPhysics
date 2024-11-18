import sys
import os
import subprocess

os.chdir('/Users/eliasgreil/visual studio/Plasma_Physics/')



subprocess.run(['python', 'main_optimization.py'])
subprocess.run(['python', 'run_poincare_simsopt.py'])